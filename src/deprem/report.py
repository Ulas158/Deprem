"""Eğitim hattı için rapor ve harita katmanları."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from .faults.budget import summarize_along_strike
from .faults.coulomb import coulomb_on_receivers
from .features.table import make_cells
from .region import TURKEY

NOTE = ("Deneysel araştırma çıktısıdır, deprem tahmini değildir. Model ve varsayımlar "
        "doğrulanmamıştır.")


def _fmt_table(rows):
    return pd.DataFrame(rows).to_markdown(index=False)


def write_report(run):
    ev = json.loads((run.dir / "evaluation.json").read_text())
    model = json.loads((run.dir / "model.json").read_text())
    etas = json.loads((run.dir / "etas_params.json").read_text())
    rows = ev["by_threshold"]
    md = [f"# Eğitim raporu: {run.cfg['run_name']}", "",
          f"> {NOTE} Sonuçlar, modelin eğitimde hiç görmediği {run.test_years[0]}–"
          f"{run.test_years[-1]} yıllarında, kabul görmüş referans model ETAS'a karşı ölçülmüştür.", "",
          "## Sade dille özet", ""]
    for r in rows:
        ig = r["bilgi_kazancı_olay_başına"]
        p = r["p_değeri"]
        if ig is None:
            verdict = "değerlendirilemedi (olay yok)"
        elif p is not None and p < 0.05 and ig > 0:
            verdict = "model ETAS'tan **istatistiksel olarak daha iyi**"
        elif p is not None and p < 0.05 and ig < 0:
            verdict = "model ETAS'tan **daha kötü**"
        else:
            verdict = "model ile ETAS arasında **anlamlı fark yok**"
        md.append(f"* **{r['eşik']}** ({r['gözlenen_olay']} olay): {verdict} "
                  f"(olay başına bilgi kazancı {ig if ig is not None else '–'}, p = "
                  f"{p:.3g})" if p is not None else f"* **{r['eşik']}**: {verdict}")
    md += ["", "Bilgi kazancı pozitifse model depremlerin olduğu hücrelere ETAS'tan daha yüksek "
           "olasılık vermiştir. p < 0,05 farkın şans eseri olma ihtimalinin düşük olduğunu "
           "gösterir; tek bir sınama dönemi olduğu için bu bile bir kanıt değil, ileriye dönük "
           "test edilmesi gereken bir işarettir.", "",
           "## Eşik eşik sonuçlar", "", _fmt_table(rows), "",
           f"## Yıl yıl (M≥{run.train_mag:g})", "", _fmt_table(ev["by_year"]), "",
           "## ETAS parametreleri (eğitim dönemi)", "",
           _fmt_table([{k: round(v, 4) for k, v in etas["params"].items()}]), "",
           f"Dallanma oranı (10 yıl kesilmiş): {etas['fit_info'].get('branching_ratio_truncated', float('nan')):.2f}, "
           f"b = {etas['fit_info']['b']:.2f}, hedef olay sayısı {etas['fit_info']['n_targets']}.", "",
           "## Öğrenilen model", "",
           f"Seçilen efektif sürtünme μ' = {model['mu_eff']} (doğrulama olabilirliği: " +
           ", ".join(f"{s['mu_eff']}: {s['val_loglik']:.1f}" for s in model["mu_eff_selection"]) + ").", "",
           "En büyük 15 katsayı (standartlaştırılmış özellik başına log-oran etkisi; pozitif = "
           "o özellik yüksekken ETAS'ın öngördüğünden fazla deprem):", ""]
    coef = pd.Series(model["coef"]).drop("(sabit)")
    top = coef.reindex(coef.abs().sort_values(ascending=False).index).head(15)
    md += [pd.DataFrame({"özellik": top.index, "katsayı": top.round(3).values}).to_markdown(index=False), "",
           "## Kalibrasyon (sınama dönemi)", "",
           "Her satır, modelin benzer olasılık verdiği hücre-yıllarını gruplar: iyi kalibre "
           "bir modelde ortalama olasılık ile gözlenen sıklık birbirine yakındır.", ""]
    for key, rel in ev["reliability"].items():
        df = pd.DataFrame(rel["model"])[["p_mean", "obs_freq", "n"]].round(3)
        md += [f"**{key} (model)**", "", df.to_markdown(index=False), ""]
    text = "\n".join(md)
    (run.dir / "report.md").write_text(text)
    rep = run.root / "reports" / f"training_{run.cfg['run_name']}.md"
    rep.write_text(text)


def _layer(dir_: Path, id_, label, desc, group, order, geojson, color, legend, tooltip,
           style=None, note=NOTE):
    dir_.mkdir(parents=True, exist_ok=True)
    (dir_ / f"{id_}.geojson").write_text(json.dumps(geojson, ensure_ascii=False, separators=(",", ":")))
    spec = {"id": id_, "label": label, "description": desc, "group": group, "order": order,
            "data": f"{id_}.geojson", "visible": False, "note": note,
            "updated": pd.Timestamp.now(tz="UTC").strftime("%Y-%m-%d"),
            "style": style or {"type": "auto", "opacity": 0.8, "width": 2.5},
            "color": color, "legend": legend, "tooltip": tooltip}
    (dir_ / f"{id_}.layer.json").write_text(json.dumps(spec, ensure_ascii=False, indent=2))


def _clean(v):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return None
    if isinstance(v, (np.floating, float)):
        return round(float(v), 4)
    if isinstance(v, (np.integer,)):
        return int(v)
    return v


def _lines(df, props):
    return {"type": "FeatureCollection", "metadata": {"note": NOTE},
            "features": [{"type": "Feature", "geometry": {"type": "LineString", "coordinates": r["coords"]},
                          "properties": {p: _clean(r[p]) for p in props}} for _, r in df.iterrows()]}


def write_layers(run):
    out = run.root / run.cfg["map"]["layer_dir"]
    el = run.elements()
    # 1) kayma açığı (eleman bütçesi)
    fin = pd.read_csv(run.dir / "element_budget_final.csv.gz")
    s = summarize_along_strike(el, fin[["slip_deficit_m", "last_release_year"]], bin_km=3.0)
    s = s.merge(el.groupby("fault_id")["slip_rate_mm_yr"].first(), on="fault_id")
    _layer(out, "slip_deficit", "Kayma açığı (1 km eleman)",
           "Her 1 km fay elemanının 1900'den bu yana biriktirip depremle boşaltmadığı kayma (m), "
           "3 km parçalara ve derinliğe göre ortalanmış.", "Faylar ve enerji", 10,
           _lines(s, ["fault_id", "slip_deficit_m", "last_release_year", "slip_rate_mm_yr"]),
           {"property": "slip_deficit_m", "type": "interpolate",
            "stops": [[0, "#2a78d6"], [0.5, "#8a5cf0"], [1.0, "#e0497a"], [2.0, "#ff8c42"], [3.0, "#ffe680"]]},
           {"title": "Kayma açığı", "unit": "m", "labels": ["0", "1", "2", "3+"]},
           {"title": "fault_id", "fields": [
               {"property": "slip_deficit_m", "label": "Kayma açığı", "unit": "m", "digits": 2},
               {"property": "last_release_year", "label": "Son boşalım", "digits": 0},
               {"property": "slip_rate_mm_yr", "label": "Kayma hızı", "unit": "mm/yıl", "digits": 1}]},
           style={"type": "line", "width": 3, "glow": True, "opacity": 0.9})

    # 2) Coulomb gerilme değişimi (bugün, seçilen μ')
    model = json.loads((run.dir / "model.json").read_text())
    src = pd.read_csv(run.dir / "coulomb_sources.csv")
    rec = pd.read_csv(run.dir / "coulomb_receivers.csv.gz")
    per = np.load(run.dir / "coulomb_stress_per_source.npy")
    dcff, _, _ = coulomb_on_receivers(per.sum(0).astype(float), rec["strike"].to_numpy(),
                                      rec["dip"].to_numpy(), rec["rake"].fillna(180.0).to_numpy(),
                                      model["mu_eff"])
    rec["dcff_bar"] = np.clip(dcff / 1e5, -50, 50)
    sc = summarize_along_strike(rec.assign(j_dip=0), rec[["dcff_bar"]], bin_km=3.0)
    _layer(out, "coulomb_dcff", "Coulomb gerilme değişimi",
           f"{len(src)} büyük depremin (M≥6, 1900+) fay elemanlarında oluşturduğu birikimli ΔCFF "
           f"(bar, μ'={model['mu_eff']}). Pozitif: kırılmaya yaklaşmış; negatif: gerilme gölgesi.",
           "Gerilme ve ETAS", 20, _lines(sc, ["fault_id", "dcff_bar"]),
           {"property": "dcff_bar", "type": "interpolate",
            "stops": [[-1, "#3f7be0"], [-0.1, "#7fa6e8"], [0, "#5b6170"], [0.1, "#ff9a8f"], [1, "#ff3b30"]]},
           {"title": "ΔCFF", "unit": "bar", "labels": ["≤-1", "0", "≥+1"]},
           {"title": "fault_id", "fields": [{"property": "dcff_bar", "label": "ΔCFF", "unit": "bar", "digits": 2}]},
           style={"type": "line", "width": 3, "glow": True, "opacity": 0.9})

    # 3) ETAS ve model: son sınama yılı, hücre poligonları
    pred = pd.read_csv(run.dir / "predictions.csv.gz")
    yr = run.test_years[-1]
    d = 0.5 * run.cfg["data"]["cell_deg"]
    b = model["etas_b"]
    sc4 = 10 ** (-b * (4.0 - run.train_mag))
    p = pred[pred["year"] == yr]
    feats = []
    for _, r in p.iterrows():
        lat, lon = r["lat_c"], r["lon_c"]
        ring = [[round(lon - d, 4), round(lat - d, 4)], [round(lon + d, 4), round(lat - d, 4)],
                [round(lon + d, 4), round(lat + d, 4)], [round(lon - d, 4), round(lat + d, 4)],
                [round(lon - d, 4), round(lat - d, 4)]]
        feats.append({"type": "Feature", "geometry": {"type": "Polygon", "coordinates": [ring]},
                      "properties": {"etas_m4": _clean(r["lam_etas"] * sc4),
                                     "model_m4": _clean(r["lam_model"] * sc4),
                                     "log10_ratio": _clean(np.log10(r["lam_model"] / r["lam_etas"])),
                                     "observed_m4": int(r.get("y_count_m4", 0))}})
    gj = {"type": "FeatureCollection", "metadata": {"note": NOTE, "year": yr}, "features": feats}
    tip = {"title": "Hücre", "fields": [
        {"property": "etas_m4", "label": "ETAS beklenen M≥4", "digits": 2},
        {"property": "model_m4", "label": "Model beklenen M≥4", "digits": 2},
        {"property": "observed_m4", "label": f"{yr} gözlenen M≥4", "digits": 0},
        {"property": "log10_ratio", "label": "log10(model/ETAS)", "digits": 2}]}
    _layer(out, "etas_rate", f"ETAS beklenen M≥4 ({yr})",
           f"{yr} yılı için ETAS'ın 1 Ocak itibarıyla beklediği M≥4 sayısı (0,5° hücre). "
           "Sınama dönemi: model bu yılı eğitimde görmedi.", "Gerilme ve ETAS", 21, gj,
           {"property": "etas_m4", "type": "interpolate",
            "stops": [[0, "#0b1220"], [0.1, "#1f3b73"], [0.5, "#3f7be0"], [2, "#8a5cf0"], [10, "#ffe680"]]},
           {"title": "Beklenen M≥4", "unit": "olay/yıl", "labels": ["0", "0.5", "2", "10+"]}, tip,
           style={"type": "fill", "opacity": 0.45})
    _layer(out, "model_vs_etas", f"Model / ETAS oranı ({yr})",
           "Eğitilen modelin ETAS'a göre oranı (log10). Kırmızı: model ETAS'tan fazla bekliyor.",
           "Gerilme ve ETAS", 22, gj,
           {"property": "log10_ratio", "type": "interpolate",
            "stops": [[-0.5, "#3f7be0"], [0, "#2a2f3a"], [0.5, "#ff3b30"]]},
           {"title": "log10(model/ETAS)", "unit": "", "labels": ["-0.5", "0", "+0.5"]}, tip,
           style={"type": "fill", "opacity": 0.45})

    # 4) eleman parametreleri (seçili, 3 km özet)
    from .faults.parameters import element_seismicity
    seis = element_seismicity(el, run.catalog, run.catalog["time"].max(), mc=run.cfg["data"]["mc"],
                              radius_km=5.0, independent=run.catalog["independent"])
    cols = ["seis_n_10y", "seis_b_10y", "seis_rate_z", "seis_max_mag", "seis_mean_depth"]
    vals = pd.concat([el[["slip_rate_mm_yr", "depth_rel", "end_dist_km", "strike_bend_deg"]], seis[cols]], axis=1)
    sp = summarize_along_strike(el, vals, bin_km=3.0)
    props = ["fault_id"] + list(vals.columns)
    _layer(out, "element_params", "Fay elemanı parametreleri",
           "1 km elemanların parametreleri (3 km parçalara özetlenmiş): kayma hızı, son 10 yıl "
           "olay sayısı, b-değeri, hız değişimi (z), en büyük büyüklük, ortalama derinlik.",
           "Faylar ve enerji", 11, _lines(sp, props),
           {"property": "seis_n_10y", "type": "interpolate",
            "stops": [[0, "#2a2f3a"], [5, "#3f7be0"], [50, "#8a5cf0"], [500, "#ff8c42"]]},
           {"title": "Son 10 yıl olay (5 km içinde)", "unit": "", "labels": ["0", "5", "50", "500+"]},
           {"title": "fault_id", "fields": [
               {"property": "slip_rate_mm_yr", "label": "Kayma hızı", "unit": "mm/yıl", "digits": 1},
               {"property": "seis_n_10y", "label": "Son 10 yıl olay", "digits": 0},
               {"property": "seis_b_10y", "label": "b-değeri", "digits": 2},
               {"property": "seis_rate_z", "label": "Hız değişimi z", "digits": 1},
               {"property": "seis_max_mag", "label": "En büyük M", "digits": 1},
               {"property": "seis_mean_depth", "label": "Ort. derinlik", "unit": "km", "digits": 1}]},
           style={"type": "line", "width": 2.5, "opacity": 0.9})

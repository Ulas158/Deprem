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
    md += _obs_sections(run, model)
    text = "\n".join(md)
    (run.dir / "report.md").write_text(text)
    rep = run.root / "reports" / f"training_{run.cfg['run_name']}.md"
    rep.write_text(text)


def _verdict(ig, p):
    if ig is None or (isinstance(ig, float) and not np.isfinite(ig)):
        return "değerlendirilemedi"
    if p is not None and p < 0.05:
        return "anlamlı iyileşme" if ig > 0 else "anlamlı kötüleşme"
    return "anlamlı fark yok"


def _obs_sections(run, model):
    path = run.dir / "obs_evaluation.json"
    if not path.exists():
        return []
    from .observations.registry import table as obs_table
    ev = json.loads(path.read_text())
    md = ["", "## Yeni gözlemler: jeodezik, sismik ve hipotezler", "",
          "Her gözlem için literatürdeki fiziksel ilişkiden bir **başlangıç etkisi** (öncül) ve "
          "belirsizliği tanımlandı; model bunları 2012-2019 depremleriyle eğitip güncelledi. "
          "Sınama yine modelin görmediği yıllarda, ETAS'a karşıdır.", ""]
    reg = obs_table()
    md += ["### Özellik listesi ve durumları", "",
           reg[["name", "group", "evidence", "horizon", "status", "prior_mean", "prior_sd", "source"]]
           .rename(columns={"name": "özellik", "group": "grup", "evidence": "kanıt", "horizon": "vade",
                            "status": "durum", "prior_mean": "başlangıç etkisi", "prior_sd": "belirsizlik",
                            "source": "kaynak"}).to_markdown(index=False), ""]
    cov = ev.get("coverage", {}).get("overall")
    if cov:
        md += ["Kapsama (hücre-yıllarının kaçında değer var):", "",
               pd.DataFrame({"özellik": list(cov), "kapsama": [f"%{100 * v:.0f}" for v in cov.values()]})
               .to_markdown(index=False), ""]
    pp = model.get("obs_prior_posterior", [])
    if pp:
        df = pd.DataFrame(pp)
        df["oran_etkisi (×)"] = np.exp(df["öğrenilen"]).round(2)
        md += ["### Başlangıç etkisi → eğitimle öğrenilen etki (yıllık model)", "",
               "Katsayı: özellik 1 standart sapma arttığında deprem oranının logaritmasındaki değişim. "
               "Sonsal sd, öncül sd'den çok küçükse veri bu etki hakkında bilgi taşıyor demektir; "
               "yakınsa sonuç büyük ölçüde başlangıç varsayımıdır.", "",
               df.round(3).to_markdown(index=False), ""]
    ab = ev.get("yearly_ablation")
    if ab:
        md += ["### Gözlemlerin katkısı (sınama yılları, ETAS'a karşı olay başına bilgi kazancı)", "",
               "Aynı model gözlemli ve gözlemsiz eğitildi; fark, gözlemlerin katkısıdır.", "",
               pd.DataFrame(ab).to_markdown(index=False, floatfmt=".4g"), ""]
        base = ab[0]
        full = ab[1]
        for key in [k for k in full if k.endswith("kazanç")]:
            b, f = base.get(key), full.get(key)
            if b is not None and f is not None:
                md.append(f"* **{key.replace(' kazanç', '')}**: gözlemsiz {b:+.3f} → gözlemli {f:+.3f} "
                          f"({_verdict(f, full.get(key.replace('kazanç', 'p')))} ETAS'a göre)")
        md.append("")
    stv = ev.get("shortterm")
    if stv:
        md += ["### Zayıf kanıtlı hipotezlerin kısa vadeli testi (iyonosfer, ULF, termal)", "",
               "2° × 2° kutular, 5 günlük pencereler. Soru: son 10 günde anomali varsa, sonraki 5 günde "
               "bağımsız (artçı olmayan) deprem olasılığı, depremlerin kendi kümelenmesinden "
               "beklenenden yüksek mi? Eğitim 2012-2019, sınama 2020-2025. Alarm olasılık kazancı: "
               "eğitimdeki %90'lık eşiği aşan anomaliden sonraki pencerelerde deprem payı / alarm "
               "altındaki uzay-zaman payı (1 = şans düzeyi).", ""]
        for mag, r in stv.items():
            md += [f"**Hedef {mag}**", "", pd.DataFrame(r["results"]).to_markdown(index=False, floatfmt=".4g"), ""]
            if r.get("coefs"):
                md += ["Başlangıç etkisi → öğrenilen:", "",
                       pd.DataFrame(r["coefs"]).to_markdown(index=False, floatfmt=".3g"), ""]
    return md


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
    if run.obs_enabled():
        write_obs_layers(run, out)


def _points(df, props, lat="lat", lon="lon"):
    return {"type": "FeatureCollection", "metadata": {"note": NOTE},
            "features": [{"type": "Feature",
                          "geometry": {"type": "Point", "coordinates": [round(float(r[lon]), 4), round(float(r[lat]), 4)]},
                          "properties": {p: _clean(r[p]) for p in props}} for _, r in df.iterrows()]}


def write_obs_layers(run, out):
    """Yeni gözlem katmanları: GPS istasyonları ve gerinim hızı, tekrarlayan depremler, dv/v."""
    from .observations import gnss, repeaters
    obs = run.root / "data" / "obs"
    group = "Jeodezi ve sismik gözlemler"
    if (obs / "gnss_station_year.csv").exists():
        sty = pd.read_csv(obs / "gnss_station_year.csv")
        last = sty.sort_values("year").groupby("site").tail(1).copy()
        last["v_horiz"] = np.hypot(last["ve"], last["vn"])
        _layer(out, "gnss_stations", "GPS istasyonları (hız, geçici hareket)",
               "Nevada Geodetic Lab günlük çözümleri, Avrasya sabit. Renk: son 1 yılda uzun dönem "
               "hızdan sapmanın gürültüye oranı (z). Her istasyonun verisi olan son yıl.", group, 30,
               _points(last, ["site", "year", "ve", "vn", "v_horiz", "tr_z", "tr_e", "tr_n"]),
               {"property": "tr_z", "type": "interpolate",
                "stops": [[0, "#2a78d6"], [2, "#8a5cf0"], [4, "#ff8c42"], [8, "#ffe680"]]},
               {"title": "Geçici hareket z", "unit": "", "labels": ["0", "2", "4", "8+"]},
               {"title": "site", "fields": [
                   {"property": "year", "label": "Yıl başı", "digits": 0},
                   {"property": "v_horiz", "label": "Yatay hız", "unit": "mm/yıl", "digits": 1},
                   {"property": "tr_z", "label": "Geçici z", "digits": 1},
                   {"property": "tr_e", "label": "Hız sapması D", "unit": "mm/yıl", "digits": 1},
                   {"property": "tr_n", "label": "Hız sapması K", "unit": "mm/yıl", "digits": 1}]},
               style={"type": "circle", "radius": 4.5, "outline": "#ffffff", "opacity": 0.9})
        cells = make_cells(TURKEY, 0.25)
        sr = gnss.strain_rate_at(cells["lat_c"], cells["lon_c"], last)
        cells = pd.concat([cells, sr], axis=1).dropna(subset=["strain_2inv"])
        d = 0.125
        feats = [{"type": "Feature", "geometry": {"type": "Polygon", "coordinates": [[
            [round(r.lon_c - d, 3), round(r.lat_c - d, 3)], [round(r.lon_c + d, 3), round(r.lat_c - d, 3)],
            [round(r.lon_c + d, 3), round(r.lat_c + d, 3)], [round(r.lon_c - d, 3), round(r.lat_c + d, 3)],
            [round(r.lon_c - d, 3), round(r.lat_c - d, 3)]]]},
            "properties": {"strain_ns": _clean(r.strain_2inv), "shear_ns": _clean(r.max_shear),
                           "dilat_ns": _clean(r.dilatation), "n_sta": int(r.n_sta)}} for r in cells.itertuples()]
        _layer(out, "gnss_strain", "GPS gerinim hızı",
               "İstasyon hızlarının Gauss ağırlıklı (40 km) gradyanından ikinci değişmez gerinim hızı "
               "(nanogerinim/yıl). Kabukta enerjinin ne hızla biriktiğini gösterir.", group, 31,
               {"type": "FeatureCollection", "metadata": {"note": NOTE}, "features": feats},
               {"property": "strain_ns", "type": "interpolate",
                "stops": [[0, "#0b1220"], [30, "#1f3b73"], [100, "#3f7be0"], [300, "#e0497a"], [1000, "#ffe680"]]},
               {"title": "Gerinim hızı", "unit": "nanogerinim/yıl", "labels": ["0", "100", "300", "1000+"]},
               {"title": "Hücre", "fields": [
                   {"property": "strain_ns", "label": "Gerinim hızı", "unit": "ns/yıl", "digits": 0},
                   {"property": "shear_ns", "label": "Makaslama", "unit": "ns/yıl", "digits": 0},
                   {"property": "dilat_ns", "label": "Hacim değişimi", "unit": "ns/yıl", "digits": 0},
                   {"property": "n_sta", "label": "İstasyon sayısı", "digits": 0}]},
               style={"type": "fill", "opacity": 0.45})
    fam = repeaters.find_families(run.catalog)
    if not fam.empty:
        g = fam.groupby("family").agg(lat=("latitude", "mean"), lon=("longitude", "mean"),
                                      depth_km=("depth_km", "mean"), n=("mag", "size"), mag=("mag", "mean"),
                                      first=("time", "min"), last=("time", "max"), slip_cm=("slip_cm", "sum"))
        span = ((g["last"] - g["first"]).dt.days / 365.25).clip(lower=0.5)
        g["slip_rate_cm_yr"] = g["slip_cm"] / span
        g["first"], g["last"] = g["first"].dt.year, g["last"].dt.year
        _layer(out, "repeaters", "Tekrarlayan deprem adayları",
               "Aynı noktada (≤1 km), benzer büyüklükte ve düzenli aralıklarla tekrarlayan küçük depremler. "
               "Fayın derindeki sessiz kayma hızını gösterir (Nadeau & Johnson 1998). Katalog konum hatası "
               "nedeniyle adaydır; dalga biçimi ile doğrulanmamıştır.", group, 32,
               _points(g.reset_index(), ["family", "n", "mag", "depth_km", "first", "last", "slip_rate_cm_yr"]),
               {"property": "slip_rate_cm_yr", "type": "interpolate",
                "stops": [[0, "#2a78d6"], [0.5, "#8a5cf0"], [2, "#ff8c42"], [5, "#ffe680"]]},
               {"title": "Kayma hızı", "unit": "cm/yıl", "labels": ["0", "0.5", "2", "5+"]},
               {"title": "family", "fields": [
                   {"property": "n", "label": "Olay sayısı", "digits": 0},
                   {"property": "mag", "label": "Ort. büyüklük", "digits": 1},
                   {"property": "depth_km", "label": "Derinlik", "unit": "km", "digits": 1},
                   {"property": "first", "label": "İlk", "digits": 0},
                   {"property": "last", "label": "Son", "digits": 0},
                   {"property": "slip_rate_cm_yr", "label": "Kayma hızı", "unit": "cm/yıl", "digits": 2}]},
               style={"type": "circle", "radius": 4, "outline": "#ffffff", "opacity": 0.9})
    if (obs / "dvv_series.csv.gz").exists():
        dv = pd.read_csv(obs / "dvv_series.csv.gz", parse_dates=["time"])
        rows = []
        for s, g in dv.groupby("sta"):
            g = g.dropna(subset=["dvv"])
            if len(g) < 20:
                continue
            t1 = g["time"].max()
            rec = g[g["time"] >= t1 - pd.DateOffset(months=6)]["dvv"].mean()
            prev = g[(g["time"] < t1 - pd.DateOffset(months=6)) & (g["time"] >= t1 - pd.DateOffset(months=30))]["dvv"].mean()
            rows.append({"sta": s, "lat": g["lat"].iloc[0], "lon": g["lon"].iloc[0],
                         "dvv_change_pct": 100 * (rec - prev), "n": len(g), "last": t1.strftime("%Y-%m-%d")})
        if rows:
            _layer(out, "dvv_stations", "Sismik hız değişimi (dv/v)",
                   "Ortam gürültüsünden tek istasyon hız değişimi: son 6 ay − önceki 2 yıl (%). "
                   "Negatif: kayaç hızı düşmüş (hasar veya akışkan basıncı artışı).", group, 33,
                   _points(pd.DataFrame(rows), ["sta", "dvv_change_pct", "n", "last"]),
                   {"property": "dvv_change_pct", "type": "interpolate",
                    "stops": [[-0.1, "#ff3b30"], [0, "#2a2f3a"], [0.1, "#3f7be0"]]},
                   {"title": "dv/v değişimi", "unit": "%", "labels": ["-0.1", "0", "+0.1"]},
                   {"title": "sta", "fields": [
                       {"property": "dvv_change_pct", "label": "dv/v değişimi", "unit": "%", "digits": 3},
                       {"property": "last", "label": "Son ölçüm"}]},
                   style={"type": "circle", "radius": 6, "outline": "#ffffff", "opacity": 0.9})

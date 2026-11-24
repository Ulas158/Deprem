"""Türkiye sismisite + fay moment bütçesi analizini çalıştır ve rapor üret.

Önce: python scripts/fetch_catalogs.py
Sonra: python scripts/run_analysis.py

Çıktılar:
  reports/figures/*.png
  reports/tables/*.csv
  reports/turkey_summary.md
  data/processed/feature_table.csv.gz  (öğrenme için hücre × yıl tablosu)

Bu betiğin ürettiği hiçbir sayı deprem tahmini değildir.
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from deprem import plotting as P
from deprem.analysis.seismicity import b_value_map, b_value_timeseries, gridded_counts
from deprem.catalog.clean import decluster_gardner_knopoff
from deprem.faults.gem import download as download_gem, load_faults
import json

from deprem.faults.mesh import mesh_faults
from deprem.faults.moment import add_geometry, fault_summary, zone_budget
from deprem.faults.parameters import element_seismicity, ensure_columns, registry_table
from deprem.region import SUBREGIONS as _SUB
from deprem.features.table import build_table
from deprem.forecast.base import ForecastGrid, SmoothedSeismicity, UniformPoisson
from deprem.forecast.tests import information_gain, l_test, n_test
from deprem.region import SUBREGIONS, TURKEY
from deprem.stats.gr import fit_gr, mc_gft, mc_maxc

ROOT = Path(__file__).resolve().parents[1]
UTC = "UTC"


def ts(s):
    return pd.Timestamp(s, tz=UTC)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modern-start", default="2010-01-01",
                    help="Küçük depremler için homojen dönem başlangıcı")
    ap.add_argument("--end", default=None)
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

    fig_dir = ROOT / "reports" / "figures"
    tab_dir = ROOT / "reports" / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir.mkdir(parents=True, exist_ok=True)

    cat = pd.read_csv(ROOT / "data/processed/turkey_catalog.csv.gz")
    cat["time"] = pd.to_datetime(cat["time"], utc=True, format="ISO8601")
    end = ts(args.end) if args.end else cat["time"].max().ceil("D")
    modern_start = ts(args.modern_start)
    modern = cat[(cat["time"] >= modern_start) & (cat["time"] < end)]
    years_modern = (end - modern_start).days / 365.25
    coast = P.load_coastlines(ROOT / "data/ext/ne_50m_coastline.geojson", TURKEY)
    md = []  # rapor satırları

    # ---------------------------------------------------------------- katalog özeti
    by_src = cat["source"].value_counts().to_dict()
    md += ["# Türkiye sismisite ve fay moment bütçesi: ilk prototip çıktıları", "",
           "> **Uyarı.** Bu rapordaki hiçbir sayı bir deprem tahmini değildir. Hepsi "
           "betimsel istatistik veya ileride ileriye dönük olarak sınanacak hipotez "
           "değişkenleridir.", "",
           f"Üretim zamanı: {pd.Timestamp.now(tz=UTC):%Y-%m-%d %H:%M} UTC. "
           f"Katalog bitişi: {end:%Y-%m-%d}.", "",
           "## 1. Katalog", "",
           f"* Birleşik olay sayısı: **{len(cat):,}** "
           f"({', '.join(f'{k}: {v:,}' for k, v in by_src.items())}; tekilleştirme sonrası "
           "tutulan kaydın kaynağı)",
           f"* Zaman aralığı: {cat['time'].min():%Y-%m-%d} – {cat['time'].max():%Y-%m-%d}",
           f"* Birden fazla kaynakta eşleşen olay: {(cat['n_sources'] > 1).sum():,}",
           f"* Mw bilgisi olan olay: {cat['mw'].notna().sum():,}", ""]

    # ---------------------------------------------------------------- GR
    mc_max = mc_maxc(modern["mag"])
    mc_g, r_g = mc_gft(modern["mag"], level=90)
    mc = float(np.nanmax([mc_max, mc_g]))
    fit = fit_gr(modern["mag"], years_modern, mc=mc, n_boot=200)
    P.fmd_plot(modern["mag"].to_numpy(), fit, fig_dir / "fmd_turkey.png",
               title=f"Türkiye frekans-büyüklük dağılımı ({modern_start:%Y}–{end:%Y})")
    md += ["## 2. Gutenberg-Richter", "",
           f"Dönem {modern_start:%Y}–{end:%Y} ({years_modern:.1f} yıl), Türkiye kutusu.", "",
           f"* Mc (maks. eğrilik +0,2): {mc_max:.1f}; Mc (iyi-uyum, R={r_g:.0f}%): {mc_g:.1f}; "
           f"kullanılan: **{mc:.1f}**",
           f"* b = **{fit.b:.3f}** ± {fit.sigma_b:.3f} (Shi-Bolt), önyükleme σ = {fit.b_boot_std:.3f}, "
           f"N(M≥Mc) = {fit.n:,}",
           f"* Yıllık a-değeri = {fit.a_annual:.2f}. GR ekstrapolasyonuyla yılda M≥5: "
           f"{fit.annual_rate(5):.1f}, M≥6: {fit.annual_rate(6):.2f}, M≥7: {fit.annual_rate(7):.3f} "
           "(artçılar dahil, betimsel)", "",
           "![FMD](figures/fmd_turkey.png)", ""]

    rows = []
    for key, reg in SUBREGIONS.items():
        sub = modern[reg.contains(modern)]
        if len(sub) < 100:
            continue
        # MAXC karışık/ardışık-dizi kataloglarında Mc'yi düşük tahmin eder; GFT ile üstten sınırla
        mc_s = float(np.nanmax([mc_maxc(sub["mag"]), mc_gft(sub["mag"], level=90)[0]]))
        f_s = fit_gr(sub["mag"], years_modern, mc=mc_s, n_boot=100)
        rows.append({"bölge": reg.name, "n": len(sub), "Mc": round(mc_s, 1),
                     "b": round(f_s.b, 3), "σ_b": round(f_s.sigma_b, 3),
                     "yıllık_M≥5_GR": round(f_s.annual_rate(5), 2)})
    sub_df = pd.DataFrame(rows)
    sub_df.to_csv(tab_dir / "subregion_gr.csv", index=False)
    md += ["Alt bölgeler:", "", sub_df.to_markdown(index=False), "",
           "Not: Doğu Anadolu Fay Zonu'ndaki düşük b büyük ölçüde 2023 Kahramanmaraş "
           "dizisinden gelir; büyük bir ana şoktan hemen sonra küçük artçılar kaydedilemez "
           "(kısa süreli tamlık kaybı), bu da b'yi yapay olarak düşürür. AFAD verisi ISC "
           "üzerinden ~2025 sonuna kadar gelir; 2026 için yalnızca EMSC ve AFAD'ın son iki "
           "haftası vardır, bu yüzden 2026 küçük olay sayıları eksiktir.", ""]

    # ---------------------------------------------------------------- zaman
    P.magnitude_time(cat[cat["time"] >= ts("1990-01-01")], fig_dir / "magnitude_time.png",
                     mc=mc, title="Büyüklük-zaman ve aylık olay sayısı")
    ind = decluster_gardner_knopoff(modern)
    mods = modern[modern["mag"] >= mc]
    indep = modern[ind & (modern["mag"] >= mc)]
    P.cumulative_plot({
        "Tüm olaylar (M≥Mc)": mods.set_index("time").assign(c=1)["c"].cumsum(),
        "Artçılar ayıklanmış (Gardner-Knopoff)": indep.set_index("time").assign(c=1)["c"].cumsum(),
    }, fig_dir / "cumulative.png", title="Birikimli olay sayısı")
    bt = b_value_timeseries(modern, mc, window=500, step=100)
    P.b_time_plot(bt, fig_dir / "b_time.png", title="b-değeri zaman serisi (500 olayluk kayan pencere)")
    md += ["## 3. Zaman içinde sismisite", "",
           f"* Gardner-Knopoff ayıklaması sonrası bağımsız olay oranı (M≥Mc): "
           f"{len(indep) / max(len(mods), 1):.0%}",
           f"* b(t) aralığı: {bt['b'].min():.2f} – {bt['b'].max():.2f}. Kısa süreli düşüşler "
           "çoğunlukla büyük deprem dizilerine (ör. 2011 Van, 2020 Elazığ/Sisam, 2023 "
           "Kahramanmaraş) ve ağ değişimlerine denk gelir; tek başına öncül kabul edilmemeli.", "",
           "![Büyüklük-zaman](figures/magnitude_time.png)", "",
           "![Birikimli](figures/cumulative.png)", "", "![b(t)](figures/b_time.png)", ""]

    # ---------------------------------------------------------------- faylar
    gem_path = download_gem(ROOT / "data/ext/gem_active_faults_harmonized.geojson")
    faults = add_geometry(load_faults(gem_path, TURKEY))
    faults = faults[faults["slip_rate_mm_yr"].notna()].reset_index(drop=True)
    P.fault_rate_map(faults, TURKEY, fig_dir / "fault_slip_rates.png", coast,
                     title=f"Aktif faylar ve kayma hızları ({len(faults)} segment, GEM/EMME)")
    hist_start = ts("1900-01-01")
    big = cat[cat["mw_best"] >= 5.0]
    zb = zone_budget(faults, big, {**SUBREGIONS, "turkey": TURKEY}, hist_start, end)
    zb.to_csv(tab_dir / "zone_moment_budget.csv", index=False)
    P.moment_budget_plot(zb, fig_dir / "moment_budget.png",
                         title=f"Moment birikimi ve boşalımı, {hist_start:%Y}–{end:%Y}")
    fs = fault_summary(faults, cat, hist_start, end, min_mw=6.0)
    fs.sort_values("m0_rate_Nm_yr", ascending=False).to_csv(
        tab_dir / "fault_segments.csv", index=False)
    zshow = zb[["zone", "n_faults", "m0_accumulated", "m0_released", "release_ratio",
                "largest_event_mw", "largest_event_year"]].copy()
    for c in ("m0_accumulated", "m0_released"):
        zshow[c] = zshow[c].map(lambda v: f"{v:.2e}")
    zshow["release_ratio"] = zshow["release_ratio"].map(lambda v: f"{v:.2f}")
    top = fs.sort_values("mw_equiv_since_ref", ascending=False).head(15)
    tshow = top[["fault_id", "slip_type", "length_km", "slip_rate_mm_yr", "last_big_event",
                 "years_since_ref", "mw_equiv_since_ref"]].copy()
    tshow["last_big_event"] = tshow["last_big_event"].dt.strftime("%Y-%m-%d").fillna("katalogda yok")
    tshow = tshow.round({"length_km": 0, "slip_rate_mm_yr": 1, "years_since_ref": 0,
                         "mw_equiv_since_ref": 2})
    md += ["## 4. Faylar ve moment bütçesi", "",
           f"GEM Küresel Aktif Fay veritabanından (EMME modeli) Türkiye kutusuna değen "
           f"{len(faults)} segment. Her segment için M0̇ = μ·A·v (μ = 30 GPa, A = uzunluk × "
           "sismojenik genişlik, v = kayma hızı).", "",
           "![Kayma hızları](figures/fault_slip_rates.png)", "",
           f"Zon başına {hist_start:%Y}–{end:%Y} birikimi ile katalogdaki M≥5 olayların "
           "momenti (bağlaşım = 1 varsayımı):", "", zshow.to_markdown(index=False), "",
           "![Moment bütçesi](figures/moment_budget.png)", "",
           "**Okuma notu.** `release_ratio` 1'den küçükse zon katalog süresince biriktirdiğinden "
           "az moment boşaltmış görünür; bu bir 'gecikmiş deprem' kanıtı değildir. 1900 öncesi "
           "depremler (ör. 1509, 1766, 1894 Marmara; 1822, 1872 Doğu Anadolu) bilinmez, "
           "kayma hızları 2–4 kat belirsizdir ve aseismik krip (ör. İsmetpaşa) bağlaşımı "
           "düşürür.", "",
           "Katalog içi biriken momenti en yüksek 15 segment (son M≥6'dan veya 1900'den bu yana; "
           "segment düzeyinde atama kabadır):", "", tshow.to_markdown(index=False), ""]

    # ---------------------------------------------------------------- mekân
    h, la, lo = gridded_counts(modern[modern["mag"] >= mc], TURKEY, 0.2)
    P.density_map(h, la, lo, TURKEY, fig_dir / "density_map.png", coast, faults,
                  title=f"Olay yoğunluğu, M≥{mc:.1f}, {modern_start:%Y}–{end:%Y} (0,2° hücre)",
                  label="Olay sayısı")
    P.epicenter_map(cat[cat["mag"] >= 4.0], TURKEY, fig_dir / "epicenters_m4.png", coast, faults,
                    title="M≥4 episantrlar (tüm katalog)")
    bmap, mcmap, la, lo = b_value_map(modern, TURKEY, cell_deg=0.25, n_nearest=300,
                                      max_radius_km=50, min_events=100)
    P.b_map(bmap, la, lo, TURKEY, fig_dir / "b_value_map.png", coast, faults,
            title="b-değeri haritası (en yakın 300 olay, r ≤ 50 km, yerel Mc)")
    md += ["## 5. Mekânda sismisite", "",
           "![Yoğunluk](figures/density_map.png)", "", "![M4 episantrlar](figures/epicenters_m4.png)", "",
           "![b haritası](figures/b_value_map.png)", "",
           "Not: b haritası artçı ayıklanmamış katalogla yapıldı; ayıklama küçük olayları "
           "orantısız çıkardığı için b'yi yapay olarak düşürüyor. Bazı yüksek-b kümeleri taş "
           "ocağı/maden patlamalarından kaynaklanıyor olabilir (EMSC kataloğunda olay türü "
           "yok; AFAD/ISC verisi patlamaları ayrı işaretliyor ve onları zaten dışarıda bırakıyoruz).", "",
           "Düşük b-değerli bölgeler bazı çalışmalarda yüksek diferansiyel gerilmeyle "
           "ilişkilendirilmiştir (Schorlemmer vd. 2005). Bu bir hipotezdir; harita bu "
           "hipotezi sınamak için girdi olarak tutulur.", ""]

    # ---------------------------------------------------------------- referans modeller
    fit_start, fit_end, test_end = ts("2010-01-01"), ts("2020-01-01"), ts("2025-01-01")
    grid = ForecastGrid(TURKEY, 0.5, np.array([4.0, 4.5, 5.0, 5.5, 6.0, 10.0]))
    models = {m.name: m.fit(cat, fit_start, fit_end) for m in
              (UniformPoisson(4.0, fit.b), SmoothedSeismicity(4.0, fit.b, bandwidth_km=30))}
    obs = grid.count(cat[(cat["time"] >= fit_end) & (cat["time"] < test_end) & (cat["mag"] >= 4.0)])
    exp = {k: m.expected(grid, fit_end, test_end) for k, m in models.items()}
    trow = []
    for k, e in exp.items():
        n = n_test(e, obs)
        lt = l_test(e, obs, n_sim=300)
        trow.append({"model": k, "N_beklenen": round(n.detail["n_forecast"], 1),
                     "N_gözlenen": n.detail["n_observed"], "N_testi_q": round(n.quantile, 4),
                     "logL": round(lt.statistic, 1), "L_testi_γ": round(lt.quantile, 3)})
    ig = information_gain(exp["smoothed_seismicity"], exp["uniform_poisson"], obs)
    md += ["## 6. Referans modeller ve test altyapısı", "",
           f"Eğitim {fit_start:%Y}–{fit_end:%Y}, sınama {fit_end:%Y}–{test_end:%Y}, M≥4, 0,5° hücre. "
           "Bu iki model ileride ETAS ve fay/özellik tabanlı modellerin geçmesi gereken "
           "asgari çıtadır.", "", pd.DataFrame(trow).to_markdown(index=False), "",
           f"Düzleştirilmiş sismisitenin düzgün modele göre olay başına bilgi kazancı: "
           f"**{ig.statistic:.2f}** nat (p = {ig.quantile:.2g}). Her iki model de zamanda "
           "sabit olduğundan 2023 Kahramanmaraş dizisini (kümelenme) yakalayamaz; N-testindeki "
           "düşük q bu yüzdendir ve ETAS'ın neden gerekli olduğunu gösterir.", ""]

    # ---------------------------------------------------------------- fay elemanları
    refine = [(_SUB["marmara"], 0.5)]
    elements = mesh_faults(faults, ds=1.0, refine=refine)
    cat_ind = decluster_gardner_knopoff(cat)
    seis = element_seismicity(elements, cat, end, mc=mc, radius_km=5.0, independent=cat_ind)
    elements = ensure_columns(pd.concat([elements, seis], axis=1))
    elements.to_csv(ROOT / "data/processed/fault_elements.csv.gz", index=False)
    reg = registry_table()
    reg.to_csv(tab_dir / "element_parameters.csv", index=False)
    filled = {n: float(elements[n].notna().mean()) for n in reg["name"]}
    reg["dolu_oran"] = reg["name"].map(filled).map(lambda v: f"{v:.0%}")
    md += ["## 7. Fay kabuk elemanları", "",
           f"{len(faults)} segment 1 km × 1 km elemanlara bölündü (Marmara'da 0,5 km): "
           f"**{len(elements):,} eleman**, toplam fay alanı {elements['area_km2'].sum():,.0f} km². "
           "Her elemanın parametre satırı `data/processed/fault_elements.csv.gz` dosyasında. "
           "Sismisite parametreleri elemana 5 km (3B) içindeki olaylardan, katalog sonuna kadar "
           "hesaplandı.", "",
           f"Kayıtlı {len(reg)} parametre (henüz verisi olmayanlar boş alan olarak tanımlı):", "",
           reg[["name", "unit", "status", "dolu_oran", "description"]].to_markdown(index=False), ""]

    # harita için segment düzeyinde moment açığı katmanı
    feats = []
    for _, r in fs.merge(faults[["fault_id", "coords"]], on="fault_id").iterrows():
        props = {k: (None if pd.isna(r[k]) else (float(r[k]) if isinstance(r[k], (int, float, np.floating)) else str(r[k])))
                 for k in ("fault_id", "slip_type", "slip_rate_mm_yr", "length_km", "m0_rate_Nm_yr",
                           "years_since_ref", "m0_since_ref", "mw_equiv_since_ref")}
        props["mw_equiv"] = props["mw_equiv_since_ref"]  # web haritasının beklediği alan
        props["last_big_event"] = None if pd.isna(r["last_big_event"]) else f"{r['last_big_event']:%Y-%m-%d}"
        feats.append({"type": "Feature", "properties": props,
                      "geometry": {"type": "LineString", "coordinates": [list(c) for c in r["coords"]]}})
    web = ROOT / "web" / "data"
    web.mkdir(parents=True, exist_ok=True)
    (web / "moment_deficit.geojson").write_text(json.dumps({
        "type": "FeatureCollection",
        "metadata": {"generated": f"{pd.Timestamp.now(tz=UTC):%Y-%m-%dT%H:%MZ}",
                     "reference_start": f"{hist_start:%Y-%m-%d}", "catalog_end": f"{end:%Y-%m-%d}",
                     "note": "Katalog içi biriken moment; tahmin değildir. 1900 öncesi depremler hesaba katılmamıştır."},
        "features": feats}, ensure_ascii=False))

    # ---------------------------------------------------------------- özellik tablosu
    cat = cat.assign(independent=cat_ind)
    starts = [ts(f"{y}-01-01") for y in range(2010, end.year)]
    table = build_table(cat, TURKEY, 0.5, starts, pd.DateOffset(years=1),
                        target_mags=[3.0, 4.0, 5.0, 6.0], mc=mc, faults=faults)
    table.to_csv(ROOT / "data/processed/feature_table.csv.gz", index=False)
    feat_cols = [c for c in table.columns if not c.startswith("y_")
                 and c not in ("window_start", "lat_c", "lon_c", "cell_id")]
    pos = {m: int(table[f"y_any_m{m:g}"].sum()) for m in (3, 4, 5, 6)}
    md += ["## 8. Öğrenme için özellik tablosu", "",
           f"`data/processed/feature_table.csv.gz`: {table['cell_id'].nunique()} hücre × "
           f"{len(starts)} yıl = {len(table):,} satır, {len(feat_cols)} özellik "
           f"({', '.join(feat_cols)}).", "",
           "Pozitif (en az bir olay olan) hücre-yıl sayısı: "
           + ", ".join(f"M≥{m}: {v:,}" for m, v in pos.items()) + ".", "",
           "Özellikler yalnızca pencere başlangıcından önceki veriden hesaplanır. Yeni "
           "değişkenler `deprem/features/table.py` içindeki `FEATURES` listesine eklenir.", ""]

    (ROOT / "reports" / "turkey_summary.md").write_text("\n".join(md))
    logging.info("Rapor yazıldı: reports/turkey_summary.md")


if __name__ == "__main__":
    main()

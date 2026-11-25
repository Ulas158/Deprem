"""Gaz/sıvı ve dış etkiler: katalog üzerinde hipotez testleri + harita katmanları.

Önce: python scripts/fetch_catalogs.py
Sonra: python scripts/run_external.py

Çıktılar:
  reports/external_summary.md              bölgesel Schuster testleri, sürüler, uzak tetikleme
  reports/tables/external_*.csv
  web/data/layers/ext_*.layer.json + .geojson   haritada aç/kapa katmanlar
  data/processed/fault_elements_external.csv.gz (fay elemanları varsa)

Bu betiğin ürettiği hiçbir sayı deprem tahmini değildir.
"""

from __future__ import annotations

import argparse
import io
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from deprem.catalog.clean import decluster_gardner_knopoff
from deprem.external.dynamic import remote_triggering
from deprem.external.layers import (GROUP_FLUID, GROUP_LOADS, cells_geojson, points_geojson,
                                    site_layer_specs, write_layer)
from deprem.external.params import DATA_SOURCES
from deprem.external.periodic import event_phases, schuster
from deprem.external.sites import load_sites
from deprem.external.swarms import cluster_stats, link_clusters
from deprem.features.table import assign_cell, make_cells
from deprem.region import SUBREGIONS, TURKEY
from deprem.stats.gr import mc_maxc

ROOT = Path(__file__).resolve().parents[1]
LAYERS = ROOT / "web" / "data" / "layers"
USGS = "https://earthquake.usgs.gov/fdsnws/event/1/query"


def fetch_remote(path: Path, start="1990-01-01", min_mag=6.5) -> pd.DataFrame:
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        r = requests.get(USGS, params={"format": "csv", "starttime": start, "minmagnitude": min_mag,
                                       "orderby": "time-asc", "limit": 20000}, timeout=180)
        r.raise_for_status()
        path.write_text(r.text)
    df = pd.read_csv(path)
    df["time"] = pd.to_datetime(df["time"], utc=True, format="ISO8601")
    return df[["time", "latitude", "longitude", "depth", "mag", "place"]]


def fmt_p(lp):
    return "—" if not np.isfinite(lp) else f"10^{lp:.1f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modern-start", default="2008-01-01")
    ap.add_argument("--cell", type=float, default=0.5)
    ap.add_argument("--no-tidal", action="store_true", help="Yarı günlük gelgit fazını atla (hızlı)")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    log = logging.getLogger("external")

    tab = ROOT / "reports" / "tables"
    tab.mkdir(parents=True, exist_ok=True)
    cat = pd.read_csv(ROOT / "data/processed/turkey_catalog.csv.gz")
    cat["time"] = pd.to_datetime(cat["time"], utc=True, format="ISO8601")
    cat = cat[TURKEY.contains(cat)].reset_index(drop=True)
    modern = cat[cat["time"] >= pd.Timestamp(args.modern_start, tz="UTC")].copy()
    mc = float(np.ceil(mc_maxc(modern["mag"]) * 10) / 10)
    modern = modern[modern["mag"] >= mc].reset_index(drop=True)
    log.info("Mc=%.1f, %d olay", mc, len(modern))
    modern["independent"] = decluster_gardner_knopoff(modern).to_numpy()
    ind = modern[modern["independent"]].copy()
    log.info("bağımsız olay: %d; fazlar hesaplanıyor", len(ind))
    ph = event_phases(ind, tidal=not args.no_tidal)
    ind = ind.join(ph)

    md = ["# Gaz, sıvı ve dış etkiler: ilk sınama sonuçları", "",
          "> **Uyarı.** Bu rapordaki hiçbir sayı bir deprem tahmini değildir. Hepsi geçmiş veride "
          "sınanan hipotezlerdir. Küçük p değerleri tek başına nedensellik göstermez; çok sayıda "
          "bölge ve faz denendiği için yanlış pozitif beklenir (Bonferroni: 0,05 / test sayısı).", "",
          f"Katalog: {args.modern_start[:4]}–{modern['time'].max():%Y-%m-%d}, Mc = {mc:.1f}, "
          f"{len(modern):,} olay, artçılar ayıklandıktan sonra (Gardner-Knopoff) {len(ind):,} bağımsız olay.", "",
          "Bilimsel arka plan, mekanizmalar ve kanıt düzeyleri: "
          "[docs/akiskan_ve_dis_etkiler.md](../docs/akiskan_ve_dis_etkiler.md).", ""]

    # ------------------------------------------------------------ 1. periyodik yükler
    PH = [("phase_hour", "saat"), ("phase_lunar", "ay_12sa"), ("phase_semidiurnal", "gelgit_birlesik"),
          ("phase_fortnight", "gelgit_14g"), ("phase_annual", "mevsim")]
    rows = []
    regions = {"turkiye": TURKEY, **SUBREGIONS}
    for key, reg in regions.items():
        sub = ind[reg.contains(ind)]
        for mmin in (mc, 3.0, 4.0):
            s = sub[sub["mag"] >= mmin]
            r = {"bölge": reg.name, "M≥": mmin, "N": len(s)}
            for col, name in PH:
                lp, mean, n = schuster(s[col])
                r[f"{name}_log10p"] = round(lp, 2) if np.isfinite(lp) else np.nan
                r[f"{name}_faz"] = round(mean, 0) if np.isfinite(mean) else np.nan
            rows.append(r)
    per = pd.DataFrame(rows)
    per.to_csv(tab / "external_periodic.csv", index=False)
    phys = ["ay_12sa_log10p", "gelgit_14g_log10p", "mevsim_log10p"]
    n_tests = per[phys].notna().sum().sum()
    thr = np.log10(0.05 / max(n_tests, 1))
    tr = per[per["bölge"] == TURKEY.name].set_index("M≥")
    md += ["## 1. Gelgit ve mevsim (Schuster testi)", "",
           "Her olay için şu fazlar hesaplandı: **saat** (yerel güneş saati; fiziksel değil, katalog "
           "kalitesi kontrolü), **ay_12sa** (yalnızca Ay'ın yarı günlük gelgiti; asıl gelgit testi), "
           "**gelgit_birlesik** (Ay + Güneş gelgit eğrisinde faz, Tanaka 2002), **gelgit_14g** "
           "(dolunay/yeni ay 0°) ve **mevsim** (1 Ocak 0°). Fazlar rastgeleyse p büyüktür.", "",
           "### Önemli bulgu: sahte gelgit sinyali", "",
           f"Küçük depremlerde (M≥{mc:.1f}) günün saatine çok güçlü bir yığılma var "
           f"(log10 p = {tr.loc[mc, 'saat_log10p']:.0f}); M≥3'te bu kayboluyor "
           f"(log10 p = {tr.loc[3.0, 'saat_log10p']:.1f}). Bu, gece/gündüz algılama farkından ve/veya "
           "katalogda kalan taş ocağı patlamalarından kaynaklanır, depremlerin fiziğinden değil. Güneş "
           "gelgiti tam 12 saatlik olduğundan aynı etki Ay+Güneş birleşik gelgit fazında da "
           f"görünüyor (M≥{mc:.1f}: log10 p = {tr.loc[mc, 'gelgit_birlesik_log10p']:.0f}). Güneş'ten "
           "bağımsız olan **Ay gelgiti testinde ise sinyal yok** "
           f"(M≥{mc:.1f}: {tr.loc[mc, 'ay_12sa_log10p']:.2f}, M≥3: {tr.loc[3.0, 'ay_12sa_log10p']:.2f}, "
           f"M≥4: {tr.loc[4.0, 'ay_12sa_log10p']:.2f}). Sonuç: Türkiye kataloğunda gelgit etkisine dair "
           "kanıt bulunamadı; birleşik gelgit fazı küçük depremlerde kullanılmamalı. Modelde gelgit "
           "özelliği Ay fazından, saat fazı ise ayrı bir kalite kontrol özelliği olarak hesaplanıyor.", "",
           f"Fiziksel testler (ay_12sa, gelgit_14g, mevsim) toplam {n_tests} kez yapıldığı için "
           f"anlamlılık eşiği log10 p < {thr:.1f} alınmalıdır (Bonferroni).", "",
           per.to_markdown(index=False), ""]
    sig = per[(per[phys] < thr).any(axis=1)]
    md += [f"Fiziksel testlerde eşiği geçen satır sayısı: **{len(sig)}**"
           + (" (" + ", ".join(f"{r['bölge']} M≥{r['M≥']:g}" for _, r in sig.iterrows()) + "). "
              "Bunlar doğrulanmış değildir; sonraki dönem verisinde tekrar etmesi ve katalog tamlığının "
              "mevsimsel değişimiyle açıklanamaması gerekir." if len(sig)
              else ". Hiçbir bölgede gelgit veya mevsim etkisi istatistiksel olarak anlamlı değil."), ""]

    # yıl yıl tutarlılık: gerçek bir mevsim etkisi her yıl benzer fazda görünmeli
    s3 = ind[ind["mag"] >= 3.0]
    yr = []
    for y, g in s3.groupby(s3["time"].dt.year):
        lp, mean, n = schuster(g["phase_annual"])
        yr.append({"yıl": int(y), "N": n, "mevsim_log10p": round(lp, 1), "faz": round(mean, 0)})
    yr = pd.DataFrame(yr)
    lp_all = schuster(s3["phase_annual"])[0]
    lp_ex = schuster(s3[s3["time"].dt.year != 2023]["phase_annual"])[0]
    n_strong = int((yr["mevsim_log10p"] < -2).sum())
    md += ["### Mevsim sinyali yıl yıl", "",
           f"Türkiye M≥3 bağımsız olaylarda mevsim testi log10 p = {lp_all:.1f}; 2023 çıkarılınca {lp_ex:.1f}. "
           f"Ancak sinyal yıllara eşit dağılmıyor: {len(yr)} yılın yalnızca {n_strong} tanesinde log10 p < −2. "
           "Bu, gerçek bir mevsimsel yükleme etkisinden çok birkaç yıldaki büyük deprem dizilerinin "
           "artçı ayıklamadan kaçan kısmından ve kataloğun yıllar içinde değişen büyüklük ölçeğinden "
           "(M≥3 olay sayısı 2008'de ~1270, sonraki yıllarda ~300) kaynaklanıyor olabilir. "
           "Şimdilik **kanıt sayılmaz**; daha homojen bir katalog ve ETAS tabanlı ayıklamayla yeniden sınanmalı.", "",
           yr.to_markdown(index=False), ""]

    # ------------------------------------------------------------ 2. sürüler
    log.info("kümeleme")
    cid = link_clusters(modern)
    cs = cluster_stats(modern, cid)
    cs.to_csv(tab / "external_clusters.csv", index=False)
    sw = cs[cs["is_swarm"]]
    md += ["## 2. Deprem sürüleri ve akışkan göçü", "",
           f"Uzay-zaman bağlantısıyla (10 km, 5 gün) en az 5 olaylı **{len(cs):,}** küme bulundu; "
           f"bunların **{len(sw):,}** tanesi sürü ölçütünü sağlıyor (en büyük olay başta değil ve ikinci "
           "en büyükten en fazla 0,5 büyük). Sürülerde olayların zamanla yayılmasından hidrolik yayınım "
           "D hesaplandı (r² = 4πDt). Akışkan kaynaklı sürülerde D tipik olarak 0,01–10 m²/s'dir.", ""]
    if len(sw):
        d = sw["D_m2_s"].dropna()
        md += [f"* D hesaplanabilen sürü: {len(d)}; ortanca D = {d.median():.2g} m²/s "
               f"(%25–%75: {d.quantile(.25):.2g}–{d.quantile(.75):.2g}). "
               "Konum hataları (birkaç km) küçük sürülerde D'yi büyütür; değerler üst sınır gibi okunmalı.",
               "", "En büyük 15 sürü:", "",
               sw.sort_values("n", ascending=False).head(15)[
                   ["t_start", "n", "duration_days", "lat", "lon", "depth_km", "mag_max", "D_m2_s", "radius_km"]
               ].assign(t_start=lambda d: d["t_start"].dt.strftime("%Y-%m-%d")).round(3).to_markdown(index=False), ""]

    # ------------------------------------------------------------ 3. uzak tetikleme
    remote = fetch_remote(ROOT / "data/raw/usgs_global_m65.csv")
    cells = make_cells(TURKEY, args.cell)
    modern["cell_id"] = assign_cell(modern, TURKEY, cells)
    ind["cell_id"] = assign_cell(ind, TURKEY, cells)
    lat_c = cells.set_index("cell_id")["lat_c"]
    lon_c = cells.set_index("cell_id")["lon_c"]
    rem = remote[remote["time"] >= modern["time"].min()]
    # artçı dizileri önceki 30 günü şişirip β'yı yapay olarak düşürür: bağımsız olaylar kullanılır
    dyn = remote_triggering(ind[ind["cell_id"] >= 0], rem, lat_c, lon_c, "cell_id")
    tot = remote_triggering(ind.assign(g=0), rem, pd.Series([39.0], index=[0]),
                            pd.Series([35.0], index=[0]), "g")
    md += ["## 3. Uzak büyük depremlerin dinamik tetiklemesi", "",
           f"{len(rem)} küresel M≥6,5 depremden, Türkiye'ye 1000 km'den uzak ve genlik vekili yeterli "
           f"olanlar ({int(tot['dyn_n_remote'].iloc[0])} olay) için, her birinden sonraki 5 gündeki yerel "
           "bağımsız (artçısı ayıklanmış) olay sayısı önceki 30 günle karşılaştırıldı "
           "(β istatistiği; β > 2 kabaca anlamlı artış, β < −2 azalma).", "",
           f"* Türkiye geneli: β = **{tot['dyn_beta'].iloc[0]:.2f}** "
           f"(5 gün içinde {int(tot['dyn_n_after'].iloc[0])} olay, toplam {int(tot['dyn_n_total'].iloc[0])})",
           f"* β > 2 olan hücre sayısı (≥20 olaylı hücreler): "
           f"{int(((dyn['dyn_beta'] > 2) & (dyn['dyn_n_total'] >= 20)).sum())} / {int((dyn['dyn_n_total'] >= 20).sum())}"
           " (çoklu test nedeniyle birkaç tanesi tesadüfen beklenir)", ""]

    # ------------------------------------------------------------ 4. hücre katmanları
    grid = cells.set_index("cell_id")[["lat_c", "lon_c"]].copy()
    for col, name in (("phase_lunar", "tidal_log10p"), ("phase_fortnight", "fortnight_log10p"),
                      ("phase_annual", "seasonal_log10p"), ("phase_hour", "hour_log10p")):
        g = ind[ind["cell_id"] >= 0].groupby("cell_id")[col].agg(lambda s: schuster(s.to_numpy())[0])
        grid[name] = g
        grid[name.replace("log10p", "neglog10p")] = -g
    grid["n_independent"] = ind[ind["cell_id"] >= 0].groupby("cell_id").size()
    grid = grid.join(dyn)
    grid = grid.reset_index()
    grid.to_csv(tab / "external_cells.csv", index=False)

    def cell_layer(lid, label, desc, prop, title, stops, extra_fields, sel, note, order):
        g = grid[sel]
        fields = [{"property": prop, "label": title, "digits": 2}] + extra_fields
        write_layer(LAYERS, {
            "id": lid, "label": label, "description": desc, "group": GROUP_LOADS, "order": order,
            "style": {"type": "fill", "opacity": 0.55},
            "color": {"property": prop, "type": "interpolate", "stops": stops, "missing": "rgba(0,0,0,0)"},
            "legend": {"title": title},
            "tooltip": {"title": "cell_id", "fields": fields}, "note": note,
        }, cells_geojson(g, args.cell, ["cell_id", prop, "n_independent", "dyn_n_total", "dyn_n_remote",
                                        "tidal_log10p", "fortnight_log10p", "seasonal_log10p", "dyn_beta"], note))

    nf = [{"property": "n_independent", "label": "Bağımsız olay", "digits": 0}]
    blues = [[0, "#1b2433"], [1.3, "#3f7be0"], [2, "#8a5cf0"], [3, "#e0497a"], [5, "#ffe680"]]
    cell_layer("ext_gelgit", "Gelgit duyarlılığı (Ay, 12 sa)", "Schuster −log10 p, 0,5° hücre",
               "tidal_neglog10p", "−log10 p (gelgit)", blues, nf, grid["n_independent"] >= 10,
               "Yüksek değer olayların gelgit fazına yığıldığını gösterir; çoklu test nedeniyle tek başına anlamlı değildir.", 500)
    cell_layer("ext_gelgit14", "Gelgit duyarlılığı (14 gün)", "Kuvvetli/zayıf gelgit fazı",
               "fortnight_neglog10p", "−log10 p (14 g)", blues, nf, grid["n_independent"] >= 10,
               "Dolunay/yeni ay dönemine yığılma testi; tahmin değildir.", 510)
    cell_layer("ext_mevsim", "Mevsimsel duyarlılık", "Yıllık faz Schuster −log10 p",
               "seasonal_neglog10p", "−log10 p (mevsim)", blues, nf, grid["n_independent"] >= 10,
               "Mevsimsel yük (yağış, kar, yeraltı suyu) vekili; katalog tamlığının mevsimsel değişimi de aynı izi bırakabilir.", 520)
    cell_layer("ext_saat", "Gün içi saat yığılması (kontrol)", "Gece/gündüz algılama farkı veya patlamalar",
               "hour_neglog10p", "−log10 p (saat)", blues, nf, grid["n_independent"] >= 10,
               "Fiziksel bir etki değil, katalog kalitesi göstergesidir: yüksek değer taş ocağı patlamaları "
               "veya gündüz gürültüsü olabilir. Bu hücrelerde gelgit/mevsim sonuçlarına güvenmeyin.", 540)
    cell_layer("ext_uzak_tetikleme", "Uzak deprem tetiklemesi", "Küresel M≥6,5 sonrası 5 gün β",
               "dyn_beta", "β (uzak tetikleme)",
               [[-3, "#3f7be0"], [0, "#1b2433"], [2, "#e0497a"], [4, "#ffe680"]],
               [{"property": "dyn_n_total", "label": "Pencerelerdeki olay", "digits": 0},
                {"property": "dyn_n_remote", "label": "Kullanılan uzak deprem", "digits": 0}],
               grid["dyn_n_total"] >= 20,
               "β > 2 kabaca anlamlı artış; jeotermal ve akışkanca zengin bölgelerde beklenir.", 530)

    # ------------------------------------------------------------ 5. sürü katmanı
    swl = sw.assign(year=sw["t_start"].dt.year, log10_D=np.log10(sw["D_m2_s"]))
    write_layer(LAYERS, {
        "id": "ext_suruler", "label": "Deprem sürüleri (akışkan göçü)", "group": GROUP_FLUID, "order": 410,
        "description": "Katalogdan bulunan sürüler; renk hidrolik yayınım",
        "style": {"type": "circle", "outline": "#0b111a", "opacity": 0.85},
        "color": {"property": "log10_D", "type": "interpolate", "missing": "#8796ab",
                  "stops": [[-2, "#3f7be0"], [-1, "#8a5cf0"], [0, "#e0497a"], [1, "#ff8c42"], [2, "#ffe680"]]},
        "size": {"property": "n", "stops": [[5, 4], [50, 9], [500, 16]]},
        "legend": {"title": "log10 D", "unit": "m²/s"},
        "tooltip": {"title": "t_start", "fields": [
            {"property": "n", "label": "Olay sayısı", "digits": 0},
            {"property": "duration_days", "label": "Süre", "unit": "gün", "digits": 1},
            {"property": "mag_max", "label": "En büyük", "digits": 1},
            {"property": "depth_km", "label": "Ortalama derinlik", "unit": "km", "digits": 1},
            {"property": "D_m2_s", "label": "Hidrolik yayınım", "unit": "m²/s", "digits": 3},
            {"property": "radius_km", "label": "Yayılım yarıçapı", "unit": "km", "digits": 1}]},
        "note": "Sürü sınıflaması basit bir kurala dayanır; D konum hatalarına duyarlıdır.",
    }, points_geojson(swl, ["t_start", "year", "n", "duration_days", "mag_max", "depth_km", "D_m2_s",
                            "log10_D", "radius_km"], note="Katalogdan otomatik sürü tespiti"))

    # ------------------------------------------------------------ 6. saha katmanları
    sites = load_sites()
    for spec, cats in site_layer_specs():
        s = sites[sites["category"].isin(cats)]
        write_layer(LAYERS, spec, points_geojson(s, list(sites.columns), note=spec["note"]))
    # veri bağlanınca dolacak katmanlar (dosya yokken haritada "veri yok" görünür)
    for lid, label, desc, prop, unit in (
            ("ext_isi_akisi", "Isı akısı", "IHFC Küresel Isı Akısı Veritabanı", "value", "mW/m²"),
            ("ext_radon", "Radon istasyonları", "Toprak gazı / su radonu (izinle alınacak)", "value", "Bq/m³"),
            ("ext_sicak_su", "Sıcak su kaynakları", "MTA envanteri (elle derlenecek)", "value", "°C"),
            ("ext_yeralti_suyu", "Yeraltı suyu kuyuları", "DSİ gözlem kuyuları (izinle)", "value", "m")):
        write_layer(LAYERS, {"id": lid, "label": label, "description": desc, "group": GROUP_FLUID, "order": 450,
                             "style": {"type": "circle", "radius": 5},
                             "color": {"property": prop, "type": "interpolate",
                                       "stops": [[0, "#3f7be0"], [1, "#ffe680"]]},
                             "legend": {"title": label, "unit": unit},
                             "note": "Veri henüz bağlanmadı; bkz. docs/akiskan_ve_dis_etkiler.md."}, None)

    md += ["## 4. İnsan kaynaklı sahalar", "",
           f"{len(sites)} saha (baraj, maden, jeotermal, petrol/gaz, gaz deposu, doğal gaz çıkışı) "
           "başlangıç/bitiş yılı ve durumuyla `external_data/insan_kaynakli_sahalar.csv` dosyasında. "
           "Konumların çoğu yaklaşık ve doğrulanmamış; özellik tablosuna yalnızca o yıl faal olan sahalar girer.", "",
           sites["category"].value_counts().rename_axis("tür").reset_index(name="sayı").to_markdown(index=False), "",
           "## 5. Veri kaynakları ve erişim", "", DATA_SOURCES.to_markdown(index=False), ""]

    # ------------------------------------------------------------ 7. fay elemanları
    el_path = ROOT / "data/processed/fault_elements.csv.gz"
    if el_path.exists():
        from deprem.external.elements import element_external
        el = pd.read_csv(el_path)
        ext = element_external(el, modern, modern["time"].max() + pd.Timedelta(seconds=1), mc,
                               sites=sites, phases=ph)
        from deprem.external.params import PARAMS
        from deprem.faults.fill import fill_missing, fill_summary
        el = el.drop(columns=[c for c in ext.columns if c in el.columns]).join(ext)
        # uzak tetikleme β'sı hücre düzeyinde hesaplandı; elemana bulunduğu hücreden aktarılır
        g = grid.set_index("cell_id")
        ecell = assign_cell(el.rename(columns={"lat": "latitude", "lon": "longitude"}), TURKEY, cells)
        ok = (ecell >= 0)
        beta = np.full(len(el), np.nan)
        b = g["dyn_beta"].where(g["dyn_n_total"] >= 20)
        beta[ok] = b.reindex(ecell[ok]).to_numpy()
        el["dyn_trigger_beta"] = beta
        # "50 km içinde baraj yok" anlamlı bir boşluktur, doldurulmaz
        skip = {"station_dist_km", "reservoir_years_since_fill"}
        cols = [n for n, *_ in PARAMS if n not in skip] + [
            "heat_flow_mw_m2", "vp_vs_ratio", "curie_depth_km", "crustal_thickness_km"]
        el = fill_missing(el, cols)
        el.to_csv(ROOT / "data/processed/fault_elements_external.csv.gz", index=False)
        summ = fill_summary(el, cols)
        summ.to_csv(tab / "external_element_fill.csv", index=False)
        md += ["## 6. Fay elemanlarına aktarım ve eksik değer doldurma", "",
               f"{len(el):,} fay elemanının (1 km × 1 km) her biri için gaz/sıvı ve dış etki parametreleri "
               "hesaplandı. Eksik değerler `deprem.faults.fill.fill_missing` ile önce aynı fay boyunca sağ/sol "
               "komşulardan, sonra komşu faylardan, sonra bölgesel ortancadan, hiç ölçüm yoksa literatür "
               "varsayılanından dolduruldu; her değerin kaynağı `<ad>__fill`, belirsizliği `<ad>__sd` "
               "sütunundadır. Doldurulan değerler eğitimde öncül olarak kullanılıp düzeltilecektir. "
               "Hiç ölçümü olmayan ve makul bir literatür değeri bulunmayan büyüklükler (ör. radon) uydurulmaz, "
               "'veri yok' kalır. Oranlar:", "", summ.to_markdown(index=False), ""]
    (ROOT / "reports" / "external_summary.md").write_text("\n".join(md) + "\n")
    log.info("bitti")


if __name__ == "__main__":
    main()

"""Fay fay biriken enerji raporu: reports/fault_energy.md

    python scripts/faultwalk.py partial_deficit partial_bval partial_chain partial_energy etas partial
    python scripts/faultwalk.py --energy-nulls 10
    python scripts/fault_energy_report.py

Yöntem: src/deprem/fault_energy.py. HİÇBİR ÇIKTI DEPREM TAHMİNİ DEĞİLDİR.
"""

from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from faultwalk_report import boot, holm, load, window  # noqa: E402

RUN = ROOT / "runs" / "faultwalk"
REP = ROOT / "reports" / "fault_energy.md"
TAB = ROOT / "reports" / "tables"
FIG = ROOT / "reports" / "figures"
BIG = 6.5

NEW = [("partial_deficit", "moment açığı oranı (küçük depremler dahil)", "deficit"),
       ("partial_bval", "yerel b-değeri farkı (son 5 yıl)", "bval"),
       ("partial_chain", "yüklü komşu fay zinciri uzunluğu", "chain"),
       ("partial_energy", "üçü birlikte", "deficit,bval,chain"),
       ("partial_deficit_km", "moment açığı oranı, 1 km parçalarla", "deficit_km"),
       ("partial_chain_km", "yüklü kesintisiz fay boyu, 1 km parçalarla", "chain_km"),
       ("partial_energy_km", "üçü birlikte, 1 km parçalarla", "deficit_km,bval,chain_km")]
KM_PAIRS = [("partial_deficit_km", "partial_deficit", "moment açığı"),
            ("partial_chain_km", "partial_chain", "yüklü zincir"),
            ("partial_energy_km", "partial_energy", "üçü birlikte")]
PERIODS = [("1900-2025", 1900, 2026), ("1900-1959", 1900, 1960), ("1960-1999", 1960, 2000),
           ("2000-2025", 2000, 2026)]


def f3(x):
    return f"{x:+.3f}"


def mw_of(m0):
    return (np.log10(np.maximum(m0, 1.0)) - 9.1) / 1.5


# ------------------------------------------------------------------ büyüklük testi
def _logit_fit(X, y, sd=1.0):
    from scipy.optimize import minimize
    Xb = np.column_stack([np.ones(len(y)), X])

    def nll(w):
        z = Xb @ w
        return (np.logaddexp(0, z) - y * z).sum() + 0.5 * ((w[1:] / sd) ** 2).sum()

    def grad(w):
        p = 1 / (1 + np.exp(-(Xb @ w)))
        g = Xb.T @ (p - y)
        g[1:] += w[1:] / sd ** 2
        return g
    return minimize(nll, np.zeros(Xb.shape[1]), jac=grad, method="L-BFGS-B").x


def magnitude_test(ev: pd.DataFrame, feats: dict, start_n: int = 60) -> pd.DataFrame:
    """Sıralı: her deprem için, yalnız önceki depremlerle eğitilmiş lojistik model
    "bu deprem M≥6,5 olacak mı" olasılığı verir. Puan: log olasılık; sabit oran (Gutenberg-Richter
    karşılığı) ile karşılaştırılır."""
    y = (ev["mw"].to_numpy() >= BIG).astype(float)
    out = {"i": ev["i"].to_numpy(), "year": ev["year"].to_numpy(), "mw": ev["mw"].to_numpy()}
    n = len(y)
    for name, X in [("sabit", None)] + list(feats.items()):
        sc = np.full(n, np.nan)
        for k in range(start_n, n):
            if X is None:
                p = (y[:k].sum() + 0.5) / (k + 1.0)
            else:
                w = _logit_fit(X[:k], y[:k])
                p = 1 / (1 + np.exp(-(w[0] + X[k] @ w[1:])))
            p = np.clip(p, 1e-4, 1 - 1e-4)
            sc[k] = np.log(p) if y[k] else np.log(1 - p)
        out[name] = sc
    return pd.DataFrame(out).dropna()


MARMARA = dict(lat=(39.8, 41.5), lon=(26.0, 31.6))


def _in_box(lat, lon, box=MARMARA):
    lat, lon = np.asarray(lat), np.asarray(lon)
    return (lat >= box["lat"][0]) & (lat <= box["lat"][1]) & (lon >= box["lon"][0]) & (lon <= box["lon"][1])


def marmara_section(W, E, rng) -> list:
    """Marmara: 1 km parçalarla Kuzey Anadolu Fayı kuzey kolu boyunca yük profili ve model sınaması."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from deprem.faultwalk_data import year_start_days
    L = ["## 6. Marmara bölgesi\n"]
    km, seg = E.km, E.seg
    el = km.el
    sys_el = seg["system"].to_numpy()[el["seg"]]
    m = _in_box(el["lat"], el["lon"]) & (sys_el == "Kuzey Anadolu Fayı")
    tau = E.tau[:-1]
    s_now = len(tau) - 1

    def step_at(date):
        return int(np.searchsorted(tau, year_start_days(int(date[:4])) + (pd.Timestamp(date) - pd.Timestamp(date[:4] + "-01-01")).days, side="right") - 1)

    snaps = [("2026 başı", s_now), ("16 Ağustos 1999 (İzmit'ten bir gün önce)", step_at("1999-08-16")),
             ("8 Ağustos 1912 (Mürefte-Şarköy'den bir gün önce)", step_at("1912-08-08"))]
    o = np.argsort(el["lon"].to_numpy()[m])
    lon = el["lon"].to_numpy()[m][o]
    lat = el["lat"].to_numpy()[m][o]
    fig, axs = plt.subplots(len(snaps), 1, figsize=(10, 7.5), sharex=True)
    cat = W.cat
    big = cat[(cat["mw"] >= 6.5) & _in_box(cat["latitude"], cat["longitude"])]
    from scipy.spatial import cKDTree
    from deprem.faultwalk_data import xy_km
    tx, ty = xy_km(el["lat"].to_numpy()[m], el["lon"].to_numpy()[m])
    bx, by = xy_km(big["latitude"].to_numpy(), big["longitude"].to_numpy())
    dd, _ = cKDTree(np.column_stack([tx, ty])).query(np.column_stack([bx, by]))
    big = big[dd <= 20]                                   # yalnız kuzey kol üzerindekiler
    for ax, (lab, s) in zip(axs, snaps):
        R = km.R[s][m][o]
        ax.scatter(lon, R, s=3, c=np.where(R >= 0.5, "#c0392b", "#7f8c8d"))
        ax.axhline(0.5, color="#555", lw=0.6, ls="--")
        ax.set_ylabel("doluluk R")
        ax.set_title(lab, fontsize=10, loc="left")
        t_s = tau[s]
        for _, r in big[big["t"] > t_s].head(3).iterrows():
            half = 10 ** (-2.44 + 0.59 * r["mw"]) / 2 / (111.195 * np.cos(np.radians(40.7)))
            ax.axvspan(r["longitude"] - half, r["longitude"] + half, color="#f39c12", alpha=0.15)
            ax.text(r["longitude"], ax.get_ylim()[1] * 0.9, f"{r['time'].year} M{r['mw']:.1f}", fontsize=7, ha="center")
    axs[-1].set_xlabel("boylam (batı → doğu: Saros, Ganos, Marmara Denizi, İzmit, Düzce)")
    fig.suptitle("Kuzey Anadolu Fayı, Marmara: 1 km parçalarda biriken yük (kırmızı: R ≥ 0,5). "
                 "Turuncu: sonraki büyük depremler (kaba kırık boyu)", fontsize=9)
    fig.tight_layout()
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "fault_energy_marmara.png", dpi=130)
    plt.close(fig)
    pd.DataFrame({"lat": lat.round(4), "lon": lon.round(4), "R_2026": km.R[s_now][m][o].round(3),
                  "yüklü_boy_km_2026": km.chain_km[s_now][m][o].round(0)}).to_csv(
        TAB / "fault_energy_marmara_1km.csv", index=False)
    Rn = km.R[s_now][m][o]
    L.append(f"Marmara kutusunda (enlem {MARMARA['lat'][0]}-{MARMARA['lat'][1]}, boylam {MARMARA['lon'][0]}-"
             f"{MARMARA['lon'][1]}) Kuzey Anadolu Fayı kuzey kolu {m.sum()} adet 1 km parçaya bölündü. 2026 başında "
             f"bunların %{100 * (Rn >= 0.5).mean():.0f}'i yüklü (R ≥ 0,5); en uzun kesintisiz yüklü parça "
             f"{km.chain_km[s_now][m].max():.0f} km. Tablo: `reports/tables/fault_energy_marmara_1km.csv`.\n")
    L.append("![Marmara yük profili](figures/fault_energy_marmara.png)\n")
    # segment tablosu
    fm = _in_box(seg["mid_lat"], seg["mid_lon"]) & (seg["slip_rate_mm_yr"] >= 2)
    el_seg = el["seg"].to_numpy()
    rows = []
    for k in np.nonzero(fm.to_numpy())[0]:
        e = el_seg == k
        rows.append({"segment": seg["fault_id"].iloc[k], "sistem": seg["system"].iloc[k],
                     "orta boylam": round(seg["mid_lon"].iloc[k], 2), "uzunluk_km": round(seg["length_km"].iloc[k]),
                     "kayma_hızı": seg["slip_rate_mm_yr"].iloc[k], "segment R": round(E.R[s_now, k], 2),
                     "1 km R ortalama": round(float(np.maximum(km.R[s_now, e], 0).mean()), 2),
                     "1 km yüklü pay %": round(100 * float((km.R[s_now, e] >= 0.5).mean())),
                     "b farkı z": round(E.bz[s_now, k], 2), "b olay sayısı": int(E.bn[s_now, k])})
    L.append("Marmara'daki segmentler (kayma hızı ≥ 2 mm/yıl), 2026 başı. Segment R tüm segmentin ortalaması; "
             "1 km sütunları, segmentin depremlerle kısmen boşalmış parçalarını ayrı sayar.\n")
    L.append(pd.DataFrame(rows).sort_values("orta boylam").to_markdown(index=False) + "\n")
    L.append("Dikkat: katalog 1700'de başlıyor. Marmara'da 1509 ve daha önceki büyük depremler hesaba girmiyor; "
             "1766 depremlerinin kırık uzunlukları tarihsel kaynaklarda tartışmalı, burada yalnız büyüklükten kaba "
             "tahmin edildi. Bu yüzden 1 km çözünürlük, yükün nerede biriktiğini değil, bu varsayımlarla nerede "
             "birikmiş görüneceğini gösterir.\n")
    # model sınaması Marmara alt kümesi
    part = load("partial")
    base = load("etas")
    if part is not None and base is not None:
        rows, ps = [], []
        for n, lab, _ in NEW:
            d = load(n)
            if d is None:
                continue
            for ref_name, ref in (("ETAS", base), ("ana model", part)):
                mm = d.merge(ref[["i", "score"]], on="i", suffixes=("", "_r"))
                mm = window(mm[_in_box(mm["lat"], mm["lon"])], 1900, 2026)
                if len(mm) < 5:
                    continue
                mean, lo, hi, p = boot((mm["score"] - mm["score_r"]).to_numpy(), mm["year"].to_numpy(), rng)
                rows.append({"model": "Ana model + " + lab, "karşı": ref_name, "Marmara 1900-2025": f3(mean),
                             "%95 GA": f"[{lo:+.3f}, {hi:+.3f}]", "p": p})
                ps.append(p)
        if rows:
            T = pd.DataFrame(rows)
            T["Holm p"] = [f"{x:.3f}" for x in holm(ps)]
            T["p"] = T["p"].map(lambda x: f"{x:.3f}")
            nM = int(_in_box(window(part, 1900, 2026)["lat"], window(part, 1900, 2026)["lon"]).sum())
            L.append(f"Model sınaması yalnız Marmara depremleriyle (1900-2025, {nM} deprem M≥5,5). Model tüm Türkiye "
                     "verisiyle eğitildi, puan yalnız Marmara'daki depremlerden. Bu kadar az depremle küçük farklar "
                     "şansla ayırt edilemez.\n")
            L.append(T.to_markdown(index=False) + "\n")
    return L


def main():
    from deprem.fault_energy import CHAIN_THR, energy_joule, length_to_mw
    from deprem.faultwalk import VARIANTS
    rng = np.random.default_rng(0)
    W = pickle.load(open(RUN / "walkdata.pkl", "rb"))
    E = pickle.load(open(RUN / "energy.pkl", "rb"))
    seg = E.seg
    S = E.tau.size - 1
    last = S - 1
    useg = W.units["seg"].to_numpy()
    fay_u = W.units["kind"].to_numpy() == "fay"
    L = []
    L.append("# Fay fay biriken enerji, küçük depremler ve komşu fay zinciri\n")
    L.append("> Deneysel araştırma çıktısıdır, deprem tahmini değildir. Hiçbir gösterge kanıtsız doğrulanmış "
             "sayılmaz.\n")

    # ---------------------------------------------------------- 1. bugünkü durum
    sr = seg["slip_rate_mm_yr"].fillna(0.2)
    T = pd.DataFrame({"segment": seg["fault_id"], "ad": seg.get("name"), "sistem": seg["system"],
                      "uzunluk_km": seg["length_km"].round(0), "kayma_hızı_mm_yıl": sr.round(1),
                      "açık_Mw_eşdeğeri": mw_of(E.deficit[last]).round(2),
                      "enerji_PJ": (energy_joule(np.maximum(E.deficit[last], 0)) / 1e15).round(1),
                      "doluluk_R": E.R[last].round(2),
                      "yüklü_zincir_km": E.chain_km[last].round(0),
                      "zincirin_olası_Mw": length_to_mw(E.chain_km[last]).round(1),
                      "b_farkı_z": E.bz[last].round(2), "b_olay_sayısı": E.bn[last],
                      "küçük_depremlerin_boşalttığı_pay_%": np.where(E.released_all > 0, 100 * E.released_small /
                                                                     np.maximum(E.released_all, 1), np.nan).round(2),
                      "lat": seg["mid_lat"].round(3), "lon": seg["mid_lon"].round(3)})
    TAB.mkdir(parents=True, exist_ok=True)
    T.to_csv(TAB / "fault_energy_2026.csv", index=False)
    from deprem.fault_energy import load_all_events
    ev_all = load_all_events(ROOT)
    tot_small = E.released_small.sum()
    tot_all = E.released_all.sum()
    load_rate = (3e10 * seg["area_km2"] * 1e6 * sr * 1e-3).to_numpy()
    yrs_small = tot_small / load_rate.sum()
    L.append("## 1. Küçük depremler fayı boşaltıyor mu?\n")
    L.append(f"1700'den 2025 sonuna kadar fay segmentlerine düşen toplam deprem momentinin yalnız "
             f"**%{100 * tot_small / tot_all:.2f}**'i M<5,5 depremlerden geliyor. Kayıtlı bütün küçük depremler "
             f"({int((ev_all['mw'] < 5.5).sum()):,} deprem) toplamda, Türkiye faylarının yalnız ≈{yrs_small:.1f} "
             "yılda biriktirdiği yükü boşaltmış. Yani küçük depremler fayı boşaltmıyor, yük büyük depremlere "
             "kalıyor. Küçük depremler enerjiyi değil, gerilmenin nasıl değiştiğine dair ipucunu (b-değeri) "
             "taşıyabilir. (Küçük depremlerin kaydı 1998'den önce eksik; eksik olanlar eklense de bu oran "
             "birkaç kattan fazla değişmez, çünkü her büyüklük basamağı bir öncekinin ~30 katı enerji taşır.)\n")

    # en dolu segmentler
    sel = T[(T["kayma_hızı_mm_yıl"] >= 2) & (T["uzunluk_km"] >= 10)].sort_values("doluluk_R", ascending=False)
    L.append("## 2. Fay fay biriken yük (2026 başı)\n")
    L.append("Doluluk R = biriken moment açığı / fayın kaldırabileceği (10 MPa gerilme düşümüyle). Saat 1700'de "
             "sıfırdan başlıyor; 1700'den beri kırılmamış faylarda gerçek değer daha yüksek olabilir. R kalibre "
             "edilmemiş, göreli bir ölçüdür: R = 1, fayın taşıyabileceği en büyük gerilme düşümüne denk yük demektir, "
             "fayın o an kırılacağı anlamına gelmez. Kayma hızı ≥ 2 mm/yıl, uzunluk ≥ 10 km segmentler, en dolu 20. "
             "Tamamı: `reports/tables/fault_energy_2026.csv`.\n")
    cols = ["segment", "sistem", "uzunluk_km", "kayma_hızı_mm_yıl", "doluluk_R", "açık_Mw_eşdeğeri",
            "enerji_PJ", "yüklü_zincir_km", "zincirin_olası_Mw", "b_farkı_z"]
    L.append(sel[cols].head(20).to_markdown(index=False) + "\n")
    L.append("enerji_PJ: biriken elastik enerji (petajoule; 1 PJ ≈ 240 kiloton TNT), E ≈ 3 MPa/(2μ) × moment "
             "açığı. açık_Mw_eşdeğeri: biriken açığın tamamı tek depremde boşalsa büyüklüğü.\n")

    # ---------------------------------------------------------- 3. büyük depremler
    L.append("## 3. Büyük depremler (M≥6,5, 1900 sonrası): kırılmadan hemen önce fay ne durumdaydı?\n")
    L.append(f"Yüklü = R ≥ {CHAIN_THR}. Zincir: depremin düştüğü segmente uç uca bağlı (≤ 5 km) yüklü segmentlerin "
             "toplam boyu. Gerçek kırık boyu büyüklükten (Wells & Coppersmith 1994) kaba tahmin.\n")
    from scipy.spatial import cKDTree
    from deprem.faultwalk_data import xy_km
    km = E.km
    exk, eyk = xy_km(km.el["lat"].to_numpy(), km.el["lon"].to_numpy())
    ktree = cKDTree(np.column_stack([exk, eyk]))
    rows = []
    for i in np.nonzero((W.m_ev >= BIG) & (W.cat["time"].dt.year.to_numpy() >= 1900))[0]:
        u = W.ev_unit[i]
        if not fay_u[u]:
            rows.append({"tarih": W.cat["time"].iloc[i].date(), "Mw": W.m_ev[i], "segment": "fay dışı"})
            continue
        k = useg[u]
        s = max(W.ev_step[i], 0)
        Lk = seg["length_km"].iloc[k]
        rows.append({"tarih": W.cat["time"].iloc[i].date(), "Mw": W.m_ev[i], "segment": seg["fault_id"].iloc[k],
                     "sistem": seg["system"].iloc[k], "doluluk_R": round(E.R[s, k], 2),
                     "segment_boyu_km": round(Lk), "segmentin_tek_başına_Mw": round(float(length_to_mw(Lk)), 1),
                     "yüklü_zincir_km": round(E.chain_km[s, k]),
                     "zincirin_olası_Mw": round(float(length_to_mw(E.chain_km[s, k])), 1) if E.chain_km[s, k] > 0 else "",
                     "gerçek_kırık_km (tahmini)": round(10 ** (-2.44 + 0.59 * W.m_ev[i])),
                     "b_farkı_z": round(E.bz[s, k], 2)})
        qx, qy = xy_km(W.cat["latitude"].iloc[i], W.cat["longitude"].iloc[i])
        dk, e = ktree.query([float(qx), float(qy)])
        if dk <= 25:
            rows[-1]["1km_doluluk"] = round(float(km.R[s, e]), 2)
            rows[-1]["1km_yüklü_boy_km"] = round(float(km.chain_km[s, e]))
    BT = pd.DataFrame(rows)
    L.append(BT.to_markdown(index=False) + "\n")
    bf = BT[BT["segment"] != "fay dışı"]
    n_loaded = int((bf["doluluk_R"] >= CHAIN_THR).sum())
    n_km = int((BT.get("1km_doluluk", pd.Series(dtype=float)) >= CHAIN_THR).sum())
    n_kmall = int(BT.get("1km_doluluk", pd.Series(dtype=float)).notna().sum())
    L.append(f"{len(bf)} fay üstü büyük depremin {n_loaded}'inde depremin başladığı segment yüklüydü. 1 km "
             f"parçalarla bakıldığında (merkez üssüne en yakın parça, 25 km içinde) {n_kmall} depremin {n_km}'inde "
             "parça yüklüydü. Büyük depremlerin çoğu, bu hesapla 'boş' görünen yavaş faylarda oldu: bu faylarda "
             "deprem tekrarı binlerce yıl, saatimiz ise 1700'de başlıyor; fayın 1700'deki yükünü bilmiyoruz.\n")

    # ---------------------------------------------------------- 4. büyüklük testi
    L.append("## 4. \"Biri sınırdaysa yeter mi, yoksa boydan boya yüklü olmalı mı?\" sınaması\n")
    km = E.km
    el_seg = km.el["seg"].to_numpy()
    seg_el = {k: np.nonzero(el_seg == k)[0] for k in np.unique(el_seg)}
    ev = []
    for i in range(len(W.m_ev)):  # noqa: B007
        y = W.cat["time"].iloc[i].year
        u = W.ev_unit[i]
        if y < 1900 or not fay_u[u]:
            continue
        k = useg[u]
        s = max(W.ev_step[i], 0)
        Lk = seg["length_km"].iloc[k]
        loaded = E.R[s, k] >= CHAIN_THR
        ev.append({"i": i, "year": y, "mw": W.m_ev[i],
                   "seed": np.log10(max(Lk if loaded else 5.0, 5.0) / 100),
                   "chain": np.log10(max(E.chain_km[s, k], 5.0) / 100),
                   "R": np.clip(np.log(max(E.R[s, k], 1e-3)), -3, 1.5),
                   "seglen": np.log10(max(Lk, 5.0) / 100),
                   "chain_km": np.log10(max(float(km.chain_km[s, seg_el[k]].max()), 5.0) / 100)})
    ev = pd.DataFrame(ev)
    feats = {"segment boyu (yükten bağımsız)": ev[["seglen"]].to_numpy(),
             "yalnız kendi segmenti (yüklüyse boyu)": ev[["seed"]].to_numpy(),
             "yüklü komşu zinciri": ev[["chain"]].to_numpy(),
             "zincir + kendi doluluğu": ev[["chain", "R"]].to_numpy(),
             "1 km parçalarla yüklü kesintisiz boy": ev[["chain_km"]].to_numpy()}
    mt = magnitude_test(ev, feats)
    mrows, mps = [], []
    for name in feats:
        d = (mt[name] - mt["sabit"]).to_numpy()
        mean, lo, hi, p = boot(d, mt["year"].to_numpy(), rng)
        big = mt[mt["mw"] >= BIG]
        mrows.append({"bilgi": name, "olay başına kazanç": f3(mean), "%95 GA": f"[{lo:+.3f}, {hi:+.3f}]", "p": p,
                      "M≥6,5 olaylarda": f3((big[name] - big["sabit"]).mean())})
        mps.append(p)
    MT = pd.DataFrame(mrows)
    MT["Holm p"] = [f"{x:.3f}" for x in holm(mps)]
    MT["p"] = MT["p"].map(lambda x: f"{x:.3f}")
    L.append(f"Fay üstündeki her M≥5,5 deprem için (1900-2025, ilk 60'ı yalnız eğitim, {len(mt)} deprem puanlandı, "
             f"{int((mt['mw'] >= BIG).sum())}'i M≥6,5): yalnız önceki depremlerle eğitilmiş bir model \"bu deprem "
             "M≥6,5 olacak mı\" olasılığı veriyor. Karşılaştırma: büyüklüğü hiç bilmeyen sabit oran "
             "(Gutenberg-Richter). Pozitif = o bilgi büyük depremleri ayırt etmeye yardım ediyor.\n")
    L.append(MT.to_markdown(index=False) + "\n")

    # ---------------------------------------------------------- 5. sıralı model
    base = load("etas")
    part = load("partial")
    L.append("## 5. Göstergeler deprem modeline eklenince (ne zaman ve nerede)\n")
    if base is not None and part is not None:
        rows, ps = [], []
        for n, lab, _ in NEW:
            d = load(n)
            if d is None:
                continue
            for ref_name, ref in (("ETAS", base), ("ana model", part)):
                m = d.merge(ref[["i", "score"]], on="i", suffixes=("", "_r"))
                r = {"model": "Ana model + " + lab, "karşı": ref_name}
                for plab, a, b in PERIODS:
                    w = window(m, a, b)
                    diff = (w["score"] - w["score_r"]).to_numpy()
                    if plab == "1900-2025":
                        mean, lo, hi, p = boot(diff, w["year"].to_numpy(), rng)
                        r["1900-2025"] = f3(mean)
                        r["%95 GA"] = f"[{lo:+.3f}, {hi:+.3f}]"
                        r["p"] = p
                        ps.append(p)
                    else:
                        r[plab] = f3(diff.mean()) if len(diff) else ""
                w = window(m, 1900, 2026)
                bg = w[w["mw"] >= BIG]
                r["M≥6,5"] = f3((bg["score"] - bg["score_r"]).mean())
                rows.append(r)
        RT = pd.DataFrame(rows)
        RT["Holm p"] = [f"{x:.3f}" for x in holm(ps)]
        RT["p"] = RT["p"].map(lambda x: f"{x:.3f}")
        nm = len(window(base, 1900, 2026))
        n2k = len(window(base, 2000, 2026))
        L.append(f"Her deprem, yalnız önceki verilerle yeniden eğitilmiş modelle puanlandı ({nm} deprem 1900-2025; "
                 f"{n2k}'i 2000-2025). Olay başına log olasılık farkı; GA yıl blokları bootstrap'ı; Holm düzeltmesi "
                 "tablodaki tüm karşılaştırmalar için. b-değeri girdisi küçük deprem kaydı yeterli olan yıllarda "
                 "(pratikte 1998 sonrası) doludur, öncesinde etkisizdir; bu yüzden 2000-2025 sütununa ayrıca bakın.\n")
        L.append(RT.to_markdown(index=False) + "\n")
        # 2000+ b-değeri için ayrı bootstrap
        d = load("partial_bval")
        if d is not None:
            m = window(d.merge(part[["i", "score"]], on="i", suffixes=("", "_r")), 2000, 2026)
            mean, lo, hi, p = boot((m["score"] - m["score_r"]).to_numpy(), m["year"].to_numpy(), rng)
            L.append(f"b-değeri, yalnız 2000-2025 ({len(m)} deprem), ana modele göre: {f3(mean)} "
                     f"[{lo:+.3f}, {hi:+.3f}], p = {p:.3f} (düzeltmesiz).\n")
        # öğrenilen katsayılar
        krows = []
        for n, lab, cv in NEW:
            pj = RUN / f"{n}.params.json"
            if not pj.exists():
                continue
            snap = json.load(open(pj))[-1]
            for k, c in enumerate(cv.split(",")):
                key = ("k_cov", "k_cov1", "k_cov2")[k]
                krows.append({"model": lab, "girdi": c, "katsayı (2025 sonu)": round(snap.get(key, np.nan), 3)})
        if krows:
            L.append("Öğrenilen katsayılar (+ = gösterge yüksekken o fayda deprem daha olası). Girdiler: moment açığı "
                     "log R; b farkı z; zincir log10(km/100).\n")
            L.append(pd.DataFrame(krows).to_markdown(index=False) + "\n")

    # 1 km mi segment mi?
    if part is not None:
        krows = []
        for a, b, lab in KM_PAIRS:
            da, db = load(a), load(b)
            if da is None or db is None:
                continue
            m = window(da.merge(db[["i", "score"]], on="i", suffixes=("", "_r")), 1900, 2026)
            mean, lo, hi, p = boot((m["score"] - m["score_r"]).to_numpy(), m["year"].to_numpy(), rng)
            bg = m[m["mw"] >= BIG]
            krows.append({"gösterge": lab, "1 km − segment (1900-2025)": f3(mean), "%95 GA": f"[{lo:+.3f}, {hi:+.3f}]",
                          "p": f"{p:.3f}", "M≥6,5": f3((bg["score"] - bg["score_r"]).mean())})
        if krows:
            L.append("### Fayı 1 km'lik parçalara bölmek, segment bütününe göre işe yarıyor mu?\n")
            L.append(f"Her fay {len(km.el):,} parçaya bölündü ({km.el['len_km'].iloc[0]:.0f} km). Her parçada yük, "
                     "kayma hızıyla her yıl artıyor; bir deprem o parçaya ulaştığında (merkezden kırık yarı boyu + 5 km) "
                     "momenti kadar yük düşüyor, yük sıfırın altına inmiyor (\"son kırılmadan beri biriken\"). "
                     "Pozitif = 1 km sürümü daha iyi.\n")
            L.append(pd.DataFrame(krows).to_markdown(index=False) + "\n")

    # şans kontrolü
    nulls = sorted(RUN.glob("null_energy_*.csv"))
    pe, pp = load("partial_energy_e20"), load("partial_e20")
    if nulls and pe is not None and pp is not None:
        def gain(d):
            m = window(d.merge(pp[["i", "score"]], on="i", suffixes=("", "_r")), 1900, 2026)
            return (m["score"] - m["score_r"]).mean()
        real = gain(pe)
        nv = np.array([gain(load(p.stem)) for p in nulls])
        L.append("### Şans kontrolü\n")
        L.append(f"Üç gösterge fay birimleri arasında rastgele karıştırıldı (her fay başka bir fayın göstergesini "
                 f"aldı; model her 20 depremde bir eğitildi). Gerçek katkı {f3(real)}; {len(nv)} karıştırmada "
                 f"ortalama {f3(nv.mean())}, en yüksek {f3(nv.max())}. Şans p = "
                 f"{(1 + (nv >= real).sum()) / (1 + len(nv)):.2f}.\n")

    L.extend(marmara_section(W, E, rng))

    L.append("## Yöntem ve sınırlar\n")
    L.append("- Yöntem: `src/deprem/fault_energy.py`; modele ekleme: `src/deprem/faultwalk.py` (`partial_*` sürümleri).\n"
             "- Kayma hızları GEM aktif fay veritabanından (jeolojik ve jeodezik derlemesi); GPS kilitlenme oranı bu "
             "sürümde kullanılmadı, krip yapan faylarda yük olduğundan fazla görünür.\n"
             "- Deprem momenti, merkez üssünden kırık yarı boyu + 5 km içindeki fay izlerine dağıtıldı; tek yönlü "
             "kırılan büyük depremlerde (ör. 1939, 2023) boşalımın bir kısmı yanlış segmente düşer.\n"
             "- 1700-1899 kataloğu eksik; o dönemde boşalım az sayılır, doluluk olduğundan yüksek görünür.\n"
             "- Kapasite tek bir gerilme düşümüyle (10 MPa) hesaplandı; gerçek değer fay fay 1-30 MPa arasında "
             "değişebilir. Gerilme düşümünü fay fay ölçmek için dalga biçimi spektrumları gerekir; bu sürümde yok.\n"
             "- Yüklülük eşiği (R ≥ 0,5) ve b penceresi (5 yıl, 30 km, ≥ 50 olay) sonuçlara bakılmadan önceden "
             "seçildi.\n")
    REP.write_text("\n".join(L), encoding="utf-8")
    print(REP)


if __name__ == "__main__":
    main()

"""Fay fay, depremden depreme öğrenen modelin raporu: reports/faultwalk.md

    python scripts/faultwalk_report.py
"""

from __future__ import annotations

import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "runs" / "faultwalk"
REP = ROOT / "reports" / "faultwalk.md"
FIG = ROOT / "reports" / "figures"

PERIODS = [("1700-1899 (eksik katalog)", 1700, 1900), ("1900-1959", 1900, 1960), ("1960-1999", 1960, 2000),
           ("2000-2025", 2000, 2026)]
MAIN = ("1900-2025", 1900, 2026)
COMPARE = ["poisson", "pooled", "partial", "separate", "fixed", "no_interact", "no_cross", "no_cfs", "no_slow",
           "no_far", "no_load", "partial_dvv"]
B = 4000


ABLATE = [("no_far", "100 km ötesindeki fayların etkisi"),
          ("no_interact", "faylar arası etkinin tümü (sistem, Coulomb, yavaş)"),
          ("no_slow", "yavaş (yıllar süren) etki"),
          ("no_cross", "farklı büyük sistemler arası ayrı çarpan"),
          ("no_cfs", "Coulomb gerilme yönü"),
          ("no_load", "fay doluluğu (son büyük depremden beri biriken)"),
          ("pooled", "fay fay ayar (ortak ayara karşı)"),
          ("fixed", "sürekli güncelleme (1900'de dondurmaya karşı)"),]


def load(name):
    p = RUN / f"{name}.csv"
    if not p.exists():
        return None
    d = pd.read_csv(p)
    d["time"] = pd.to_datetime(d["time"], utc=True, format="ISO8601")
    d["year"] = d["time"].dt.year
    return d


def boot(diff: np.ndarray, years: np.ndarray, rng, n=B):
    """Yıl blokları bootstrap: ortalama, %95 GA, tek yönlü p (fark ≤ 0)."""
    uy = np.unique(years)
    by = [diff[years == y] for y in uy]
    sums = np.array([b.sum() for b in by])
    cnts = np.array([b.size for b in by])
    idx = rng.integers(0, uy.size, (n, uy.size))
    means = sums[idx].sum(1) / cnts[idx].sum(1)
    m = diff.mean()
    # p: merkezlenmiş bootstrap (sıfır hipotezi: ortalama fark 0)
    p = float((means - m >= m).mean()) if m > 0 else 1.0
    return m, np.percentile(means, 2.5), np.percentile(means, 97.5), max(p, 1.0 / n)


def holm(ps):
    ps = np.asarray(ps, float)
    o = np.argsort(ps)
    adj = np.empty_like(ps)
    run = 0.0
    for k, i in enumerate(o):
        run = max(run, min(1.0, (ps.size - k) * ps[i]))
        adj[i] = run
    return adj


def window(d, a, b):
    return d[(d["year"] >= a) & (d["year"] < b)]


def f3(x):
    return f"{x:+.3f}"


def main():
    rng = np.random.default_rng(0)
    base = load("etas")
    labels = {}
    from deprem.faultwalk import VARIANTS
    for k, v in VARIANTS.items():
        labels[k] = v.label
    labels["fixed"] = "Ana model, ayar 1900'de donduruldu (sabit ayarlı)"
    labels["partial_dvv"] = "Ana model + Doğu Marmara dv/v"
    res = {n: load(n) for n in COMPARE}
    res = {n: d for n, d in res.items() if d is not None}
    W = pickle.load(open(RUN / "walkdata.pkl", "rb"))
    L = []
    L.append("# Depremden depreme, fay fay öğrenen model\n")
    L.append("> Deneysel araştırma çıktısıdır, deprem tahmini değildir.\n")
    n_all = len(base)
    L.append(f"Katalog: {n_all} deprem (M≥5,5, Türkiye kutusu, {base['time'].min().year}-{base['time'].max().year}; "
             "1700-1899 tarihsel katalog, 1900 sonrası USGS/EMSC/AFAD birleşik).\n")
    L.append("Her deprem için model **yalnız o depremden önceki** depremlerle yeniden eğitildi, sonra o depremin "
             "ne zaman ve hangi 0,5° hücrede olacağına verdiği olasılıkla puanlandı. Puan farkı olay başına doğal "
             "log birimindedir: +0,10 ≈ gerçekleşen depremlere ETAS'tan %10 daha yüksek olasılık vermek.\n")

    # ana tablo
    L.append("## Sonuçlar (ETAS'a göre olay başına kazanç, aynı sıralı düzende)\n")
    rows, ps = [], []
    for n, d in res.items():
        m = d.merge(base[["i", "score", "time_score", "loc_score"]], on="i", suffixes=("", "_e"))
        r = {"model": labels.get(n, n)}
        w = window(m, MAIN[1], MAIN[2])
        diff = (w["score"] - w["score_e"]).to_numpy()
        mean, lo, hi, p = boot(diff, w["year"].to_numpy(), rng)
        r["1900-2025 kazanç"] = f3(mean)
        r["%95 GA"] = f"[{lo:+.3f}, {hi:+.3f}]"
        r["p (tek yönlü)"] = p
        ps.append(p)
        if "time_score" in w and w["time_score"].notna().all() and "time_score_e" in w:
            r["zaman"] = f3((w["time_score"] - w["time_score_e"]).mean())
            r["konum"] = f3((w["loc_score"] - w["loc_score_e"]).mean())
        for lab, a, b in PERIODS:
            ww = window(m, a, b)
            r[lab] = f3((ww["score"] - ww["score_e"]).mean()) if len(ww) else ""
        big = w[w["mw"] >= 6.5]
        r["M≥6,5 olaylar"] = f3((big["score"] - big["score_e"]).mean()) if len(big) else ""
        rows.append(r)
    tab = pd.DataFrame(rows)
    tab["Holm p"] = holm(ps)
    tab["p (tek yönlü)"] = tab["p (tek yönlü)"].map(lambda x: f"{x:.3f}")
    tab["Holm p"] = tab["Holm p"].map(lambda x: f"{x:.3f}")
    nmain = len(window(base, MAIN[1], MAIN[2]))
    L.append(f"1900-2025: {nmain} deprem. Güven aralığı yıl blokları bootstrap'ı, p değerleri Holm ile çoklu teste göre düzeltildi.\n")
    L.append(tab.to_markdown(index=False) + "\n")

    # bileşenlerin katkısı: ana model − (o bileşen çıkarılmış sürüm), doğrudan karşılaştırma
    if "partial" in res:
        L.append("## Her parçanın katkısı (ana model − o parça çıkarılmış sürüm, 1900-2025)\n")
        L.append("Pozitif = o parça işe yarıyor. Aynı depremler üzerinde eşli fark; GA yıl blokları bootstrap'ı.\n")
        pm = window(res["partial"], MAIN[1], MAIN[2])
        crows, cps = [], []
        for n, lab in ABLATE:
            if n not in res:
                continue
            m = pm.merge(res[n][["i", "score"]], on="i", suffixes=("", "_o"))
            diff = (m["score"] - m["score_o"]).to_numpy()
            mean, lo, hi, p = boot(diff, m["year"].to_numpy(), rng)
            crows.append({"parça": lab, "katkı": f3(mean), "%95 GA": f"[{lo:+.3f}, {hi:+.3f}]", "p": p})
            cps.append(p)
        if crows:
            ct = pd.DataFrame(crows)
            ct["Holm p"] = [f"{x:.3f}" for x in holm(cps)]
            ct["p"] = ct["p"].map(lambda x: f"{x:.3f}")
            L.append(ct.to_markdown(index=False) + "\n")

    # dv/v katkısı
    if "partial_dvv" in res and "partial" in res:
        L.append("## Doğu Marmara dv/v (kaya hız değişimi) eklenince\n")
        m = res["partial_dvv"].merge(res["partial"][["i", "score"]], on="i", suffixes=("", "_o"))
        w = m[m["year"] >= 2009]
        diff = (w["score"] - w["score_o"]).to_numpy()
        mean, lo, hi, p = boot(diff, w["year"].to_numpy(), rng)
        em = w[w["lat"].between(40.0, 41.2) & w["lon"].between(29.0, 31.6)]
        L.append(f"2009-2025, tüm Türkiye ({len(w)} deprem): olay başına {f3(mean)} [{lo:+.3f}, {hi:+.3f}], p = {p:.3f}. "
                 f"Doğu Marmara kutusunda bu dönemde M≥5,5 deprem sayısı: {len(em)}"
                 + (" (" + ", ".join(f"{r.time:%Y-%m-%d} M{r.mw:.1f}: {r.score - r.score_o:+.3f}" for r in em.itertuples()) + ")"
                    if len(em) else "") + ".\n")
        pj = RUN / "partial_dvv.params.json"
        if pj.exists():
            k = json.load(open(pj))[-1].get("k_cov")
            if k is not None:
                L.append(f"Öğrenilen dv/v katsayısı (2025 sonu): {k:+.3f} (dv/v 1 standart sapma yükselince bölgedeki fayların "
                         f"deprem hızı ×{np.exp(k):.2f}).\n")
        L.append("Bu büyüklük eşiğinde (M≥5,5) dv/v'nin kapsadığı yer ve zamanda çok az deprem olduğu için bu test "
                 "neredeyse hiç güç taşımaz; sonuç ne çıkarsa çıksın dv/v için kanıt sayılmamalı.\n")

    # şans kontrolleri
    L.append("## Şans kontrolü\n")
    e = None
    for f in sorted(RUN.glob("partial_e*.csv")):
        e = f.stem.split("_e")[-1]
    if e is not None:
        pe, po, nc = load(f"partial_e{e}"), load(f"pooled_e{e}"), load(f"no_cfs_e{e}")

        def g(a, b_):
            m = a.merge(b_[["i", "score"]], on="i", suffixes=("", "_b"))
            w = window(m, MAIN[1], MAIN[2])
            return (w["score"] - w["score_b"]).mean()
        nf = [load(f.stem) for f in sorted(RUN.glob("null_fault_*.csv"))]
        ncf = [load(f.stem) for f in sorted(RUN.glob("null_cfs_*.csv"))]
        real_f = g(pe, po)
        null_f = np.array([g(x, po) for x in nf])
        real_c = g(pe, nc)
        null_c = np.array([g(x, nc) for x in ncf])
        L.append(f"Bu testlerde model her {e} depremde bir yeniden eğitildi (maliyet için); gerçek ve karıştırılmış "
                 "sürümler aynı düzende.\n")
        if len(null_f):
            pf = (1 + (null_f >= real_f).sum()) / (1 + len(null_f))
            L.append(f"- **Fay fay ayarın katkısı** (fay fay − ortak, 1900-2025): gerçek {real_f:+.3f}; fay birimleri "
                     f"hücreler arasında rastgele karıştırıldığında {len(null_f)} denemede ortalama {null_f.mean():+.3f} "
                     f"(en yüksek {null_f.max():+.3f}). Şans p = {pf:.2f}.")
        if len(null_c):
            pc = (1 + (null_c >= real_c).sum()) / (1 + len(null_c))
            L.append(f"- **Coulomb yönünün katkısı** (ana − Coulomb'suz): gerçek {real_c:+.3f}; Coulomb işaretleri "
                     f"depremler arasında karıştırıldığında {len(null_c)} denemede ortalama {null_c.mean():+.3f} "
                     f"(en yüksek {null_c.max():+.3f}). Şans p = {pc:.2f}.")
        L.append("")

    # öğrenilen etkileşim
    L.append("## Faylar birbirini nasıl etkiliyor? (ana modelin 2025 sonu ayarı)\n")
    snaps = json.load(open(RUN / "partial.params.json")) if (RUN / "partial.params.json").exists() else []
    if snaps:
        s = snaps[-1]
        q = 1 + np.exp(s["log_qm1"])
        q2 = 1 + np.exp(s["log_q2m1"])
        rows = []
        for mw in (6.0, 7.0):
            D2 = np.exp(2 * s["log_D0"] + s["gam_m"] * (mw - 5.5))
            Ds = np.exp(2 * s["log_D2"] + s["gam_m"] * (mw - 5.5))
            K = np.exp(s["logK"] + min(max(s["alpha"], 0), 3.5) * (mw - 5.5))
            K2 = np.exp(s["logK2"] + min(max(s["alpha2"], 0), 3.5) * (mw - 5.5))
            for dkm in (10, 50, 100, 200, 400, 800):
                f1 = K * (q - 1) / (np.pi * D2) * (1 + dkm ** 2 / D2) ** (-q) * 2500
                f2 = K2 * (q2 - 1) / (np.pi * Ds) * (1 + dkm ** 2 / Ds) ** (-q2) * 2500
                rows.append({"büyüklük": mw, "uzaklık km": dkm, "hızlı etki (2500 km² alanda beklenen ek deprem)": f"{f1:.2e}",
                             "yavaş etki": f"{f2:.2e}"})
        L.append(f"- Hızlı etki mesafeyle (1 + d²/D²)^−q ile azalıyor, öğrenilen q = {q:.2f}; yavaş etki q = {q2:.2f}, "
                 f"süresi τ ≈ {np.exp(s['log_tau2']) / 365.25:.1f} yıl. Uzaklık sınırı konmadı.")
        L.append(f"- Aynı fay sistemi (veya aynı bölge) içinde etki çarpanı ×{np.exp(s['g_same']):.2f}; farklı büyük sistemler arasında "
                 f"×{np.exp(s['g_cross']):.2f} (hızlı), ×{np.exp(s['g_cross2']):.2f} (yavaş).")
        L.append(f"- Coulomb gerilmesi artan yerde etki çarpanı ×{np.exp(s['g_cfs']):.2f} (hızlı), ×{np.exp(s['g_cfs2']):.2f} "
                 "(yavaş); azalan yerde bunun tersi.")
        L.append(f"- Fay doluluğu katsayısı (ortak kısım) b0 = {s['b0']:+.3f} (+ ise dolu faylarda deprem daha olası).")
        L.append(f"- Öğrenilen katalog tespit oranları: 1700-1799 %{100 / (1 + np.exp(-s['era0'])):.0f}, 1800-1849 "
                 f"%{100 / (1 + np.exp(-s['era1'])):.0f}, 1850-1899 %{100 / (1 + np.exp(-s['era2'])):.0f} (1900 sonrası %100 varsayıldı).\n")
        L.append(pd.DataFrame(rows).to_markdown(index=False) + "\n")
        # fay başına ayar
        units = W.units.copy()
        units["a (kısmi)"] = np.array(s["A"])[units["sys"]] + np.array(s["D"])
        units["b (kısmi)"] = s["b0"] + np.array(s["B"])[units["sys"]] + np.array(s["E"])
        if (RUN / "separate.params.json").exists():
            s2 = json.load(open(RUN / "separate.params.json"))[-1]
            units["a (tamamen ayrı)"] = np.array(s2["D"])
        cnt = pd.Series(W.ev_unit).value_counts()
        units["deprem sayısı"] = units.index.map(cnt).fillna(0).astype(int)
        L.append("### Fay başına öğrenilen ayar (en çok depremli 15 birim)\n")
        L.append("a: fayın uzun dönem haritaya göre deprem üretme çarpanının logu (+ = haritanın beklediğinden çok). "
                 "Kısmi havuzlamada az depremli faylar ortak değere çekilir.\n")
        cols = ["unit", "system", "deprem sayısı", "slip_rate", "a (kısmi)", "b (kısmi)"] + \
            (["a (tamamen ayrı)"] if "a (tamamen ayrı)" in units else [])
        L.append(units.sort_values("deprem sayısı", ascending=False)[cols].head(15).round(3).to_markdown(index=False) + "\n")
        sysrows = [{"sistem": nm, "A (log çarpan)": round(s["A"][k], 3), "B (doluluk)": round(s["B"][k], 3)}
                   for k, nm in enumerate(W.systems)]
        L.append(pd.DataFrame(sysrows).to_markdown(index=False) + "\n")

    # büyük depremler
    if "partial" in res:
        L.append("## Büyük depremler (M≥7, 1900 sonrası): olaydan hemen önce verilen olasılık\n")
        m = res["partial"].merge(base[["i", "score", "p_cell", "rank_frac"]], on="i", suffixes=("", "_e"))
        big = m[(m["mw"] >= 7.0) & (m["year"] >= 1900)]
        t = pd.DataFrame({"tarih": big["time"].dt.date, "Mw": big["mw"].round(1), "sistem": big["system"],
                          "hücre olasılığı (model)": big["p_cell"].map(lambda x: f"%{100 * x:.2f}"),
                          "hücre olasılığı (ETAS)": big["p_cell_e"].map(lambda x: f"%{100 * x:.2f}"),
                          "en olası hücreler içindeki dilim (model)": big["rank_frac"].map(lambda x: f"%{100 * x:.0f}"),
                          "kazanç": (big["score"] - big["score_e"]).round(3)})
        L.append("Hücre olasılığı: bir sonraki M≥5,5 depremin o 0,5° hücrede olma olasılığı (546 hücre; eşit dağılım %0,18). "
                 "Dilim: kaç hücre bu hücreden daha olası görülmüş (%5 = en olası %5 içinde).\n")
        L.append(t.to_markdown(index=False) + "\n")

    # grafik
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 4.5))
        for n, col in [("partial", "#1f6f8b"), ("pooled", "#999999"), ("separate", "#c26b2e"), ("fixed", "#6b4c9a"),
                       ("no_interact", "#3c9d5d")]:
            if n not in res:
                continue
            m = res[n].merge(base[["i", "score"]], on="i", suffixes=("", "_e"))
            m = m[m["year"] >= 1900]
            ax.plot(m["time"], np.cumsum(m["score"] - m["score_e"]), label=labels.get(n, n), color=col, lw=1.4)
        ax.axhline(0, color="k", lw=0.6)
        ax.set_ylabel("ETAS'a göre birikimli kazanç (nat)")
        ax.set_title("Depremden depreme öğrenme: 1900'den itibaren birikimli puan farkı")
        ax.legend(fontsize=7, loc="best")
        fig.tight_layout()
        fig.savefig(FIG / "faultwalk_cumulative.png", dpi=130)
        L.append("![birikimli kazanç](figures/faultwalk_cumulative.png)\n")
    except Exception as ex:  # noqa: BLE001
        L.append(f"(grafik çizilemedi: {ex})\n")

    L.append("## Yöntem ve sınırlar\n")
    L.append("- Yöntem ayrıntısı: `src/deprem/faultwalk.py`, veri: `src/deprem/faultwalk_data.py`.")
    L.append("- ETAS burada aynı ızgara ve aynı sıralı güncellemeyle kuruldu; önceki raporlardaki ETAS'la birebir aynı kod değil, "
             "aynı formül. Böylece karşılaştırma aynı koşullarda.")
    L.append("- 1700-1899 dönemi eksik ve konumları ±30-50 km hatalı; bu dönemde model öğrenir ama ana sonuç 1900-2025'tir. "
             "Eksiklik, dönem başına öğrenilen tespit oranıyla hesaba katıldı; doğu ve deniz alanlarındaki eksiklik ayrıca modellenmedi.")
    L.append("- Fay doluluğunda saat, 1700'den beri kırılmamış faylarda 1700'de başlar (gerçek doluluk daha yüksek olabilir).")
    L.append("- Coulomb hesabında her depremin mekanizması en yakın GEM fayından varsayıldı; alıcı geometrisi hücreye en yakın büyük faydır.")
    L.append("- Büyüklük tahmini tüm modellerde aynı (Gutenberg-Richter); fark yalnız zaman ve konumda.")
    REP.write_text("\n".join(L), encoding="utf-8")
    print(REP)


if __name__ == "__main__":
    main()

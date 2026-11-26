"""Mekanik kırılma modeli sonuçlarını özetle: tablolar, şekiller ve rapor sayıları.

Kullanım: python scripts/mekanik_rapor.py
Girdi: runs/mekanik/oynatma_*.pkl, gps*.csv, ayar_denemeleri.csv
Çıktı: reports/tables/mekanik_*.csv, reports/figures/mekanik_*.png, runs/mekanik/ozet.json
"""
from __future__ import annotations

import json
import pickle
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from deprem.mekanik import deney, tarih  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "runs" / "mekanik"
TAB = ROOT / "reports" / "tables"
FIG = ROOT / "reports" / "figures"
XB = np.arange(0, 456, 1.0)
XC = XB[:-1] + 0.5


def auc(icer, dis):
    """İçerideki rastgele bir km'nin dışarıdakinden daha dolu olma olasılığı (0.5 = şans)."""
    icer, dis = np.asarray(icer), np.asarray(dis)
    icer, dis = icer[np.isfinite(icer)], dis[np.isfinite(dis)]
    if len(icer) == 0 or len(dis) == 0:
        return np.nan
    r = pd.Series(np.concatenate([icer, dis])).rank().to_numpy()
    return (r[:len(icer)].sum() - len(icer) * (len(icer) + 1) / 2) / (len(icer) * len(dis))


def km_dizisi(dol: pd.DataFrame):
    out = np.full(len(XC), np.nan)
    k = np.clip(np.floor(dol["x_km"].to_numpy()).astype(int), 0, len(XC) - 1)
    out[k] = dol["doluluk"].to_numpy()
    return out


def zaman_temelli(ev: pd.DataFrame, t: float, vpl_km: np.ndarray):
    """Basit kestirme: son kırılmadan beri geçen süre x kayma hızı (1500 öncesi bilinmiyor -> 1500)."""
    son = np.full(len(XC), 1500.0)
    for e in ev[ev.t < t].itertuples():
        sel = (XC >= e.x0) & (XC <= e.x1)
        son[sel] = np.maximum(son[sel], e.t)
    return (t - son) * vpl_km


def main():
    TAB.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    v = deney.Veri(ROOT)
    cells = v.cells
    ev_tab = tarih.tablo(cells)
    from deprem.mekanik import hazirlik
    vpl = hazirlik.levha_hizi(cells)
    k = np.clip(np.floor(cells["x_km"].to_numpy()).astype(int), 0, len(XC) - 1)
    vpl_km = np.bincount(k, vpl, len(XC)) / np.maximum(np.bincount(k, minlength=len(XC)), 1)

    uyeler = []
    for f in sorted(RUN.glob("oynatma_*.pkl")):
        if f.name.endswith(".parca.pkl"):
            continue
        uyeler += pickle.loads(f.read_bytes())["uyeler"]
    if not uyeler:
        for f in sorted(RUN.glob("oynatma_*.parca.pkl")):
            uyeler += pickle.loads(f.read_bytes())["uyeler"]
    n = len(uyeler)
    ozet = {"uye_sayisi": n, "senaryolar": sorted({u["senaryo"] for u in uyeler})}

    # 1) zorlanan tarihsel depremlerde modelin büyüklüğü
    z = pd.DataFrame([dict(r, uye=u["uye"], senaryo=u["senaryo"]) for u in uyeler for r in u["zorla"]])
    zt = z.groupby(["tarih"]).agg(mw_gozlem=("mw_gozlem", "first"), mw_model_ort=("mw_model", "mean"),
                                  mw_model_min=("mw_model", "min"), mw_model_maks=("mw_model", "max"),
                                  kayma_ort_m=("ort_kayma", "mean"), sahte_once=("sahte_once", "mean")).reset_index()
    zt.to_csv(TAB / "mekanik_tarihsel_buyukluk.csv", index=False)
    ozet["mw_fark_ort"] = float((zt.mw_model_ort - zt.mw_gozlem).mean())
    ozet["mw_fark_rms"] = float(np.sqrt(((zt.mw_model_ort - zt.mw_gozlem) ** 2).mean()))
    ozet["sahte_buyuk_deprem_ort"] = float(np.mean([u["sahte"] for u in uyeler]))
    ozet["sahte_buyuk_deprem_yil"] = float(np.mean([u["sahte"] for u in uyeler]) / (2026.9 - 1500))

    # 2) sınama: büyük depremden hemen önce kırılacak yer daha mı doluydu?
    testler = ev_tab[(ev_tab.mw >= deney.BUYUK_MW) & (ev_tab.t >= 1700)]
    rows = []
    for e in testler.itertuples():
        icer = (XC >= e.x0) & (XC <= e.x1)
        aucs = []
        for u in uyeler:
            d = km_dizisi(u["once"][e.tarih])
            aucs.append(auc(d[icer], d[~icer]))
        zb = zaman_temelli(ev_tab, e.t, vpl_km)
        rows.append({"tarih": e.tarih, "ad": e.ad, "mw": e.mw, "guven": e.guven,
                     "auc_model_ort": float(np.nanmean(aucs)), "auc_model_min": float(np.nanmin(aucs)),
                     "auc_model_maks": float(np.nanmax(aucs)),
                     "auc_zaman_kestirme": float(auc(zb[icer], zb[~icer]))})
    st = pd.DataFrame(rows)
    st.to_csv(TAB / "mekanik_sinama.csv", index=False)
    ozet["auc_model_ort"] = float(st.auc_model_ort.mean())
    ozet["auc_zaman_ort"] = float(st.auc_zaman_kestirme.mean())
    ozet["model_zamani_gecti"] = int((st.auc_model_ort > st.auc_zaman_kestirme).sum())
    ozet["test_sayisi"] = len(st)
    # şans sınaması: her testte kırık alanı fay boyunca rastgele kaydırılırsa AUC dağılımı
    rng = np.random.default_rng(0)
    perm = []
    for _ in range(300):
        vals = []
        for e in testler.itertuples():
            L = e.x1 - e.x0
            x0 = rng.uniform(0, 455 - L)
            icer = (XC >= x0) & (XC <= x0 + L)
            u = uyeler[rng.integers(n)]
            d = km_dizisi(u["once"][e.tarih])
            vals.append(auc(d[icer], d[~icer]))
        perm.append(np.nanmean(vals))
    ozet["auc_sans_p"] = float((np.array(perm) >= ozet["auc_model_ort"]).mean())
    ozet["auc_sans_ort"] = float(np.mean(perm))

    # 3) bugün: km başına doluluk ve 30 yıllık olasılık
    D = np.array([km_dizisi(u["bugun"]) for u in uyeler])
    P = np.array([u["ileri"]["p"] for u in uyeler])
    bug = pd.DataFrame({"x_km": XC, "doluluk_ort": np.nanmean(D, 0), "doluluk_p10": np.nanpercentile(D, 10, 0),
                        "doluluk_p90": np.nanpercentile(D, 90, 0), "p30_ort": P.mean(0)})
    lon = cells[cells.z_ust == 0].sort_values("x_km")
    bug["lon"] = np.interp(XC, lon.x_km, lon.lon)
    bug["lat"] = np.interp(XC, lon.x_km, lon.lat)
    bug.to_csv(TAB / "mekanik_bugun_km.csv", index=False)
    ilk = pd.concat([u["ileri"]["ilk"].assign(uye=u["uye"]) for u in uyeler if len(u["ileri"]["ilk"])],
                    ignore_index=True) if any(len(u["ileri"]["ilk"]) for u in uyeler) else pd.DataFrame()
    deniz = (bug.lon >= 27.0) & (bug.lon <= 29.3)  # Marmara Denizi içi (Ganos doğusu - Hersek)
    xd0, xd1 = float(bug.x_km[deniz].min()), float(bug.x_km[deniz].max())
    if len(ilk):
        ilk.to_csv(TAB / "mekanik_ileri_ilk_buyuk.csv", index=False)
    p_tum, p_deniz, p_deniz7 = [], [], []
    for u in uyeler:
        r, nk = u["ileri"]["ilk"], u["ileri"].get("n", 6)
        p_tum.append(len(r) / nk)
        if len(r):
            ov = (r.x1 >= xd0) & (r.x0 <= xd1)
            p_deniz.append(ov.sum() / nk)
            p_deniz7.append((ov & (r.mw >= 7.0)).sum() / nk)
        else:
            p_deniz.append(0.0)
            p_deniz7.append(0.0)
    ozet["p30_tum_fay"] = [float(np.mean(p_tum)), float(np.min(p_tum)), float(np.max(p_tum))]
    ozet["p30_marmara_denizi"] = [float(np.mean(p_deniz)), float(np.min(p_deniz)), float(np.max(p_deniz))]
    ozet["p30_marmara_denizi_m7"] = [float(np.mean(p_deniz7)), float(np.min(p_deniz7)), float(np.max(p_deniz7))]
    ozet["ileri_kosu_toplam"] = int(sum(u["ileri"].get("n", 6) for u in uyeler))
    if len(ilk):
        ozet["ileri_ilk_mw_medyan"] = float(ilk.mw.median())
    ozet["deniz_x"] = [xd0, xd1]
    json.dump(ozet, open(RUN / "ozet_ara.json", "w"), indent=1, default=float)

    # şekiller
    fig, ax = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True)
    ax[0].fill_between(XC, bug.doluluk_p10, bug.doluluk_p90, color="#9ecae1", label="topluluk %10-90")
    ax[0].plot(XC, bug.doluluk_ort, color="#08519c", lw=1.5, label="ortalama")
    ax[0].set_ylabel("doluluk (gerilme / dayanım)")
    ax[0].set_title("Bugün (2026) fay boyunca doluluk: 1 km çözünürlük")
    ax[0].legend(loc="upper right", fontsize=8)
    for e in ev_tab.itertuples():
        ax[1].plot([e.x0, e.x1], [e.t, e.t], lw=3 if e.mw >= 7 else 1.5, color="#d94801")
        ax[1].text(e.x1 + 2, e.t, e.tarih[:4], fontsize=7, va="center")
    ax[1].set_ylabel("yıl")
    ax[1].set_xlabel("fay boyunca uzaklık (km, batıdan doğuya: Saros → Düzce)")
    ax[1].set_title("Tarihsel kırılmalar (modelde zorlanan)")
    for xv, lab in [(0, "Saros"), (60, "Ganos"), (165, "Orta Marmara"), (255, "Adalar-Çınarcık"),
                    (300, "İzmit"), (410, "Düzce")]:
        ax[0].text(xv, ax[0].get_ylim()[1] * 0.95, lab, fontsize=7, color="gray")
    fig.tight_layout()
    fig.savefig(FIG / "mekanik_bugun_doluluk.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.errorbar(zt.mw_gozlem, zt.mw_model_ort, yerr=[zt.mw_model_ort - zt.mw_model_min,
                                                     zt.mw_model_maks - zt.mw_model_ort],
                fmt="o", color="#08519c")
    for r in zt.itertuples():
        ax.text(r.mw_gozlem + 0.02, r.mw_model_ort, r.tarih[:4], fontsize=7)
    lim = [5.8, 8.0]
    ax.plot(lim, lim, color="gray", lw=1)
    ax.set_xlabel("gözlenen Mw")
    ax.set_ylabel("modelin ürettiği Mw (aynı kırık alanında)")
    ax.set_title("Tarihsel depremler: biriken enerji gözlenen büyüklüğü veriyor mu?")
    fig.tight_layout()
    fig.savefig(FIG / "mekanik_tarihsel_mw.png", dpi=130)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(11, 3.2))
    ax.plot(XC, bug.p30_ort, color="#a50f15")
    ax.set_ylabel("30 yılda M≥6,8\nkırığına girme olasılığı")
    ax.set_xlabel("fay boyunca uzaklık (km)")
    ax.set_title("Modelin 30 yıllık ileri bakışı (topluluk ortalaması; doğrulanmış tahmin değildir)")
    fig.tight_layout()
    fig.savefig(FIG / "mekanik_30yil_olasilik.png", dpi=130)
    plt.close(fig)

    gps = sorted(RUN.glob("gps.csv"))
    if gps:
        g = pd.read_csv(gps[0]).sort_values("rms")
        best = g.iloc[0]
        s = pd.read_csv(RUN / f"gps_seri_{best.tM_ref:g}_{best.krip_carpan:g}.csv")
        fig, axs = plt.subplots(1, 2, figsize=(10, 3.5), sharey=True)
        for a, site in zip(axs, ["TUBI", "ISTA"]):
            for bil, c in (("e", "#08519c"), ("n", "#d94801")):
                q = s[(s.site == site) & (s.bil == bil)].sort_values("a")
                tm = (q.a + q.b) / 2
                a.plot(tm, q.anomali_g, "o", color=c, label=f"GPS {bil}")
                a.plot(tm, q.anomali_m, "-", color=c, label=f"model {bil}")
            a.set_title(f"{site}: 1999 sonrası hız fazlası")
            a.set_xlabel("yıl")
        axs[0].set_ylabel("mm/yıl (2015-2023 hızına göre)")
        axs[0].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(FIG / "mekanik_gps_1999.png", dpi=130)
        plt.close(fig)
        ozet["gps_en_iyi"] = {"tM_ref": float(best.tM_ref), "krip_carpan": float(best.krip_carpan),
                              "rms": float(best.rms), "rms_en_kotu": float(g.rms.max())}
    json.dump(ozet, open(RUN / "ozet.json", "w"), indent=1, default=float)
    print(json.dumps(ozet, indent=1, default=float))
    print(zt.round(2).to_string())
    print(st.round(3).to_string())


if __name__ == "__main__":
    main()

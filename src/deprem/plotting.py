"""Statik rapor grafikleri (matplotlib, harita kütüphanesi gerektirmez).

Renk kuralları: büyüklük/yoğunluk tek tonlu sıralı skala; b-değeri ve
moment oranı nötr orta noktalı ıraksak skala; kategoriler sabit sırada.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import requests  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, LogNorm, TwoSlopeNorm  # noqa: E402

from .region import Region  # noqa: E402

INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"
SURFACE = "#fcfcfb"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#e8f1fb", "#9cc3ee", "#2a78d6", "#0d3f7f"])
DIV = LinearSegmentedColormap.from_list("div", ["#2a78d6", "#d9d8d3", "#e34948"])

COAST_URL = ("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/"
             "geojson/ne_50m_coastline.geojson")

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
    "axes.spines.right": False, "font.size": 9, "axes.titlesize": 11, "axes.titleweight": "bold",
    "axes.titlecolor": INK, "legend.frameon": False,
})


def load_coastlines(cache: Path, region: Region):
    cache = Path(cache)
    if not cache.exists():
        try:
            cache.parent.mkdir(parents=True, exist_ok=True)
            r = requests.get(COAST_URL, timeout=120)
            r.raise_for_status()
            cache.write_bytes(r.content)
        except Exception:
            return []
    lines = []
    for f in json.loads(cache.read_text())["features"]:
        g = f["geometry"]
        parts = [g["coordinates"]] if g["type"] == "LineString" else g["coordinates"]
        for p in parts:
            c = np.asarray(p)
            if ((c[:, 0] > region.min_lon - 1) & (c[:, 0] < region.max_lon + 1)
                    & (c[:, 1] > region.min_lat - 1) & (c[:, 1] < region.max_lat + 1)).any():
                lines.append(c)
    return lines


def _basemap(ax, region: Region, coast=None, faults=None, fault_color="#8a8984"):
    for c in coast or []:
        ax.plot(c[:, 0], c[:, 1], color=INK2, lw=0.5, zorder=1)
    if faults is not None:
        for coords in faults["coords"]:
            c = np.asarray(coords)
            ax.plot(c[:, 0], c[:, 1], color=fault_color, lw=0.7, zorder=2)
    ax.set_xlim(region.min_lon, region.max_lon)
    ax.set_ylim(region.min_lat, region.max_lat)
    ax.set_aspect(1 / np.cos(np.radians((region.min_lat + region.max_lat) / 2)))
    ax.set_xlabel("Boylam (°D)")
    ax.set_ylabel("Enlem (°K)")
    ax.grid(False)


def epicenter_map(df, region, out, coast=None, faults=None, min_mag=None, title=None):
    d = df if min_mag is None else df[df["mag"] >= min_mag]
    d = d.sort_values("mag")
    fig, ax = plt.subplots(figsize=(11, 5.2))
    _basemap(ax, region, coast, faults)
    sc = ax.scatter(d["longitude"], d["latitude"], s=1.5 * 2.2 ** (d["mag"] - 2),
                    c=d["depth_km"], cmap=SEQ, vmin=0, vmax=40, lw=0, alpha=0.7, zorder=3)
    big = d[d["mag"] >= 6.5]
    ax.scatter(big["longitude"], big["latitude"], s=60, facecolor="none",
               edgecolor=INK, lw=0.8, zorder=4, label="M ≥ 6,5")
    fig.colorbar(sc, ax=ax, label="Derinlik (km)", shrink=0.8)
    if len(big):
        ax.legend(loc="lower left")
    ax.set_title(title or f"Episantrlar ({len(d):,} olay)")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def density_map(h, lats, lons, region, out, coast=None, faults=None, title="", label=""):
    fig, ax = plt.subplots(figsize=(11, 5.2))
    hm = np.where(h > 0, h, np.nan)
    pc = ax.pcolormesh(lons, lats, hm, cmap=SEQ, norm=LogNorm(vmin=1, vmax=np.nanmax(hm)),
                       zorder=0)
    _basemap(ax, region, coast, faults)
    fig.colorbar(pc, ax=ax, label=label, shrink=0.8)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def b_map(b, lats, lons, region, out, coast=None, faults=None, title="", vcenter=None):
    """Iraksak skala, bölgesel medyan b etrafında: göreli düşük/yüksek b gösterilir."""
    vc = float(np.nanmedian(b)) if vcenter is None else vcenter
    fig, ax = plt.subplots(figsize=(11, 5.2))
    pc = ax.pcolormesh(lons, lats, b, cmap=DIV.reversed(),
                       norm=TwoSlopeNorm(vcenter=vc, vmin=vc - 0.4, vmax=vc + 0.4), zorder=0)
    _basemap(ax, region, coast, faults, fault_color=INK)
    fig.colorbar(pc, ax=ax, label=f"b-değeri (orta = medyan {vc:.2f}; kırmızı: göreli düşük b)",
                 shrink=0.8)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fmd_plot(m, fit, out, dm=0.1, title=""):
    from .stats.gr import fmd
    bins, inc, cum = fmd(m, dm)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.semilogy(bins, cum, "o", ms=4, color=SERIES[0], label="Birikimli N(≥M)")
    ax.semilogy(bins, np.where(inc > 0, inc, np.nan), "s", ms=3, color=SERIES[1],
                label="Artımlı N(M)")
    mm = np.linspace(fit.mc, bins.max(), 50)
    ax.semilogy(mm, fit.n * 10 ** (-fit.b * (mm - fit.mc)), "-", lw=2, color=INK,
                label=f"GR uyumu: b = {fit.b:.2f} ± {fit.sigma_b:.2f}")
    ax.axvline(fit.mc, color=INK2, lw=1, ls="--")
    ax.annotate(f"Mc = {fit.mc:.1f}", (fit.mc, cum.max()), xytext=(5, -10),
                textcoords="offset points", color=INK2)
    ax.set_xlabel("Büyüklük")
    ax.set_ylabel("Olay sayısı")
    ax.set_title(title)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def magnitude_time(df, out, mc=None, title=""):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 6), sharex=True,
                                   gridspec_kw={"height_ratios": [2, 1]})
    ax1.scatter(df["time"], df["mag"], s=0.8 * 2 ** (df["mag"] - 2), color=SERIES[0],
                lw=0, alpha=0.5)
    big = df[df["mag"] >= 6.5]
    for _, r in big.iterrows():
        ax1.annotate(f"{r['time'].year} M{r['mag']:.1f}", (r["time"], r["mag"]),
                     xytext=(3, 3), textcoords="offset points", fontsize=7, color=INK2)
    ax1.set_ylabel("Büyüklük")
    ax1.set_title(title)
    d = df if mc is None else df[df["mag"] >= mc]
    s = d.set_index("time").resample("MS").size()
    ax2.plot(s.index, s.values, lw=1.2, color=SERIES[0])
    ax2.set_yscale("log")
    ax2.set_ylabel(f"Aylık sayı{'' if mc is None else f' (M≥{mc:.1f})'}")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def cumulative_plot(series: dict, out, title="", ylabel="Birikimli olay sayısı"):
    fig, ax = plt.subplots(figsize=(11, 4.5))
    for k, (name, s) in enumerate(series.items()):
        ax.plot(s.index, s.values, lw=2, color=SERIES[k % len(SERIES)], label=name)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def b_time_plot(bt: pd.DataFrame, out, title=""):
    fig, ax = plt.subplots(figsize=(11, 4))
    ax.fill_between(bt["time"], bt["b"] - bt["sigma_b"], bt["b"] + bt["sigma_b"],
                    color=SERIES[0], alpha=0.2, lw=0)
    ax.plot(bt["time"], bt["b"], lw=2, color=SERIES[0])
    ax.set_ylabel("b-değeri (±1σ)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def fault_rate_map(faults, region, out, coast=None, title=""):
    fig, ax = plt.subplots(figsize=(11, 5.2))
    _basemap(ax, region, coast, None)
    sr = faults["slip_rate_mm_yr"].fillna(0)
    norm = LogNorm(vmin=0.1, vmax=30)
    for coords, v in zip(faults["coords"], sr):
        c = np.asarray(coords)
        ax.plot(c[:, 0], c[:, 1], color=SEQ(norm(max(v, 0.1))), lw=0.8 + 2.2 * norm(max(v, 0.1)),
                solid_capstyle="round", zorder=3)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=SEQ)
    fig.colorbar(sm, ax=ax, label="Kayma hızı (mm/yıl, EMME)", shrink=0.8)
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)


def moment_budget_plot(zb: pd.DataFrame, out, title=""):
    fig, ax = plt.subplots(figsize=(9, 4.2))
    y = np.arange(len(zb))
    acc = zb["m0_accumulated"].to_numpy()
    lo = zb["m0_acc_min"].to_numpy()
    hi = zb["m0_acc_max"].to_numpy()
    rel = zb["m0_released"].to_numpy()
    ax.hlines(y, lo, hi, color=SERIES[0], lw=6, alpha=0.3, label="Birikim (kayma hızı min-maks)")
    ax.plot(acc, y, "o", color=SERIES[0], ms=8, label="Birikim (en iyi kayma hızı)")
    ax.plot(rel, y, "D", color=SERIES[1], ms=8, label="Katalogda boşalan moment")
    ax.set_xscale("log")
    ax.set_yticks(y, zb["zone"])
    ax.invert_yaxis()
    ax.set_xlabel("Sismik moment (N·m)")
    ax.set_title(title)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)

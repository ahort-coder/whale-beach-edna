"""
Sediment core figure: lithology (data-derived), sand/silt/clay, mean grain size, LOI.
Formatted for Nature Geoscience (183 mm double-column width, 5-7 pt sans-serif, vector PDF).

Inputs (on Desktop): GSData/*.csv (laser grain-size, cumulative %Pass), loi_data.csv
Grain-size classes (Wentworth): sand > 62.5 um, silt 3.9-62.5 um, clay < 3.9 um.
Mean grain size: Folk & Ward (1957) graphic mean, Mz = (phi16 + phi50 + phi84)/3, reported in um.
"""
import os, re, glob
import numpy as np, pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D

DESK = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
GSDIR = os.path.join(DESK, "GSData")
LOIFILE = os.path.join(DESK, "loi_data.csv")
SKIP = {"WB2509B26_30.csv"}          # duplicate run; WB2509B_26_30.csv used instead
CORES = ["WB2505B", "WB2507B", "WB2508B", "WB2509B"]
MAXD = 50

# ---------------- data ----------------
def load_gs(f):
    t = open(f, encoding="utf-8-sig").read().replace("\r\n", "\n").replace("\r", "\n")
    rows = [r.split(",") for r in t.strip().split("\n")][1:]
    s, p = [], []
    for r in rows:
        for k in (0, 3):                       # two column blocks per row
            if len(r) > k + 2 and r[k].strip():
                s.append(float(r[k])); p.append(float(r[k + 2]))
    o = np.argsort(s)
    return np.array(s)[o], np.array(p)[o]      # size ascending, % finer

def gs_metrics(size, pfiner):
    ls = np.log(size)
    at = lambda x: float(np.interp(np.log(x), ls, pfiner))
    sand = 100 - at(62.5); clay = at(3.9); silt = 100 - sand - clay
    pp, i = np.unique(pfiner, return_index=True)
    d = lambda pf: float(np.exp(np.interp(pf, pp, ls[i])))
    phi = lambda x: -np.log2(x / 1000)
    mz = (phi(d(84)) + phi(d(50)) + phi(d(16))) / 3   # phi16 = 16% coarser = 84% finer
    return sand, silt, clay, 1000 * 2 ** -mz

rows = []
for f in sorted(glob.glob(os.path.join(GSDIR, "*.csv"))):
    b = os.path.basename(f)
    if b in SKIP: continue
    m = re.match(r"(WB\d{4})B_([\d.]+)_([\d.]+)\.csv", b)
    sand, silt, clay, mz = gs_metrics(*load_gs(f))
    rows.append(dict(core=m.group(1) + "B", top=float(m.group(2)), bot=float(m.group(3)),
                     sand=sand, silt=silt, clay=clay, mean_um=mz))
gs = pd.DataFrame(rows).sort_values(["core", "top"]).reset_index(drop=True)

L = pd.read_csv(LOIFILE, encoding="utf-8-sig").iloc[:, :5]
L.columns = ["sample", "dry", "ash", "crucible", "LOI"]
L["core"] = L["sample"].str.split().str[0].str.rstrip("B") + "B"  # WB2507 -> WB2507B
L[["top", "bot"]] = L["sample"].str.split().str[1].str.split("-", expand=True).astype(float)
loi = L[["core", "top", "bot", "LOI"]].sort_values(["core", "top"]).reset_index(drop=True)

# ------------- data-derived lithology (on LOI intervals, the finer sampling) -------------
def sand_for(core, top, bot):
    g = gs[(gs.core == core) & (gs.top < bot) & (gs.bot > top)]
    if g.empty: return np.nan
    w = np.minimum(g.bot, bot) - np.maximum(g.top, top)
    return float(np.average(g.sand, weights=w))

def classify(l, sand):
    if (not np.isnan(sand) and sand >= 50 and l < 35) or (np.isnan(sand) and l < 10):
        return "sand"
    if l >= 40: return "peat"
    if l >= 20: return "orgmud"
    return "minmud"

loi["sand"] = [sand_for(r.core, r.top, r.bot) for r in loi.itertuples()]
loi["litho"] = [classify(r.LOI, r.sand) for r in loi.itertuples()]
loi["sandy"] = loi.sand.between(25, 50, inclusive="left") | ((loi.sand >= 50) & (loi.litho != "sand"))

# ---------------- style ----------------
MM = 1 / 25.4
mpl.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 6, "axes.labelsize": 6, "xtick.labelsize": 5, "ytick.labelsize": 5,
    "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2, "xtick.minor.size": 1.2, "xtick.minor.width": 0.4,
    "xtick.direction": "out", "ytick.direction": "out",
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "hatch.linewidth": 0.35, "hatch.color": "#3d3528",
})
LITHO = {  # fill, hatch, legend label
    "peat":   ("#5b4030", None, "Peat / organic-rich mud (LOI ≥ 40%)"),
    "orgmud": ("#9a7a5c", None, "Organic mud (LOI 20–40%)"),
    "minmud": ("#c9bfae", None, "Minerogenic mud (LOI < 20%)"),
    "sand":   ("#ecd391", "....", "Sand (> 50% sand)"),
}
TEX = {"sand": "#e3c46f", "silt": "#a9a090", "clay": "#4f5a66"}
INK, GRID, BAND = "#1a1a1a", "#d9d9d9", "#f6ecc9"

# Event-deposit attribution (top, bottom in cm; confident=True -> filled symbol, False -> open/tentative)
EVENTS = {
    "WB2505B": [("noreaster", 23, 26, True),   ("h1821", 42, 50, True)],
    "WB2507B": [("noreaster", 22, 25, True),   ("h1821", 48, 50, False)],  # 48.6% sand, LOI 4.7%
    "WB2508B": [("noreaster", 24.5, 27, False), ("h1821", 47, 50, False)], # no GS at 24.5-27; 46% sand at base
    "WB2509B": [("noreaster", 24, 25, True),   ("h1821", 42, 50, True)],
}
EVSYM = {"noreaster": dict(marker="D", ms=4.2, label="1962 nor'easter overwash"),
         "h1821":     dict(marker="*", ms=7.5, label="1821 hurricane overwash")}

fig = plt.figure(figsize=(183 * MM, 112 * MM))
left, right, bottom, top = 0.055, 0.995, 0.24, 0.93
gap = 0.022
gw = (right - left - 3 * gap) / 4
wr = np.array([0.42, 1, 1, 1]); inner = 0.006
unit = (gw - 3 * inner) / wr.sum()

def step_profile(ax, df, col, **kw):
    """Interval-true step profile; breaks where intervals are not contiguous."""
    segs, cur = [], []
    for r in df.itertuples():
        if cur and abs(cur[-1][1] - r.top) > 1e-6:
            segs.append(cur); cur = []
        cur += [(getattr(r, col), r.top), (getattr(r, col), r.bot)]
    if cur: segs.append(cur)
    for s in segs:
        x, y = zip(*s); ax.plot(x, y, color=INK, lw=0.8, solid_joinstyle="miter", **kw)
    ax.plot(df[col], (df.top + df.bot) / 2, "o", ms=1.8, mfc=INK, mec="none")

for ci, core in enumerate(CORES):
    x0 = left + ci * (gw + gap)
    axes, x = [], x0
    for w in wr:
        ax = fig.add_axes([x, bottom, unit * w, top - bottom]); axes.append(ax); x += unit * w + inner
    aL, aT, aM, aO = axes
    g, l = gs[gs.core == core], loi[loi.core == core]

    # sand-layer bands across data panels
    for r in l[l.litho == "sand"].itertuples():
        for ax in (aM, aO):
            ax.axhspan(r.top, r.bot, color=BAND, lw=0, zorder=0)

    # lithology column
    for r in l.itertuples():
        fc, h, _ = LITHO[r.litho]
        aL.add_patch(Rectangle((0, r.top), 1, r.bot - r.top, fc=fc, ec="white", lw=0.3))
        hh = h if h else ("..." if r.sandy else None)
        if hh:
            hc = "#f3ead8" if r.litho == "peat" else "#3d3528"
            aL.add_patch(Rectangle((0, r.top), 1, r.bot - r.top, fill=False, ec=hc, lw=0, hatch=hh))
    aL.add_patch(Rectangle((0, 0), 1, l.bot.max(), fill=False, ec=INK, lw=0.6))
    for ev, t0, t1, ok in EVENTS.get(core, []):
        aL.plot(0.5, (t0 + t1) / 2, ls="none", marker=EVSYM[ev]["marker"], ms=EVSYM[ev]["ms"],
                mfc=INK if ok else "white", mec=INK if not ok else "white", mew=0.6 if not ok else 0.4,
                zorder=5, clip_on=False)
    aL.set_xlim(0, 1); aL.set_xticks([])
    for s in ("top", "right", "bottom"): aL.spines[s].set_visible(False)
    aL.spines["left"].set_visible(False)
    aL.tick_params(axis="y", length=2, pad=1)

    # texture: stacked horizontal bars
    for r in g.itertuples():
        x = 0
        for k in ("sand", "silt", "clay"):
            v = getattr(r, k)
            aT.add_patch(Rectangle((x, r.top), v, r.bot - r.top, fc=TEX[k], ec="white", lw=0.3))
            x += v
    aT.set_xlim(0, 100); aT.set_xticks([0, 50, 100])

    # mean grain size
    step_profile(aM, g, "mean_um")
    aM.set_xscale("log"); aM.set_xlim(5, 400)
    aM.set_xticks([10, 100]); aM.set_xticklabels(["10", "100"])
    aM.axvline(62.5, color="#8a8a8a", lw=0.4, ls=(0, (2, 1.5)), zorder=1)

    # LOI
    step_profile(aO, l, "LOI")
    aO.set_xlim(0, 70); aO.set_xticks([0, 30, 60])

    for ax in axes:
        ax.set_ylim(MAXD, 0)
        ax.set_yticks(range(0, MAXD + 1, 10))
        ax.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(5))
    for ax in (aT, aM, aO):
        ax.tick_params(axis="y", which="both", left=False, labelleft=False)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="x", pad=1.5)
        ax.xaxis.grid(True, which="major", color=GRID, lw=0.35); ax.set_axisbelow(True)
    aT.set_xlabel("Texture\n(%)", labelpad=2)
    aM.set_xlabel("Mean grain\nsize (µm)", labelpad=2)
    aO.set_xlabel("LOI\n(%)", labelpad=2)
    if ci == 0:
        aL.set_ylabel("Depth (cm)", labelpad=2)

    fig.text(x0, 0.965, "abcd"[ci], fontsize=8, fontweight="bold", va="bottom", ha="left")
    fig.text(x0 + 0.012, 0.965, core, fontsize=6.5, va="bottom", ha="left")

# legend
lh = [Patch(fc=LITHO[k][0], hatch=LITHO[k][1], ec="#7a7a7a", lw=0.3, label=LITHO[k][2])
      for k in ("peat", "orgmud", "minmud", "sand")]
lh.append(Patch(fc="#c9bfae", hatch="...", ec="#7a7a7a", lw=0.3, label="Sandy (25–50% sand)"))
lh.append(Patch(fc=BAND, ec="none", label="Sand layer (shaded on profiles)"))
th = [Patch(fc=TEX[k], ec="none", label=t) for k, t in
      (("sand", "Sand (> 62.5 µm)"), ("silt", "Silt (3.9–62.5 µm)"), ("clay", "Clay (< 3.9 µm)"))]
kw = dict(frameon=False, fontsize=5.5, handlelength=1.6, handleheight=1.0, columnspacing=1.2,
          handletextpad=0.5, borderaxespad=0)
leg1 = fig.legend(handles=lh, loc="upper left", bbox_to_anchor=(left, 0.135), ncol=3, title="Lithology",
                  title_fontsize=6, alignment="left", **kw)
fig.legend(handles=th, loc="upper left", bbox_to_anchor=(0.63, 0.135), ncol=1, title="Texture",
           title_fontsize=6, alignment="left", **kw)

eh = [Line2D([], [], ls="none", marker=v["marker"], ms=v["ms"], mfc=INK, mec="white", mew=0.4, label=v["label"])
      for v in EVSYM.values()]
eh.append(Line2D([], [], ls="none", marker="o", ms=3.5, mfc="white", mec=INK, mew=0.6, label="Open symbol: tentative"))
fig.legend(handles=eh, loc="upper left", bbox_to_anchor=(0.80, 0.135), ncol=1, title="Event deposits",
           title_fontsize=6, alignment="left", **kw)

out = os.path.dirname(os.path.abspath(__file__))
for ext, extra in (("pdf", {}), ("svg", {}), ("png", {"dpi": 600}), ("tiff", {"dpi": 600, "pil_kwargs": {"compression": "tiff_lzw"}})):
    fig.savefig(os.path.join(out, f"WB_core_figure.{ext}"), **extra)

gs.round(2).to_csv(os.path.join(out, "grain_size_summary.csv"), index=False)
loi.round(2).to_csv(os.path.join(out, "loi_lithology_summary.csv"), index=False)
print("done")

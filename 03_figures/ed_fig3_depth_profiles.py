"""Extended Data Fig. 2: downcore class composition of each eDNA core, deposit types marked.

Per-sample relative abundance (% of reads in each sample). The ten most abundant
assigned classes (mean per-sample relative abundance across all 64 samples) are
shown; all other classes and ASVs unassigned at class rank are grouped as "Others".
Deposit types follow Supplementary Table 2 (depth of the top of each sampled interval).
"""
import sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from load import load

out = sys.argv[1] if len(sys.argv) > 1 else "depth_profiles.png"
CORES = ["WB2501", "WB2506", "WB2508", "WB2509"]
NE = {"WB2501": (20, 26), "WB2506": (24, 27), "WB2508": (20, 30), "WB2509": (20, 24)}
def dtype(core, top):
    lo, hi = NE[core]
    if lo <= top <= hi: return "Nor’easter overwash"
    if top >= 40: return "Hurricane overwash"
    return "Salt-marsh peat"

rel = {}
for c in CORES:
    X, T = load(c)
    g = X.groupby(T["Class"].fillna("Unassigned")).sum()
    r = g / g.sum(0) * 100
    r = r[sorted(r.columns, key=lambda s: float(s.split("-")[0]))]
    rel[c] = r
allm = pd.concat([rel[c] for c in CORES], axis=1).fillna(0)
TOP = [k for k in allm.mean(1).sort_values(ascending=False).index if k != "Unassigned"][:10]

COL = {"Nanoarchaeia": "#2A7AD5", "Thermoplasmata": "#1CB17B", "Bathyarchaeia": "#EDA200",
       "Anaerolineae": "#EC6836", "Dehalococcoidia": "#E87BA4", "Microgenomatia": "#008000",
       "Bacteroidia": "#6A3D9A", "Lokiarchaeia": "#8C564B", "Clostridia": "#17BECF",
       "Gammaproteobacteria": "#B2182B", "Others": "#BFBEB9"}
spare = iter(["#4B3BA8", "#666600", "#FF7F00"])
for k in TOP:
    if k not in COL: COL[k] = next(spare)
DT = {"Salt-marsh peat": dict(facecolor="#D8C7A0", hatch=None),
      "Nor’easter overwash": dict(facecolor="#2B2B2B", hatch=None),
      "Hurricane overwash": dict(facecolor="white", hatch="xxxx")}

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "hatch.linewidth": 0.6,
                     "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6})
fig = plt.figure(figsize=(7.2, 8.4), dpi=300)
gs = fig.add_gridspec(2, 2, left=0.08, right=0.985, top=0.965, bottom=0.17, hspace=0.18, wspace=0.22)
tables = {}
for n, c in enumerate(CORES):
    sub = gs[n // 2, n % 2].subgridspec(1, 2, width_ratios=[0.045, 1], wspace=0.03)
    axs, ax = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
    r = rel[c]
    tab = r.reindex(TOP).fillna(0).T
    tab["Others"] = 100 - tab.sum(1)
    tables[c] = tab
    tops = [float(s.split("-")[0]) for s in r.columns]
    labels = ["–".join(str(int(float(v))) for v in s.split("-")) for s in r.columns]
    y = np.arange(len(labels))
    left = np.zeros(len(y))
    for k in TOP + ["Others"]:
        ax.barh(y, tab[k].values, left=left, height=0.86, color=COL[k], edgecolor="white", linewidth=0.3)
        left += tab[k].values
    types = [dtype(c, t) for t in tops]
    for yi, t in zip(y, types):
        axs.bar(0, 0.86, bottom=yi - 0.43, width=1, edgecolor="black", linewidth=0.4, **DT[t])
    for yi in range(1, len(types)):
        if types[yi] != types[yi - 1]:
            ax.axhline(yi - 0.5, color="black", lw=0.8, ls=(0, (3, 1.5)))
    for a in (ax, axs):
        a.set_ylim(len(y) - 0.5, -0.5)
    axs.set_xlim(-0.5, 0.5); axs.set_xticks([])
    axs.set_yticks(y, labels); axs.tick_params(axis="y", length=2, pad=1.5)
    for s in ("top", "right", "bottom"): axs.spines[s].set_visible(False)
    ax.set_yticks([]); ax.set_xlim(0, 100); ax.set_xticks([0, 25, 50, 75, 100])
    for s in ("top", "right", "left"): ax.spines[s].set_visible(False)
    ax.set_xlabel("Relative abundance (%)")
    if n % 2 == 0: axs.set_ylabel("Depth (cm)")
    ax.set_title(c, fontsize=8, loc="left", pad=3)
    axs.text(-0.5, 1.012, "abcd"[n], transform=axs.transAxes, fontsize=10, fontweight="bold",
             va="bottom", ha="right")

h1 = [Patch(facecolor=COL[k], label=k) for k in TOP + ["Others"]]
h2 = [Patch(edgecolor="black", linewidth=0.4, label=k, **v) for k, v in DT.items()]
fig.legend(handles=h1, title="Class", loc="lower left", bbox_to_anchor=(0.06, 0.005), ncol=4,
           frameon=False, fontsize=7, title_fontsize=7.5, alignment="left", handlelength=1.4, columnspacing=1.2)
fig.legend(handles=h2, title="Deposit type", loc="lower left", bbox_to_anchor=(0.75, 0.005), ncol=1,
           frameon=False, fontsize=7, title_fontsize=7.5, alignment="left", handlelength=1.4)
fig.savefig(out, dpi=300)
pd.concat(tables, names=["core", "depth_cm"]).round(2).to_csv(out.replace(".png", "_values.csv"))
print("top classes:", TOP)

"""Figure S1 (class level): top five classes per core, pooled relative abundance.

Same method as phylum_composition.png: for each core, reads are summed across
all samples, each class's share of the core's total reads is computed, the top
five assigned classes per core are shown, and everything else (other classes
plus ASVs unassigned at class rank) is grouped as "Others".
"""
import sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from load import load

CORES = ["WB2501", "WB2506", "WB2508", "WB2509"]
TOPN = 5
out = sys.argv[1] if len(sys.argv) > 1 else "class_composition.png"

pct = {}
for c in CORES:
    X, T = load(c)
    g = X.groupby(T["Class"].fillna("NA")).sum().sum(1)
    pct[c] = g / g.sum() * 100

top = {c: [k for k in pct[c].sort_values(ascending=False).index if k != "NA"][:TOPN] for c in CORES}
classes = sorted({k for v in top.values() for k in v},
                 key=lambda k: -np.mean([pct[c].get(k, 0) for c in CORES]))
tab = pd.DataFrame({c: {k: (pct[c][k] if k in top[c] else 0.0) for k in classes} for c in CORES}).T
tab["Others"] = 100 - tab.sum(1)
tab.round(2).to_csv(out.replace(".png", "_values.csv"))

# colours follow the parent phylum's colour in the phylum-level figure
COL = {"Nanoarchaeia": "#2A7AD5",     # Nanoarchaeota
       "Anaerolineae": "#EC6836",     # Chloroflexi
       "Dehalococcoidia": "#E87BA4",  # Chloroflexi (second class)
       "Thermoplasmata": "#1CB17B",   # Thermoplasmatota
       "Bathyarchaeia": "#EDA200",    # Crenarchaeota
       "Microgenomatia": "#008000",   # Patescibacteria
       "Others": "#BFBEB9"}
DARK = {"Nanoarchaeia", "Anaerolineae", "Microgenomatia"}
extra = ["#4B3BA8", "#8C564B", "#17BECF"]
for k in classes:
    if k not in COL:
        COL[k] = extra.pop(0); DARK.add(k)

plt.rcParams.update({"font.size": 13, "font.family": "DejaVu Sans"})
fig, ax = plt.subplots(figsize=(8, 6.5), dpi=300)
x = np.arange(len(CORES)); w = 0.6
bottom = np.zeros(len(CORES))
for k in classes + ["Others"]:
    v = tab[k].values
    ax.bar(x, v, w, bottom=bottom, color=COL[k], edgecolor="white", linewidth=1.5, label=k)
    for i, val in enumerate(v):
        if val >= 2.5:
            ax.text(x[i], bottom[i] + val / 2, f"{val:.1f}",
                    ha="center", va="center", fontsize=11,
                    color="white" if k in DARK else "black")
    bottom += v
ax.set_xticks(x, CORES)
ax.set_xlim(-0.5, len(CORES) - 0.5)
ax.set_ylim(0, 100)
ax.set_ylabel("Relative abundance (%)"); ax.set_xlabel("Core")
for s in ("top", "right"): ax.spines[s].set_visible(False)
for s in ("left", "bottom"): ax.spines[s].set_color("#808080")
h, l = ax.get_legend_handles_labels()
ax.legend(h[::-1], l[::-1], title="Class", frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.0))
fig.subplots_adjust(left=0.10, right=0.69, top=0.97, bottom=0.11)
fig.savefig(out, dpi=300)
print(tab.round(1))

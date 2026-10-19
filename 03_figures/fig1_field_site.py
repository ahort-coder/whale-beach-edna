"""Field-site figure (Nature Geoscience style, 180 mm wide): Whale Beach, NJ.
a  contiguous U.S.        b  New Jersey        c  Strathmere / Whale Beach      d  core locations
Core colours: blue = eDNA only, yellow = sediment only, green = both, black = neither.
Basemaps: Natural Earth 10 m (public domain); panels c-d are vector graphics classified from USGS National Map imagery (public domain).
"""
import csv, json, math, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import Rectangle, FancyArrowPatch, Polygon
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
from pyproj import Transformer
from PIL import Image
from scipy import ndimage as ndi

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
CSV = os.path.expanduser("~/Desktop/Annabelle DNA NJ Sites.csv")

plt.rcParams.update({"font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 6,
                     "axes.linewidth": 0.5, "pdf.fonttype": 42, "ps.fonttype": 42})
MM = 1 / 25.4
HALO = [pe.withStroke(linewidth=1.6, foreground="white")]
DHALO = [pe.withStroke(linewidth=1.6, foreground="black")]

# ---------------------------------------------------------------- core data
EDNA = {"WB2501", "WB2506", "WB2508", "WB2509"}
SED = {"WB2505", "WB2507", "WB2508", "WB2509"}
COL = {"eDNA": "#0057e7", "sed": "#ffd400", "both": "#00c040", "none": "#000000"}
import re
cores = []
# the Name field contains an unquoted comma ("Jun 2, 2025"), so parse each record with a regex
for nm, la, lo in re.findall(r"(WB\d{4}) [^,\n]*,[^,\n]*? at \d\d:\d\d,(-?\d+\.\d+),(-?\d+\.\d+)", open(CSV, encoding="utf-8-sig").read()):
    e, s_ = nm in EDNA, nm in SED
    cores.append(dict(name=nm, lat=float(la), lon=float(lo), kind="both" if e and s_ else "eDNA" if e else "sed" if s_ else "none"))
assert len(cores) == 9, cores
cores.sort(key=lambda c: c["name"])
GAUGE = dict(name="Strathmere (NOAA 8535163)", lat=39.2000, lon=-74.6567)
SITE = (np.mean([c["lon"] for c in cores]), np.mean([c["lat"] for c in cores]))
SEAWARD = 124.0  # compass bearing (deg) perpendicular to the shoreline at Whale Beach, ocean to the ESE

# ---------------------------------------------------------------- helpers
def read_geoms(path, filt=None):
    rd = shpreader.Reader(path)
    return [rec.geometry for rec in rd.records() if filt is None or filt(rec.attributes)]
land = read_geoms(os.path.join(DATA, "ne_10m_land", "ne_10m_land.shp"))
states = [(rec.attributes["postal"], rec.geometry) for rec in
          shpreader.Reader(os.path.join(DATA, "ne_10m_admin_1", "ne_10m_admin_1_states_provinces.shp")).records()
          if rec.attributes["adm0_a3"] == "USA"]

OCEAN, LAND, STATEC, NJC = "#dbe9f4", "#f0eee8", "#9a9a9a", "#d9d9d9"

def label(ax, x, y, s, halo=True, **kw):
    t = ax.text(x, y, s, **kw)
    if halo: t.set_path_effects(HALO)
    return t

def north_arrow(ax, x=0.93, y=0.80, h=0.14, color="k", halo=True):
    a = ax.annotate("", xy=(x, y + h), xytext=(x, y), xycoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=0.9, mutation_scale=7, shrinkA=0, shrinkB=0))
    t = ax.text(x, y - 0.015, "N", transform=ax.transAxes, ha="center", va="top", fontsize=7, fontweight="bold", color=color)
    if halo: t.set_path_effects(HALO if color == "k" else DHALO)

def scale_bar(ax, length_m, unit_m_per_data, x0=0.06, y0=0.07, label_txt=None, n=2, color="k", halo=True, bg=False):
    """Alternating black/white bar, axes-fraction placement, length given in data units via unit_m_per_data."""
    xl, yl = ax.get_xlim(), ax.get_ylim()
    L = length_m / unit_m_per_data
    X0 = xl[0] + x0 * (xl[1] - xl[0]); Y0 = yl[0] + y0 * (yl[1] - yl[0]); H = 0.012 * (yl[1] - yl[0])
    if bg:
        ax.add_patch(Rectangle((X0 - 0.015 * (xl[1] - xl[0]), Y0 - 0.025 * (yl[1] - yl[0])), L + 0.03 * (xl[1] - xl[0]) + 0.06 * (xl[1] - xl[0]),
                               H + 0.085 * (yl[1] - yl[0]), fc="white", ec="none", alpha=0.75, zorder=9))
    for i in range(n):
        ax.add_patch(Rectangle((X0 + i * L / n, Y0), L / n, H, fc="k" if i % 2 == 0 else "white", ec="k", lw=0.5, zorder=10))
    txt = label_txt or f"{length_m:g} m"
    ax.text(X0, Y0 + H + 0.012 * (yl[1] - yl[0]), "0", ha="center", va="bottom", fontsize=5.5, zorder=10)
    ax.text(X0 + L, Y0 + H + 0.012 * (yl[1] - yl[0]), txt, ha="center", va="bottom", fontsize=5.5, zorder=10)

def ptr(ax, x, y, kind, size=22, z=8):
    ax.scatter([x], [y], s=size, c=COL[kind], edgecolors="white" if kind == "none" else "#222222", linewidths=0.6, zorder=z)

def panel_letter(ax, s, dark=False):
    t = ax.text(0.015, 0.985, s, transform=ax.transAxes, ha="left", va="top", fontsize=9, fontweight="bold",
                color="white" if dark else "black", zorder=20)
    t.set_path_effects(DHALO if dark else HALO)


# ---------------------------------------------------------------- graphic rendition of the imagery
WATER, MARSH, SAND, BUILT, EDGE = "#dbe9f4", "#cfe3bd", "#f3ead2", "#b4b4b4", "#6f6f6f"
def classify(path, sig=2.0, gb=14, bright=158, tex_t=14, tw=9):
    raw = np.array(Image.open(path).convert("RGB")).astype(float)
    img = ndi.gaussian_filter(raw, (sig, sig, 0)); r, g, b = img[..., 0], img[..., 1], img[..., 2]
    lum = 0.3 * r + 0.59 * g + 0.11 * b; lr = 0.3 * raw[..., 0] + 0.59 * raw[..., 1] + 0.11 * raw[..., 2]
    m = ndi.uniform_filter(lr, tw); tex = ndi.gaussian_filter(np.sqrt(np.maximum(ndi.uniform_filter(lr ** 2, tw) - m ** 2, 0)), 3)
    veg = ((g - b) > gb) & (g > r - 10)
    cls = np.zeros(lum.shape, int); cls[veg] = 1                       # 0 water, 1 marsh
    cls[(~veg) & (lum >= bright)] = 2                                 # 2 sand
    cls[(tex > tex_t) & (lum > 95) & (~veg)] = 3                      # 3 built
    return cls

def clip_seaward(cls):
    """Everything seaward of the beach's outer edge is ocean (removes surf-zone foam classified as sand/built)."""
    lab, n = ndi.label(cls == 2)
    big = lab == (1 + np.argmax(ndi.sum(np.ones_like(lab), lab, range(1, n + 1))))
    rows = np.arange(cls.shape[0]); edge = np.array([np.where(big[r])[0].max() if big[r].any() else np.nan for r in rows])
    ok = ~np.isnan(edge); edge = np.interp(rows, rows[ok], edge[ok])
    out = cls.copy()
    for r in rows: out[r, int(edge[r]) + 2:] = 0
    return out

def clean(cls, min_area, sigma, merge_built=False):
    cls = cls.copy()
    if merge_built: cls[cls == 3] = 2
    for k in np.unique(cls):                                           # drop specks, fill from nearest class
        lab, n = ndi.label(cls == k)
        if n:
            sizes = ndi.sum(np.ones_like(lab), lab, range(1, n + 1))
            lim = min_area.get(int(k), 0) if isinstance(min_area, dict) else min_area
            small = np.isin(lab, 1 + np.where(sizes < lim)[0]); cls[small] = -1
    idx = ndi.distance_transform_edt(cls < 0, return_distances=False, return_indices=True)
    cls = cls[tuple(idx)]
    P = np.stack([ndi.gaussian_filter((cls == k).astype(float), sigma) for k in range(4)])
    return P

def draw_graphic(ax, P, ext, colors, draw=(1, 2, 3), edge_lw=0.25, z=0):
    h, w = P.shape[1:]
    x = np.linspace(ext[0], ext[2], w); y = np.linspace(ext[1], ext[3], h)
    ax.set_facecolor(colors[0])
    for k in draw:
        others = np.max(np.delete(P, k, axis=0), axis=0); M = (P[k] - others)[::-1]
        ax.contourf(x, y, M, levels=[0, 2], colors=[colors[k]], antialiased=True, zorder=z + k)
        if edge_lw: ax.contour(x, y, M, levels=[0], colors=[EDGE], linewidths=edge_lw, zorder=z + 4)

# ---------------------------------------------------------------- figure grid (mm)
W, H = 180, 172
fig = plt.figure(figsize=(W * MM, H * MM))
def axmm(x, y, w, h, **kw): return fig.add_axes([x / W, y / H, w / W, h / H], **kw)
rowbot_y, hC = 20, 70
hA = 70
rowtop_y = rowbot_y + hC + 8
WL, WR, XR = 92, 76, 100
axA = axmm(4, rowtop_y, WL, hA, projection=ccrs.AlbersEqualArea(central_longitude=-96, central_latitude=37.5,
                                                                 standard_parallels=(29.5, 45.5)))
axB = axmm(XR, rowtop_y, WR, hA, projection=ccrs.UTM(18))
axC = axmm(4, rowbot_y, WL, hC)
axD = axmm(XR, rowbot_y, WR, hC)

# ---------------------------------------------------------------- a: U.S.
PC = ccrs.PlateCarree()
_pts = [(lo, la) for lo in np.linspace(-124.5, -67.5, 40) for la in (24.2, 49.8)] + [(lo, la) for lo in (-124.5, -67.5) for la in np.linspace(24.2, 49.8, 40)]
_xy = np.array([axA.projection.transform_point(lo, la, PC) for lo, la in _pts])
_cx, _cy = (_xy[:, 0].min() + _xy[:, 0].max()) / 2, (_xy[:, 1].min() + _xy[:, 1].max()) / 2
_w, _h = np.ptp(_xy[:, 0]), np.ptp(_xy[:, 1]); _asp = WL / hA
if _w / _h < _asp: _w = _h * _asp
else: _h = _w / _asp
axA.set_extent([_cx - _w / 2, _cx + _w / 2, _cy - _h / 2, _cy + _h / 2], crs=axA.projection)
axA.set_facecolor(OCEAN)
axA.add_geometries(land, PC, facecolor=LAND, edgecolor="#6f6f6f", linewidth=0.3)
for pc, g in states:
    axA.add_geometries([g], PC, facecolor=NJC if pc != "NJ" else "none", edgecolor=STATEC, linewidth=0.2)
    if pc == "NJ": axA.add_geometries([g], PC, facecolor="#e4572e", edgecolor="#8c2a10", linewidth=0.4, zorder=5)
axA.plot([SITE[0]], [SITE[1]], marker="o", ms=3.2, mfc="white", mec="k", mew=0.7, transform=PC, zorder=8)
axA.annotate("Whale Beach", xy=(SITE[0], SITE[1]), xycoords=PC._as_mpl_transform(axA), xytext=(-1.0, -0.2), textcoords="offset fontsize",
             fontsize=6, ha="right", va="top", path_effects=HALO, zorder=9)
axA.text(0.50, 0.50, "United States", transform=axA.transAxes, fontsize=8, color="#6a6a6a", ha="center", va="center", style="italic")
axA.text(0.89, 0.30, "Atlantic\nOcean", transform=axA.transAxes, fontsize=6, color="#4a6a85", ha="center", style="italic")
axA.text(0.10, 0.38, "Pacific\nOcean", transform=axA.transAxes, fontsize=6, color="#4a6a85", ha="center", style="italic")
# NJ bounding box shows panel b
x0, y0 = axA.projection.transform_point(-76.0, 38.7, PC); x1, y1 = axA.projection.transform_point(-73.6, 41.5, PC)
axA.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc="none", ec="k", lw=0.7, zorder=7))
north_arrow(axA, x=0.06, y=0.80)
scale_bar(axA, 500e3, 1.0, x0=0.40, y0=0.07, label_txt="500 km")
panel_letter(axA, "a")

# ---------------------------------------------------------------- b: New Jersey
UTM = ccrs.UTM(18)
cx, cy = UTM.transform_point(-74.62, 40.15, PC)
HB = 300e3; WB = HB * WR / hA
axB.set_xlim(cx - WB / 2, cx + WB / 2); axB.set_ylim(cy - HB / 2, cy + HB / 2)
axB.set_facecolor(OCEAN)
axB.add_geometries(land, PC, facecolor=LAND, edgecolor="#6f6f6f", linewidth=0.3)
for pc, g in states:
    axB.add_geometries([g], PC, facecolor=NJC if pc != "NJ" else "#fbd5c8", edgecolor="#8a8a8a", linewidth=0.3)
axB.set_xlim(cx - WB / 2, cx + WB / 2); axB.set_ylim(cy - HB / 2, cy + HB / 2)
for nm, lo, la, dx, dy, ha in [("New York City", -74.006, 40.713, -0.6, 0.3, "right"), ("Philadelphia", -75.165, 39.953, -0.6, 0.0, "right"),
                               ("Atlantic City", -74.423, 39.364, -0.6, -0.3, "right")]:
    x, y = UTM.transform_point(lo, la, PC)
    axB.plot(x, y, "o", ms=1.8, mfc="k", mec="white", mew=0.3, zorder=6)
    label(axB, x + dx * 4e3, y + dy * 4e3, nm, fontsize=5.2, ha=ha, va="center")
for nm, x, y, r, c in [("NJ", 0.46, 0.60, 0, "#8c2a10"), ("PA", 0.12, 0.60, 0, "#6a6a6a"), ("NY", 0.38, 0.93, 0, "#6a6a6a"),
                       ("DE", 0.12, 0.10, 0, "#6a6a6a")]:
    axB.text(x, y, nm, transform=axB.transAxes, fontsize=7, color=c, ha="center", fontweight="bold" if nm == "NJ" else "normal")
axB.text(0.80, 0.36, "Atlantic\nOcean", transform=axB.transAxes, fontsize=6, color="#4a6a85", ha="center", style="italic")
sx, sy = UTM.transform_point(SITE[0], SITE[1], PC)
# panel c extent box
cxm = Transformer.from_crs(4326, 32618, always_xy=True)
axB.plot(sx, sy, "s", ms=4.5, mfc="none", mec="#c4161c", mew=0.9, zorder=8)
axB.annotate("Whale Beach", xy=(sx, sy), xytext=(sx + 38e3, sy + 14e3), fontsize=6, fontweight="bold", ha="left", va="center",
             arrowprops=dict(arrowstyle="-", lw=0.5, color="k", shrinkA=0, shrinkB=2), path_effects=HALO, zorder=9)
north_arrow(axB, x=0.90, y=0.83)
scale_bar(axB, 50e3, 1.0, x0=0.58, y0=0.04, label_txt="50 km")
panel_letter(axB, "b")

# ---------------------------------------------------------------- c: Strathmere / Whale Beach (imagery)
Tm = Transformer.from_crs(4326, 3857, always_xy=True)
ext = json.load(open(os.path.join(DATA, "extents_3857.json")))
ext["panel_d"] = ext["panel_d_nj"]
kc = 1 / math.cos(math.radians(39.19))
bc = ext["panel_c"]
Pc = clean(classify(os.path.join(DATA, "panel_c.jpg")), min_area=60, sigma=2.2)
draw_graphic(axC, Pc, bc, [WATER, MARSH, SAND, BUILT])
axC.set_xlim(bc[0], bc[2]); axC.set_ylim(bc[1], bc[3]); axC.set_aspect("equal")
axC.set_xticks([]); axC.set_yticks([])
gx, gy = Tm.transform(GAUGE["lon"], GAUGE["lat"])
axC.scatter([gx], [gy], marker="^", s=40, c="#ff2d95", edgecolors="#222222", linewidths=0.6, zorder=9)
label(axC, gx - 80 * kc, gy + 5 * kc, "Strathmere\ntide gauge", fontsize=5.8, ha="right", va="center", fontweight="bold", halo=True)
sx3, sy3 = Tm.transform(SITE[0], SITE[1])
# panel d extent box
bd = ext["panel_d"]
axC.add_patch(Rectangle((bd[0], bd[1]), bd[2] - bd[0], bd[3] - bd[1], fc="none", ec="k", lw=0.8, zorder=8))
for c in cores:
    x, y = Tm.transform(c["lon"], c["lat"]); ptr(axC, x, y, c["kind"], size=5, z=8)
ux, uy = math.sin(math.radians(SEAWARD)), math.cos(math.radians(SEAWARD))
bw_c, bh_c = bc[2] - bc[0], bc[3] - bc[1]
# seaward arrow drawn in open water, parallel to the shore-normal at Whale Beach
a0 = (bc[0] + 0.66 * bw_c, bc[1] + 0.40 * bh_c); L = 0.17 * bw_c
a1 = (a0[0] + L * ux, a0[1] + L * uy)
axC.add_patch(FancyArrowPatch(a0, a1, arrowstyle="simple,head_width=5,head_length=6,tail_width=1.6", fc="#222222", ec="white", lw=0.5, zorder=9))
label(axC, a0[0] - 0.012 * bw_c, a0[1] + 0.0 * bh_c, "Seaward", fontsize=6.5, ha="right", va="center", fontweight="bold")
label(axC, bc[0] + 0.015 * bw_c, bc[1] + 0.74 * bh_c, "Back-barrier\nsalt marsh", fontsize=6, ha="left", va="center", style="italic", color="#1b2a10")
t_ = label(axC, bc[0] + 0.84 * bw_c, bc[1] + 0.20 * bh_c, "Atlantic\nOcean", fontsize=6.5, ha="center", va="center", style="italic", color="#4a6a85")
label(axC, bc[0] + 0.82 * bw_c, bc[1] + 0.78 * bh_c, "Strathmere", fontsize=6.2, ha="center", va="center", fontweight="bold")
tt = axC.annotate("Whale Beach\n(cores)", xy=(sx3 - 25 * kc, sy3 + 5 * kc), xytext=(sx3 - 600 * kc, sy3 + 330 * kc), fontsize=6.2, ha="center", va="center",
                  fontweight="bold", arrowprops=dict(arrowstyle="-", lw=0.6, color="k", shrinkA=0, shrinkB=3), zorder=9)
tt.set_path_effects(HALO)
north_arrow(axC, x=0.93, y=0.80, color="k")
scale_bar(axC, 500, kc, x0=0.60, y0=0.06, label_txt="500 m")
panel_letter(axC, "c")

# ---------------------------------------------------------------- d: cores (graphic) with pins
from matplotlib.path import Path as MPath
bd = ext["panel_d"]
# Features traced from the NJOGIS 2020 orthophoto (frame = 2000 x 1842 px over the panel extent): water, beach, roads, buildings
GEO = dict(
    canal=[(0, 0), (340, 0), (325, 60), (300, 150), (262, 240), (225, 320), (190, 410), (150, 500), (110, 570), (60, 650), (25, 720), (0, 740)],
    pond=[(452, 305), (520, 345), (600, 392), (690, 430), (740, 452), (722, 500), (660, 515), (560, 468), (498, 420), (455, 350)],
    ocean=[(2000, 640), (1960, 780), (1890, 900), (1830, 1020), (1760, 1160), (1690, 1300), (1620, 1450), (1540, 1610), (1470, 1730), (1400, 1842), (2000, 1842)],
    main_road=[(1845, 0), (565, 1842)],
    local_road=[(360, 345), (500, 440), (640, 540), (760, 620), (880, 690), (990, 780), (1100, 880), (1180, 925)],
    sw_road=[(0, 1485), (130, 1560), (250, 1650), (380, 1760), (520, 1842)],
    buildings=[(760, 495, 70, 60), (845, 545, 90, 70), (895, 575, 60, 50), (945, 640, 70, 70), (995, 690, 70, 60), (1075, 745, 50, 50), (1150, 820, 90, 70),
               (1290, 625, 85, 65), (1235, 725, 70, 65), (1475, 320, 120, 90), (1660, 90, 100, 80),
               (1060, 950, 80, 60), (1040, 1020, 60, 50), (960, 1090, 80, 70), (925, 1165, 70, 60), (880, 1235, 70, 55), (845, 1280, 60, 55), (720, 1470, 90, 70),
               (140, 1450, 130, 110), (215, 1520, 70, 70), (20, 1625, 60, 60), (90, 1670, 80, 70),
               (330, 200, 80, 70), (230, 350, 80, 60), (180, 440, 60, 60), (120, 510, 60, 55)])
_bw, _bh = bd[2] - bd[0], bd[3] - bd[1]
def D(pts): return [(bd[0] + x / 2000 * _bw, bd[3] - y / 1842 * _bh) for x, y in pts]
axD.set_facecolor(MARSH)
axD.add_patch(Polygon(D([(1875, 0), (2000, 0), (2000, 1842), (595, 1842)]), fc=SAND, ec="none", zorder=1))             # beach (east of the road)
for k in ("canal", "pond", "ocean"):
    axD.add_patch(Polygon(D(GEO[k]), fc=WATER, ec=EDGE, lw=0.25, zorder=2))
_ppm = 0.158 * 72 / 25.4                                                                                              # points per metre at this scale
for k, wm in (("main_road", 12), ("local_road", 6), ("sw_road", 6)):
    xy = np.array(D(GEO[k])); axD.plot(xy[:, 0], xy[:, 1], color=BUILT, lw=wm * _ppm, solid_capstyle="butt", solid_joinstyle="round", zorder=3)
_d = np.array([-1280, 1842.]); _d /= np.linalg.norm(_d); _n = np.array([_d[1], -_d[0]]) * -1
for cx, cy, w, h in GEO["buildings"]:
    c = np.array([cx, cy]); pts = [c + sx * w / 2 * _n + sy * h / 2 * _d for sx, sy in ((-1, -1), (1, -1), (1, 1), (-1, 1))]
    axD.add_patch(Polygon(D(pts), fc=BUILT, ec=EDGE, lw=0.3, zorder=4))
axD.set_xlim(bd[0], bd[2]); axD.set_ylim(bd[1], bd[3]); axD.set_aspect("equal")
axD.set_xticks([]); axD.set_yticks([])
# teardrop pin whose tip sits exactly on the coordinate
_th = np.radians(np.linspace(-27, 207, 60)); _head = np.c_[np.cos(_th), 2.2 + np.sin(_th)]
PIN = MPath(np.vstack([_head, [[0, 0]], _head[:1]]), [MPath.MOVETO] + [MPath.LINETO] * (len(_head)) + [MPath.CLOSEPOLY])
def pin(ax, x, y, kind, size=420, z=8):
    ax.scatter([x], [y], marker=PIN, s=size, c=COL[kind], edgecolors="white" if kind == "none" else "#222222", linewidths=0.6, zorder=z)
cxy = {c["name"]: Tm.transform(c["lon"], c["lat"]) for c in cores}
mx = np.mean([v[0] for v in cxy.values()]); my_ = np.mean([v[1] for v in cxy.values()])
WHALO = HALO
def place(names, xcol, ha, spacing):
    names = sorted(names, key=lambda n: -cxy[n][1]); ys = [cxy[n][1] for n in names]
    for i in range(1, len(ys)): ys[i] = min(ys[i], ys[i - 1] - spacing)
    shift = np.mean([cxy[n][1] for n in names]) - np.mean(ys); ys = [y + shift for y in ys]
    for n, yl in zip(names, ys):
        x, y = cxy[n]
        axD.plot([x, xcol], [y, yl], color="k", lw=0.4, zorder=6, solid_capstyle="round")
        t = axD.text(xcol + (3 * kc if ha == "left" else -3 * kc), yl, n, fontsize=5.8, ha=ha, va="center", fontweight="bold", color="k", zorder=9)
        t.set_path_effects(WHALO)
west = [n for n in cxy if cxy[n][0] < mx - 10 * kc or n == "WB2505"]
east = [n for n in cxy if n not in west]
place(west, mx - 62 * kc, "right", 15 * kc)
place(east, mx + 66 * kc, "left", 15 * kc)
for c in sorted(cores, key=lambda c: c["lat"]):
    x, y = cxy[c["name"]]; pin(axD, x, y, c["kind"])
bw, bh = bd[2] - bd[0], bd[3] - bd[1]
a0 = (bd[0] + 0.66 * bw, bd[1] + 0.40 * bh); L = 0.24 * bw
a1 = (a0[0] + L * ux, a0[1] + L * uy)
axD.add_patch(FancyArrowPatch(a0, a1, arrowstyle="simple,head_width=5,head_length=6,tail_width=1.8", fc="#222222", ec="white", lw=0.5, zorder=9))
t_ = axD.text(a0[0] - 0.02 * bw, a0[1] + 0.035 * bh, "Seaward", fontsize=6.5, ha="left", va="bottom", fontweight="bold", color="k", zorder=9); t_.set_path_effects(WHALO)
t_ = axD.text(bd[0] + 0.90 * bw, bd[1] + 0.07 * bh, "Atlantic Ocean", fontsize=6.3, ha="center", va="center", style="italic", color="#4a6a85", zorder=9); t_.set_path_effects(WHALO)
north_arrow(axD, x=0.94, y=0.80, color="k")
scale_bar(axD, 100, kc, x0=0.05, y0=0.05, label_txt="100 m")
panel_letter(axD, "d")

# connector lines c -> d
for (xa, ya, xb, yb) in [(bd[2], bd[3], 0, 1), (bd[2], bd[1], 0, 0)]:
    pass
fig.add_artist(Line2D([0, 0], [0, 0], alpha=0))

# ---------------------------------------------------------------- legend (bottom)
items = [Line2D([], [], marker="o", ls="", mfc=COL["eDNA"], mec="#222222", mew=0.6, ms=6, label="eDNA analysis only"),
         Line2D([], [], marker="o", ls="", mfc=COL["sed"], mec="#222222", mew=0.6, ms=6, label="Sediment analysis only"),
         Line2D([], [], marker="o", ls="", mfc=COL["both"], mec="#222222", mew=0.6, ms=6, label="eDNA + sediment analysis"),
         Line2D([], [], marker="o", ls="", mfc=COL["none"], mec="white", mew=0.6, ms=6, label="Cored; no analysis in this study"),
         Line2D([], [], marker="^", ls="", mfc="#ff2d95", mec="#222222", mew=0.6, ms=6, label="Strathmere tide gauge"),
         Line2D([], [], marker="s", ls="", mfc="none", mec="#c4161c", mew=0.9, ms=5, label="Study area (b)"),
         Patch(fc=MARSH, ec=EDGE, lw=0.4, label="Marsh / dune vegetation (c, d)"), Patch(fc=SAND, ec=EDGE, lw=0.4, label="Beach / sand (c, d)"),
         Patch(fc=BUILT, ec=EDGE, lw=0.4, label="Roads / buildings (c, d)")]
leg = fig.legend(handles=items, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3, frameon=False, fontsize=6.2,
                 columnspacing=1.6, handletextpad=0.3, labelspacing=0.5)

out = os.path.join(HERE, "Figure_FieldSite")
fig.savefig(out + ".pdf", dpi=600, bbox_inches="tight", pad_inches=0.03)
fig.savefig(out + ".png", dpi=400, bbox_inches="tight", pad_inches=0.03)
print("saved", out)

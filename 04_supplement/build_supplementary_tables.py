"""Rebuild 'Master Supplementary Document.docx' with one consistent table style / numbering and recomputed numbers.

Style follows the existing Table S1 (loss on ignition): Arial, 10 pt bold-label captions above, 9 pt three-line tables
(no vertical rules), 8 pt notes below. Table S1 and the phylum-composition figure are carried over unchanged.
"""
import csv, os
from docx import Document
from docx.shared import Pt
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.text import WD_ALIGN_PARAGRAPH

HOME = os.path.expanduser("~/Desktop")
SRC = os.path.join(HOME, "Master Supplementary Document.docx")
OUT = os.path.join(HOME, "Master Supplementary Document (updated).docx")
DIR = os.path.join(HOME, "WB Data Analysis ", "Supplementary_Tables_Update", "sqrt_mednorm")   # square-root / median-normalised analyses
rd = lambda f: list(csv.DictReader(open(os.path.join(DIR, f))))

# ------------------------------------------------------------------ number formats
def fp(p, d=3):
    p = float(p); return "< 0.001" if p < 0.001 else f"{p:.{d}f}"
f2 = lambda x: f"{float(x):.2f}"
CORES = ["WB2501", "WB2506", "WB2508", "WB2509"]
NE = "Nor’easter"
CMP = {"Hurricane vs Peat": "Hurricane vs. peat", "Noreaster vs Peat": f"{NE} vs. peat", "Hurricane vs Noreaster": f"Hurricane vs. {NE}"}
TYPE = {"Noreaster": f"{NE} overwash", "Hurricane": "Hurricane overwash", "Peat": "Salt marsh peat"}

# ------------------------------------------------------------------ low-level helpers
def set_run(r, size, bold=False, italic=False):
    r.font.name = "Arial"; r.font.size = Pt(size); r.bold = bold; r.italic = italic
    rf = r._r.get_or_add_rPr().rFonts
    for a in ("w:eastAsia", "w:cs", "w:hAnsi", "w:ascii"): rf.set(qn(a), "Arial")

def para_fmt(p, before=None, after=None, keep_next=False, keep_lines=False, line=None):
    pf = p.paragraph_format
    if before is not None: pf.space_before = Pt(before)
    if after is not None: pf.space_after = Pt(after)
    if line is not None: pf.line_spacing = line
    pf.keep_with_next = True if keep_next else None; pf.keep_together = True if keep_lines else None

def border_xml(tc_pr, top=None, bottom=None):
    b = OxmlElement("w:tcBorders")
    for side, spec in (("top", top), ("left", None), ("bottom", bottom), ("right", None)):
        e = OxmlElement(f"w:{side}")
        if spec: e.set(qn("w:val"), "single"); e.set(qn("w:sz"), str(spec[0])); e.set(qn("w:color"), spec[1])
        else: e.set(qn("w:val"), "none"); e.set(qn("w:sz"), "0"); e.set(qn("w:color"), "FFFFFF")
        e.set(qn("w:space"), "0"); b.append(e)
    tc_pr.append(b)

def style_cell(cell, width, top=None, bottom=None):
    tc_pr = cell._tc.get_or_add_tcPr()
    for tag in ("w:tcW", "w:tcBorders", "w:tcMar", "w:vAlign"):
        for e in tc_pr.findall(qn(tag)): tc_pr.remove(e)
    w = OxmlElement("w:tcW"); w.set(qn("w:w"), str(width)); w.set(qn("w:type"), "dxa"); tc_pr.insert(0, w)
    border_xml(tc_pr, top, bottom)
    m = OxmlElement("w:tcMar")
    for side, v in (("top", 30), ("left", 80), ("bottom", 30), ("right", 80)):
        e = OxmlElement(f"w:{side}"); e.set(qn("w:w"), str(v)); e.set(qn("w:type"), "dxa"); m.append(e)
    tc_pr.append(m)
    v = OxmlElement("w:vAlign"); v.set(qn("w:val"), "center"); tc_pr.append(v)

class Builder:
    def __init__(self, doc): self.doc = doc

    def caption(self, label, text, first=False, keep_next=True):
        p = self.doc.add_paragraph(); set_run(p.add_run(label + " "), 10, bold=True); set_run(p.add_run(text), 10)
        para_fmt(p, before=0 if first else 18, after=6, keep_next=keep_next, keep_lines=True, line=1.0); p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        return p

    def note(self, text):
        p = self.doc.add_paragraph(); set_run(p.add_run(text), 8)
        para_fmt(p, before=6, after=2, keep_lines=True, line=1.0); p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY; return p

    def table(self, header, body, widths, align, block_ends=(), merges=(), bold_cells=(), italic_rows=(), group_header=None):
        """header: list of header rows (each a list); body: list of rows; widths in dxa (sum 9360);
        align: per column 'l'/'c'; block_ends: body-row indices followed by a thin rule;
        group_header: optional list of (start_col, end_col, text) shown in an extra top header row."""
        hdr_rows = ([None] if group_header else []) + list(header)
        nr, nc = len(hdr_rows) + len(body), len(widths)
        t = self.doc.add_table(rows=nr, cols=nc); t.autofit = False
        tblPr = t._tbl.tblPr
        for e in list(tblPr):
            if e.tag in (qn("w:tblStyle"), qn("w:tblBorders"), qn("w:tblW"), qn("w:tblLayout")): tblPr.remove(e)
        w = OxmlElement("w:tblW"); w.set(qn("w:w"), str(sum(widths))); w.set(qn("w:type"), "dxa")
        lay = OxmlElement("w:tblLayout"); lay.set(qn("w:type"), "fixed")
        look = tblPr.find(qn("w:tblLook"))
        if look is not None: look.addprevious(w); look.addprevious(lay)
        else: tblPr.append(w); tblPr.append(lay)
        for i, c in enumerate(t.columns): c.width = __import__("docx").shared.Twips(widths[i])
        n_h = len(hdr_rows)
        def put(cell, text, size=9, bold=False, italic=False, al="l"):
            cell.text = ""; p = cell.paragraphs[0]; set_run(p.add_run(text), size, bold, italic)
            para_fmt(p, before=0, after=0, line=1.0)
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER if al == "c" else WD_ALIGN_PARAGRAPH.LEFT
            pPr = p._p.get_or_add_pPr(); old = pPr.find(qn("w:rPr"))
            if old is not None: pPr.remove(old)
            mark = OxmlElement("w:rPr"); sz = OxmlElement("w:sz"); sz.set(qn("w:val"), str(int(size * 2))); szc = OxmlElement("w:szCs"); szc.set(qn("w:val"), str(int(size * 2)))
            mark.append(sz); mark.append(szc); pPr.append(mark)
        # header
        ri = 0
        if group_header:
            for (a, b, text) in group_header:
                c = t.cell(0, a).merge(t.cell(0, b)) if b > a else t.cell(0, a)
                put(c, text, bold=True, al="c")
            ri = 1
        for row in header:
            for j, text in enumerate(row): put(t.cell(ri, j), text, bold=True, al=align[j])
            ri += 1
        for bi, row in enumerate(body):
            for j, text in enumerate(row):
                put(t.cell(n_h + bi, j), text, bold=((bi, j) in bold_cells) or (bi in italic_rows and False), italic=(bi in italic_rows), al=align[j])
        for (bi, a, b) in merges:                                   # merge cells of a body row (sub-headings)
            r = n_h + bi; label = body[bi][a]; c = t.cell(r, a).merge(t.cell(r, b)); put(c, label, bold=True, al="l")
        # borders / widths / row properties
        for r in range(nr):
            tr = t.rows[r]._tr; trPr = tr.get_or_add_trPr()
            cs = OxmlElement("w:cantSplit"); trPr.append(cs)
            if r < n_h: trPr.append(OxmlElement("w:tblHeader"))
            last = (r == nr - 1)
            top = (12, "000000") if r == 0 else None
            bottom = (6, "000000") if r == n_h - 1 else ((12, "000000") if last else None)
            if group_header and r == 0: bottom = (4, "7F7F7F")
            if r >= n_h and (r - n_h) in block_ends and not last: bottom = (4, "BFBFBF")
            seen = set(); j = 0
            for tc in tr.findall(qn("w:tc")):
                span = tc.tcPr.find(qn("w:gridSpan")) if tc.tcPr is not None else None
                n = int(span.get(qn("w:val"))) if span is not None else 1
                from docx.table import _Cell
                style_cell(_Cell(tc, t), sum(widths[j:j + n]), top, bottom); j += n
        # keep table rows with the next one so short tables do not split across pages
        for r in range(nr - 1):
            for tc in t.rows[r]._tr.findall(qn("w:tc")):
                for p in tc.findall(qn("w:p")):
                    pPr = p.find(qn("w:pPr"));
                    if pPr is None: pPr = OxmlElement("w:pPr"); p.insert(0, pPr)
                    if pPr.find(qn("w:keepNext")) is None and len(body) <= 20: pPr.insert(0, OxmlElement("w:keepNext"))
        return t

# ------------------------------------------------------------------ load source and keep reusable pieces
doc = Document(SRC); body = doc.element.body
kids = list(body.iterchildren()); sect = kids[-1]
def text(el): return "".join(t.text or "" for t in el.iter(qn("w:t")))
title = kids[0]
i_cap = next(i for i, k in enumerate(kids) if k.tag == qn("w:p") and text(k).startswith("Table S1."))
s1 = [kids[i_cap], kids[i_cap + 1], kids[i_cap + 2], kids[i_cap + 3]]         # caption, table, two notes
fig_p = next(k for k in kids if k.tag == qn("w:p") and k.xpath(".//w:drawing"))
assert text(s1[2]).startswith("Sand,") and text(s1[3]).startswith("aGrain"), [text(x)[:30] for x in s1]
for k in kids[:-1]: body.remove(k)

# compact the paragraphs in Table S1 so it matches the other tables
from docx.table import Table
from docx.text.paragraph import Paragraph
t1 = Table(s1[1], doc)
for row in t1.rows:
    for cell in row.cells:
        for p in cell.paragraphs: para_fmt(p, before=0, after=0, line=1.0)
for r in t1.rows[:1]:
    pass

B = Builder(doc)
# order: title, S1 (unchanged), S2..S6 (deposit types), figure, S7..S8 (among cores)
sect.addprevious(title)
for e in s1[:1]: sect.addprevious(e)
for r_ in Paragraph(title, doc).runs: set_run(r_, 14, bold=True)
para_fmt(Paragraph(title, doc), after=10)
cap = Paragraph(s1[0], doc); para_fmt(cap, before=0, after=6, keep_next=True, line=1.0)
for e in s1[1:]: sect.addprevious(e)
for p in (Paragraph(s1[2], doc), Paragraph(s1[3], doc)): para_fmt(p, before=3, after=2, line=1.0)

# detach the document-end sectPr so builder appends land before it automatically (python-docx keeps sectPr last)
# ---- S2 deposit definitions
defs = rd("S_deposit_definitions.csv"); rng = {"WB2501": "20–26", "WB2506": "24–27", "WB2508": "20–30", "WB2509": "20–24"}
rows, ends = [], []
for k in CORES:
    for t in ("Noreaster", "Hurricane", "Peat"):
        r = next(x for x in defs if x["core"] == k and x["type"] == t)
        rule = rng[k] if t == "Noreaster" else ("≥ 40" if t == "Hurricane" else "Remaining")
        rows.append([k if t == "Noreaster" else "", TYPE[t], rule, r["n"], r["intervals"].replace("-", "–")])
    ends.append(len(rows) - 1)
B.caption("Table S2.", f"Deposit-type definitions and sample counts for the four eDNA cores. Samples were assigned by the depth of the top of each sampled interval: {NE.lower()} overwash, "
          "top depth within the stated range (inclusive); hurricane overwash, intervals starting at 40 cm or deeper; salt marsh peat, all remaining intervals. "
          "These definitions are used in the PCoA and dbRDA figures and in Tables S3–S7.")
B.table([["Core", "Deposit type", "Top depth (cm)", "n", "Sampled intervals (cm)"]], rows, [900, 2000, 1450, 500, 4510], ["l", "l", "c", "c", "l"], block_ends=ends)

# ---- S3 average Bray-Curtis between deposit types
bc = rd("S_bc_between_groups.csv"); rows, ends = [], []
for k in CORES:
    for c in ("Hurricane vs Peat", "Noreaster vs Peat", "Hurricane vs Noreaster"):
        r = next(x for x in bc if x["core"] == k and x["comparison"] == c); d = round(float(r["dissimilarity"]), 1)
        rows.append([k if c.startswith("Hurricane vs Peat") else "", CMP[c], f"{d:.1f}%", f"{100 - d:.1f}%"])
    ends.append(len(rows) - 1)
B.caption("Table S3.", "Average Bray–Curtis dissimilarity and similarity between deposit types within each core. Values are the mean of all pairwise Bray–Curtis dissimilarities "
          "(of square-root-transformed, median-normalized counts) between samples of the two deposit types compared; similarity = 100 − dissimilarity.")
B.table([["Core", "Comparison", "Dissimilarity", "Similarity"]], rows, [1500, 3860, 2000, 2000], ["l", "l", "c", "c"], block_ends=ends)

# ---- S4 PERMANOVA (3 deposit types)
pm = rd("S_permanova_3groups.csv"); rows, ends = [], []
for k in CORES:
    r = next(x for x in pm if x["core"] == k); g = lambda c, d=4: f"{float(r[c]):.{d}f}"
    rows.append([k, "Deposit type", r["df_model"], g("SS_model"), g("MS_model"), g("F"), f"{float(r['p']):.4f}", r["perms"]])
    rows.append(["", "Residual", r["df_resid"], g("SS_resid"), g("MS_resid"), "", "", ""])
    rows.append(["", "Total", r["df_total"], g("SS_total"), "", "", "", ""]); ends.append(len(rows) - 1)
B.caption("Table S4.", "PERMANOVA of community composition among the three deposit types in each core (Bray–Curtis dissimilarity of square-root-transformed, median-normalized counts; one factor, deposit type; "
          "9,999 permutations). df, degrees of freedom; SS, sum of squares; MS, mean square.")
B.table([["Core", "Source", "df", "SS", "MS", "Pseudo-F", "P(perm)", "Perms"]], rows, [1000, 1700, 700, 1100, 1100, 1200, 1300, 1260],
        ["l", "l", "c", "c", "c", "c", "c", "c"], block_ends=ends)

# ---- S5 pairwise PERMANOVA between deposit types
pw = rd("S_permanova_pairwise.csv"); cb = rd("S_hurricane_vs_noreaster_combined.csv")[0]; rows, ends, bold = [], [], set()
for k in CORES:
    for c in ("Hurricane vs Peat", "Noreaster vs Peat", "Hurricane vs Noreaster"):
        r = next(x for x in pw if x["core"] == k and x["comparison"] == c)
        rows.append([k if c.startswith("Hurricane vs Peat") else "", CMP[c], f"{r['n1']}, {r['n2']}", f2(r["F"]), f2(r["R2"]), fp(r["p"]), fp(r["min_p"]), fp(r["p_holm"]),
                     fp(r["PERMDISP_p"]) if r["PERMDISP_p"] not in ("", "NA") else "–"])
        if float(r["p_holm"]) < 0.05: bold.add((len(rows) - 1, 7))
    ends.append(len(rows) - 1)
B.caption("Table S5.", "Pairwise PERMANOVA between deposit types within each core (Bray–Curtis dissimilarity of square-root-transformed, median-normalized counts). n, number of samples in the first and second group; "
          "p, exact permutation p-value (all possible relabelings of the samples); min. p, smallest p-value attainable with those group sizes; p (Holm), p adjusted across the three comparisons within a core "
          "(bold, < 0.05); PERMDISP p, permutation test of homogeneity of multivariate dispersion (9,999 permutations; uninformative when a group has two samples).")
B.table([["Core", "Comparison", "n", "Pseudo-F", "R²", "p", "Min. p", "p (Holm)", "PERMDISP p"]], rows, [900, 2250, 650, 1050, 650, 800, 800, 950, 1310],
        ["l", "l", "c", "c", "c", "c", "c", "c", "c"], block_ends=ends, bold_cells=bold)
B.note(f"Hurricane vs. {NE.lower()} combined across the four cores (labels permuted within each core; statistic = sum of the per-core pseudo-F = {float(cb['sum_F']):.2f}): "
       f"p < 0.001 ({int(float(cb['perms'])):,} permutations); Fisher’s combined p of the four per-core tests = {float(cb['fisher_p']):.4f}.")

# ---- S6 dbRDA
db = rd("S_dbRDA.csv"); rows = []
for k in CORES:
    r = next(x for x in db if x["core"] == k)
    rows.append([k, r["n"], f"{float(r['constrained_pct']):.1f}", f"{float(r['CAP1_pct']):.1f}", f"{float(r['CAP2_pct']):.1f}", f2(r["F"]), f"{float(r['p']):.4f}", r["perms"]])
B.caption("Table S6.", "Distance-based redundancy analysis (dbRDA) of Bray–Curtis dissimilarity (square-root-transformed, median-normalized counts) constrained by deposit type. Constrained, percentage of total variance explained by deposit type; "
          "CAP1 and CAP2, percentage of total variance on each canonical axis (as plotted); the permutation test is an ANOVA-like test of the whole model (9,999 permutations).")
B.table([["Core", "n", "Constrained (%)", "CAP1 (%)", "CAP2 (%)", "Pseudo-F", "P(perm)", "Perms"]], rows, [1200, 800, 1700, 1200, 1200, 1200, 1100, 960], ["l", "c", "c", "c", "c", "c", "c", "c"])

# ---- figure S1 (carried over) with a caption derived from the figure's own title
sect.addprevious(fig_p); pf = Paragraph(fig_p, doc); pf.alignment = WD_ALIGN_PARAGRAPH.CENTER; para_fmt(pf, before=18, after=4, keep_next=True)
B.caption("Figure S1.", "Phylum-level community composition of the four cores (top five phyla per core; relative abundance, %). Phyla outside a core’s top five are grouped as “Others”.", keep_next=False)

# ---- S7 among-core composition / dispersion (transposed for portrait layout)
ac = {r["Level"]: r for r in rd("S_across_cores_finalgroups.csv")}
lv = ["Phylum", "Class", "ASV (exact sequence)", "Phylum (rarefied to 6040)", "Class (rarefied to 6040)", "ASV (rarefied to 6040)"]
g = lambda col, fn: [fn(ac[l][col]) for l in lv]
rows = [["No. of features"] + [f"{int(ac[l]['Features']):,}" for l in lv]]
merges = []
def sub(label): rows.append([label] + [""] * 6); merges.append((len(rows) - 1, 0, 6))
sub("Core effect, all samples (n = 64)")
rows += [["  Pseudo-F"] + g("Core_F", f2), ["  R²"] + g("Core_R2", f2), ["  p"] + g("Core_p", fp)]
sub("Dispersion among cores (PERMDISP)")
rows += [["  F"] + g("PERMDISP_F", f2), ["  p"] + g("PERMDISP_p", fp)]
sub("Core effect, storm deposits only (n = 38)")
rows += [["  Pseudo-F"] + g("StormOnly_F", f2), ["  R²"] + g("StormOnly_R2", f2), ["  p"] + g("StormOnly_p", fp)]
sub("Core × deposit-type interaction")
rows += [["  R²"] + g("CoreXStorm_R2", f2), ["  p"] + g("CoreXStorm_p", fp)]
B.caption("Table S7.", "Differences in microbial eDNA community composition and dispersion among the four Whale Beach cores (WB2501, WB2506, WB2508, WB2509). Cores were collapsed separately and placed on a common feature set before pooling (n = 64 samples): "
          "amplicon sequence variants (ASVs) were summed within phylum or class (ASVs unresolved at that rank were retained as one “Unassigned” category) or matched across cores by exact sequence (ASV level; absent ASVs set to zero). "
          "Bray–Curtis dissimilarities were calculated on square-root-transformed, median-normalized counts (counts scaled to the median library size, then square-root transformed). “Core effect” is a one-factor PERMANOVA (adonis2) with core as the grouping variable, tested on all samples and on hurricane and "
          f"{NE.lower()} overwash samples only (n = 38). PERMDISP (betadisper) tests whether cores differ in multivariate dispersion; non-significant values indicate that PERMANOVA differences reflect differences in community "
          f"composition rather than unequal spread. The core × deposit-type interaction (deposit types as defined in Table S2) tests whether the deposit-type signal differs among cores. Right-hand columns repeat each analysis "
          "after rarefying all samples to the smallest library (6,040 reads) to control for sequencing effort. All p-values are from 9,999 permutations (minimum attainable p = 0.0001, reported as < 0.001). "
          "Class level is the primary analysis; phylum and ASV levels are sensitivity analyses.")
B.table([["", "Phylum", "Class", "ASV", "Phylum", "Class", "ASV"]], rows, [2760, 1100, 1100, 1100, 1100, 1100, 1100], ["l"] + ["c"] * 6, merges=merges,
        group_header=[(0, 0, ""), (1, 3, "All reads"), (4, 6, "Rarefied to 6,040 reads")])
B.note("F, pseudo-F statistic; R², proportion of variance explained. Values are rounded.")

# ---- S8 pairwise cores (long format)
pc = rd("S_pairwise_cores.csv"); rows, ends = [], []
pairs = []
for r in pc:
    if r["pair"] not in pairs: pairs.append(r["pair"])
LV = [("Phylum", "Phylum"), ("Class", "Class"), ("ASV (exact sequence)", "ASV")]
for pr in pairs:
    for n, (l, ln) in enumerate(LV):
        r = next(x for x in pc if x["pair"] == pr and x["Level"] == l)
        rows.append([pr.replace(" vs ", " vs. ") if n == 0 else "", ln, f2(r["F"]), f2(r["R2"]), fp(r["p_holm"]), fp(r["PERMDISP_p"])])
    ends.append(len(rows) - 1)
B.caption("Table S8.", "Pairwise comparisons of microbial eDNA community composition and dispersion between cores (all samples; WB2501 n = 23, WB2506 n = 13, WB2508 n = 14, WB2509 n = 14). Each pair of cores was compared with a PERMANOVA "
          "(adonis2, 9,999 permutations) on Bray–Curtis dissimilarities of square-root-transformed, median-normalized counts, calculated on features pooled across cores as described in Table S7 (phylum, class, or exact-sequence ASV level). "
          "The six PERMANOVA p-values at each level were adjusted for multiple comparisons with the Holm method (p (Holm); minimum attainable value 0.0006, reported as < 0.001). PERMDISP p is the unadjusted permutation test "
          "(9,999 permutations) of whether the two cores differ in multivariate dispersion; values > 0.05 indicate that a significant PERMANOVA reflects a shift in community composition rather than unequal spread. F, pseudo-F statistic; R², proportion of variance explained.")
B.table([["Core comparison", "Feature level", "Pseudo-F", "R²", "p (Holm)", "PERMDISP p"]], rows, [2300, 1800, 1200, 1100, 1500, 1460], ["l", "l", "c", "c", "c", "c"], block_ends=ends)

# tiny closing paragraph so Word does not push an empty 12-pt paragraph onto a blank last page
endp = doc.add_paragraph(); set_run(endp.add_run(""), 1); para_fmt(endp, before=0, after=0, line=1.0)
pPr = endp._p.get_or_add_pPr(); mk = OxmlElement("w:rPr"); sz = OxmlElement("w:sz"); sz.set(qn("w:val"), "2"); mk.append(sz); pPr.append(mk)
doc.save(OUT); print("saved", OUT)

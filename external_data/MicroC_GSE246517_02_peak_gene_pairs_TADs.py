"""
=============================================================================
Micro-C TADs (Supplementary Fig S4B): are the 10-kb distal peak-gene pairs in the same domain?
  pairs : distal peaks (TSS == "X") and each gene in their "enhancer" column (10-kb assignment) of the peak table
  gene position: the TSS-associated peak (TSS == "O") whose matched_gene contains the gene; if several, the one closest
          to the distal peak; pairs without such a peak on the same chromosome are dropped
  class : domains (MicroC_GSE246517_01) containing the midpoint of the distal peak and the midpoint of the TSS peak.
          Same TAD / Different TADs; a midpoint on a boundary bin (or outside all domains) -> Unassignable (excluded)
  distance = between the two midpoints; bins <=2, 2-5, 5-10 kb

Input:  input_peak_annotations.xlsx (sheet 1), $ATAC_OUT/MicroC_GSE246517/output_domains_w25kb.bed
Output: $ATAC_OUT/MicroC_GSE246517/output_pairs_domain_classes.tsv, output_TAD_same_different_fraction.csv
Source: 261016_revision/19_MicroC_TAD/scripts/02_classify_pairs.py (column P_linked_w25) and
        11_TAD_different_fraction_figure.R (fractions); interval-based classes of the original are left out
=============================================================================
"""
import os
import numpy as np, pandas as pd

P = os.environ.get("ATAC_PROJECT", ".")
R = os.path.join(os.environ.get("ATAC_OUT", "results"), "MicroC_GSE246517")

t = pd.read_excel(os.path.join(P, "input_peak_annotations.xlsx"), sheet_name=0)
t["mid"] = ((t.start + t.end) // 2).astype(int)
split = lambda x: [g.strip() for g in str(x).replace(",", ";").split(";") if g.strip() not in ("", "-", "nan", "NA")]
g2tss = {}
for r in t[t.TSS == "O"].itertuples():
    for g in split(r.matched_gene):
        g2tss.setdefault(g, []).append((r.chr, r.start, r.end, r.mid, r.peak_id))
rows = []
for r in t[t.TSS == "X"].itertuples():
    for g in split(r.enhancer):
        c = [x for x in g2tss.get(g, []) if x[0] == r.chr]
        if not c: continue
        x = min(c, key=lambda x: abs(x[3] - r.mid))
        rows.append(dict(peak=r.peak_id, gene=g, chr=r.chr, e_start=r.start, e_end=r.end, e_mid=r.mid,
                         t_start=x[1], t_end=x[2], t_mid=x[3], tss_peak=x[4]))
pr = pd.DataFrame(rows)
pr["dist"] = (pr.t_mid - pr.e_mid).abs()
print("distal peak-gene pairs with a TSS peak:", len(pr))

d = pd.read_csv(os.path.join(R, "output_domains_w25kb.bed"), sep="\t", header=None, names=["chr", "s", "e"])
d["id"] = np.arange(len(d)); D = {c: x.sort_values("s") for c, x in d.groupby("chr")}
def dom_of_point(chrom, pos):
    x = D.get(chrom)
    if x is None: return -1
    i = np.searchsorted(x.s.values, pos, side="right") - 1
    return int(x.id.values[i]) if i >= 0 and pos < x.e.values[i] else -1
e = np.array([dom_of_point(r.chr, r.e_mid) for r in pr.itertuples()])
g = np.array([dom_of_point(r.chr, r.t_mid) for r in pr.itertuples()])
pr["TAD_class"] = np.where((e < 0) | (g < 0), "Unassignable", np.where(e == g, "Same TAD", "Different TADs"))
pr.to_csv(os.path.join(R, "output_pairs_domain_classes.tsv"), sep="\t", index=False)

pr["bin"] = pd.cut(pr.dist, [-1, 2000, 5000, np.inf], labels=["<=2", "2-5", "5-10"]).astype(str)
out = []
for b in ["<=2", "2-5", "5-10", "All"]:
    z = pr if b == "All" else pr[pr.bin == b]; v = z.TAD_class[z.TAD_class != "Unassignable"]
    out.append(dict(distance_kb=b, n_pairs=len(z), n_assignable=len(v), pct_same=100 * (v == "Same TAD").mean(),
                    pct_different=100 * (v == "Different TADs").mean()))
s = pd.DataFrame(out); s.to_csv(os.path.join(R, "output_TAD_same_different_fraction.csv"), index=False)
print(s.to_string(index=False, float_format=lambda x: f"{x:.1f}"))

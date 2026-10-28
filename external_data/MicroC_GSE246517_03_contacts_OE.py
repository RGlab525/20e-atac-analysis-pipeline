"""
=============================================================================
Micro-C contacts between distal peaks and their linked TSSs (Supplementary Fig S4C)
  matrix : 1-kb cis-only balanced matrix (MicroC_GSE246517_01)
  pairs  : output_pairs_domain_classes.tsv (MicroC_GSE246517_02); bins of the distal-peak midpoint and the TSS-peak midpoint.
           Kept: bins >= 2 apart (ignore_diags = 2); one row per (distal bin, TSS bin) combination; valid balanced
           values at the pair and at the mirror position (2 x TSS - distal peak, same distance on the other side of
           the TSS), as in the original analysis
  expected: mean balanced contact per chromosome and distance (cooltools 0.7.1 expected_cis, ignore_diags=2, no smoothing)
  O/E    : balanced contact / expected
  control: for each pair, the mean O/E of bins at the same distance (+/- 1 bin, both sides) from the same TSS bin that
           do not contain an ATAC-seq consensus peak (14,003) and are not the linked distal bin
  test   : two-sided Wilcoxon signed-rank test on log2((O/E pair + e) / (O/E control + e)), e = 1st percentile of O/E > 0;
           all pairs and 2-4, 4-6, 6-8, 8-10 kb
Input:  $ATAC_OUT/MicroC_GSE246517/output_S2_microC_1kb.cool, output_pairs_domain_classes.tsv; input_consensus_peaks.bed
Output: $ATAC_OUT/MicroC_GSE246517/output_pair_contacts_OE.tsv, output_OE_summary.tsv
Source: 261016_revision/19_MicroC_TAD/scripts/06_pair_contacts.py (pair selection) and 09_standard_controls.py (part A)
=============================================================================
"""
import os
import numpy as np, pandas as pd, cooler, cooltools
from scipy.stats import wilcoxon

P = os.environ.get("ATAC_PROJECT", ".")
R = os.path.join(os.environ.get("ATAC_OUT", "results"), "MicroC_GSE246517")
clr = cooler.Cooler(os.path.join(R, "output_S2_microC_1kb.cool")); bs = clr.binsize; size = clr.chromsizes
CHROMS = ["chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX"]
view = pd.DataFrame({"chrom": CHROMS, "start": 0, "end": [size[c] for c in CHROMS], "name": CHROMS})
mat = clr.matrix(balance=True)

# ---- pairs ----
pr = pd.read_csv(os.path.join(R, "output_pairs_domain_classes.tsv"), sep="\t")
pr["m_mid"] = 2 * pr.t_mid - pr.e_mid
pr["eb"] = pr.e_mid // bs; pr["tb"] = pr.t_mid // bs; pr["mb"] = pr.m_mid // bs
pr["sep"] = (pr.eb - pr.tb).abs()
pr = pr[(pr.sep >= 2) & (pr.m_mid >= 0) & (pr.m_mid < pr.chr.map(size))].copy()
pr = pr.drop_duplicates(["chr", "eb", "tb"])
vals = []
for c, z in pr.groupby("chr"):
    lo = int(min(z.eb.min(), z.tb.min(), z.mb.min())) * bs; hi = int(max(z.eb.max(), z.tb.max(), z.mb.max()) + 1) * bs
    for start in range(lo, hi, 400_000):          # 400-kb windows (pairs are <= 21 kb apart)
        sub = z[(z.tb * bs >= start) & (z.tb * bs < start + 400_000)]
        if sub.empty: continue
        w0 = max(0, start - 30_000); w1 = min(size[c], start + 430_000)
        m = mat.fetch(f"{c}:{w0}-{w1}"); off = w0 // bs
        for r in sub.itertuples():
            vals.append((r.Index, m[r.eb - off, r.tb - off], m[r.mb - off, r.tb - off]))
v = pd.DataFrame(vals, columns=["i", "linked", "mirror"]).set_index("i")
pc = pr.join(v)
pc = pc[np.isfinite(pc.linked) & np.isfinite(pc.mirror)].copy()
pc["dist_bin"] = pd.cut(pc.sep * bs, [1999, 4000, 6000, 8000, 11000], labels=["2-4", "4-6", "6-8", "8-10"])

# ---- O/E ----
exp = cooltools.expected_cis(clr, view_df=view, ignore_diags=2, smooth=False, clr_weight_name="weight")
E = {c: g.set_index("dist")["balanced.avg"] for c, g in exp.groupby("region1")}
pc["exp"] = [E[c].get(s, np.nan) for c, s in zip(pc.chr, pc.sep)]
pc["OE"] = pc.linked / pc.exp

# ---- TSS-anchored control bins ----
pk = pd.read_csv(os.path.join(P, "input_consensus_peaks.bed"), sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "s", "e"])
peakbins = {c: set(np.concatenate([np.arange(s // bs, (e - 1) // bs + 1) for s, e in zip(g.s, g.e)])) for c, g in pk.groupby("chr")}
ctrl = []
for c, z in pc.groupby("chr"):
    nb = size[c] // bs
    for start in range(0, int(size[c]), 400_000):
        sub = z[(z.tb * bs >= start) & (z.tb * bs < start + 400_000)]
        if sub.empty: continue
        w0 = max(0, start - 30_000); w1 = min(size[c], start + 430_000)
        m = mat.fetch(f"{c}:{w0}-{w1}"); off = w0 // bs
        for r in sub.itertuples():
            oe = []
            for side in (-1, 1):
                for ds in (r.sep - 1, r.sep, r.sep + 1):
                    if ds < 2: continue
                    b = r.tb + side * ds
                    if b < 0 or b >= nb or b == r.eb or b in peakbins.get(c, ()) or not (w0 // bs <= b < w1 // bs): continue
                    x = m[b - off, r.tb - off]; e = E[c].get(ds, np.nan)
                    if np.isfinite(x) and np.isfinite(e) and e > 0: oe.append(x / e)
            ctrl.append((r.Index, np.mean(oe) if oe else np.nan, len(oe)))
pc = pc.join(pd.DataFrame(ctrl, columns=["i", "OE_ctrl", "n_ctrl"]).set_index("i"))
pc.drop(columns=["mirror"]).to_csv(os.path.join(R, "output_pair_contacts_OE.tsv"), sep="\t", index=False)

# ---- test ----
eps = np.nanpercentile(pc.OE[pc.OE > 0], 1)
rows = []
def summ(label, z):
    z = z[np.isfinite(z.OE) & np.isfinite(z.OE_ctrl)]
    lr = np.log2((z.OE + eps) / (z.OE_ctrl + eps))
    rows.append(dict(set=label, n=len(z), median_OE=z.OE.median(), median_OE_ctrl=z.OE_ctrl.median(),
                     median_log2_pair_vs_ctrl=np.median(lr), P_wilcoxon=wilcoxon(lr).pvalue))
summ("all", pc)
for b in ["2-4", "4-6", "6-8", "8-10"]: summ(b + " kb", pc[pc.dist_bin == b])
s = pd.DataFrame(rows); s.to_csv(os.path.join(R, "output_OE_summary.tsv"), sep="\t", index=False)
print(s.to_string(index=False, float_format=lambda x: f"{x:.3g}"))

"""
=============================================================================
Micro-C (GSE246517, S2 cells): 1-kb matrix, balancing, insulation boundaries and domains (Supplementary Fig S4B, S4C)
  200-bp matrix -> 1-kb bins (cooler.coarsen_cooler, factor 5)
  cis-only iterative correction (cooler 0.10.4): ignore_diags=2, min_nnz=10, mad_max=5, tol=1e-5, max_iters=500
  diamond insulation (cooltools 0.7.1), 25-kb window, min_frac_valid_pixels=0.66, min_dist_bad_bin=2,
  boundaries = insulation minima with boundary strength above the Li threshold
  domain = interval between two consecutive boundaries
  chromosomes: chr2L, chr2R, chr3L, chr3R, chr4, chrX

Input:  public_data/GSE246517/input_microc_200bp.cool (MicroC_GSE246517_00_download.sh)
Output: $ATAC_OUT/MicroC_GSE246517/: output_S2_microC_1kb.cool (with balancing weights), output_insulation_1kb.tsv,
        output_boundaries_w25kb.bed, output_domains_w25kb.bed, output_domain_summary.tsv
Source: 261016_revision/19_MicroC_TAD/scripts/01_call_domains_cooltools.py (the original also ran 10- and 50-kb
        windows; each window is computed independently, so the 25-kb result is the same)
=============================================================================
"""
import os
import numpy as np, pandas as pd, cooler, cooltools

P = os.environ.get("ATAC_PROJECT", ".")
R = os.path.join(os.environ.get("ATAC_OUT", "results"), "MicroC_GSE246517"); os.makedirs(R, exist_ok=True)
SRC = os.path.join(P, "public_data/GSE246517/input_microc_200bp.cool")
C1 = os.path.join(R, "output_S2_microC_1kb.cool")
CHROMS = ["chr2L", "chr2R", "chr3L", "chr3R", "chr4", "chrX"]
WINDOW = 25000

if not os.path.exists(C1):
    cooler.coarsen_cooler(SRC, C1, factor=5, chunksize=10_000_000)
clr = cooler.Cooler(C1)
if "weight" not in clr.bins().columns:
    bias, stats = cooler.balance_cooler(clr, ignore_diags=2, min_nnz=10, mad_max=5, tol=1e-5, max_iters=500,
                                        cis_only=True, store=True)
    print("balance converged:", stats.get("converged"))
clr = cooler.Cooler(C1)

view = pd.DataFrame({"chrom": CHROMS, "start": 0, "end": [clr.chromsizes[c] for c in CHROMS], "name": CHROMS})
ins = cooltools.insulation(clr, [WINDOW], view_df=view, clr_weight_name="weight", min_frac_valid_pixels=0.66,
                           min_dist_bad_bin=2, threshold="Li", ignore_diags=2, verbose=False)
ins.to_csv(os.path.join(R, "output_insulation_1kb.tsv"), sep="\t", index=False)

b = ins[ins[f"is_boundary_{WINDOW}"]][["chrom", "start", "end", f"boundary_strength_{WINDOW}"]].copy()
b.to_csv(os.path.join(R, "output_boundaries_w25kb.bed"), sep="\t", header=False, index=False)
dom = []
for c in CHROMS:
    x = b[b.chrom == c].sort_values("start")
    for (_, s), (_, e) in zip(x.iloc[:-1].iterrows(), x.iloc[1:].iterrows()):
        dom.append((c, int(s.end), int(e.start)))
d = pd.DataFrame(dom, columns=["chrom", "start", "end"])
d = d[d.end > d.start]
d.to_csv(os.path.join(R, "output_domains_w25kb.bed"), sep="\t", header=False, index=False)
s = pd.DataFrame([dict(window_kb=WINDOW // 1000, n_boundaries=len(b), n_domains=len(d),
                       median_domain_kb=round((d.end - d.start).median() / 1000, 1),
                       genome_in_domains_pct=round(100 * (d.end - d.start).sum() / sum(clr.chromsizes[c] for c in CHROMS), 1))])
s.to_csv(os.path.join(R, "output_domain_summary.tsv"), sep="\t", index=False); print(s.to_string(index=False))

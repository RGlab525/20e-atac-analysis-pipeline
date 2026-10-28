"""
=============================================================================
R1 Major 2 (1 h RNA-seq), step 2: gene value = length-weighted mean coverage over the exon union, per sample
  GSE304853 RNA-seq: DMSO 1 h vs 20E 0.3 uM 1 h (control dsRNA), 2 replicates; authors' bigWig files
Input:  public_data/GSE304853/GSM91582{68,69,72,73}_*.bw; $ATAC_OUT/RNA_1h_GSE304853/output_gene_exons_union_dm6.bed (step 1)
Output: $ATAC_OUT/RNA_1h_GSE304853/output_RNA_1h_gene_exon_coverage.csv
Source: 261016_revision/13_R1_Major2_RNA_1h/scripts/02_gene_coverage_from_bigwig.py (paths only changed)
Requires: numpy, pyBigWig
=============================================================================
"""
import os, csv, numpy as np, pyBigWig
from collections import defaultdict
P = os.environ.get("ATAC_PROJECT", "."); O = os.path.join(os.environ.get("ATAC_OUT", "results"), "RNA_1h_GSE304853")
D = os.path.join(P, "public_data/GSE304853")
S = {"DMSO_1": "input_GSM9158268_1_RNA_DMSO_SKi_1.bw", "DMSO_2": "input_GSM9158269_2_RNA_DMSO_SKi_2.bw", "E20_1": "input_GSM9158272_5_RNA_20E1h_SKi_1.bw", "E20_2": "input_GSM9158273_6_RNA_20E1h_SKi_2.bw"}
bw = {k: pyBigWig.open(os.path.join(D, f)) for k, f in S.items()}
ex = defaultdict(list); sym = {}
for c, s, e, fb, sy in (l.rstrip("\n").split("\t") for l in open(os.path.join(O, "output_gene_exons_union_dm6.bed"))):
    ex[fb].append((c, int(s), int(e))); sym[fb] = sy
chroms = bw["DMSO_1"].chroms()
with open(os.path.join(O, "output_RNA_1h_gene_exon_coverage.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["fbgn", "symbol", "exon_bp"] + list(S))
    for fb, iv in ex.items():
        iv = [(c, s, e) for c, s, e in iv if c in chroms]
        if not iv: continue
        L = sum(e - s for _, s, e in iv); vals = []
        for k in S:
            tot = 0.0
            for c, s, e in iv:
                v = bw[k].stats(c, s, min(e, chroms[c]), type="mean")[0]; tot += (v or 0) * (e - s)
            vals.append(tot / L)
        w.writerow([fb, sym[fb], L] + [round(v, 5) for v in vals])
print("done")

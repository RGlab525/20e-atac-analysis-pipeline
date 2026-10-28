"""
=============================================================================
R3 Major 8: H3K27ac and EcR ChIP-seq signal (GSE304853; 20E 0.3 uM 1 h vs DMSO, 2 replicates) at every consensus ATAC peak
  Signal = mean of the authors' bigWig (dm6) over the peak +-500 bp (H3K27ac marks the nucleosomes flanking an open site)
  E1 = peak6200, P = peak6201 (CR43626 TSS peak), E2 = peak6202 (Fig 6B)

Input  (relative to ATAC_PROJECT): input_consensus_peaks.bed, public_data/GSE304853/*.bw (GSE304853_00_download.sh)
Output ($ATAC_OUT/GSE304853/): output_peaks_H3K27ac_EcR_20E_vs_DMSO.csv   (input of GSE304853_02_H3K27ac_EcR_summary.R)
Source: 261016_revision/12_R3_Major8_H3K27ac/scripts/01_H3K27ac_EcR_at_peaks.py (paths only changed)
Requires: python 3, numpy, pyBigWig
=============================================================================
"""
import os, csv, pyBigWig
P = os.environ.get("ATAC_PROJECT", "."); O = os.path.join(os.environ.get("ATAC_OUT", "results"), "GSE304853")
D = os.path.join(P, "public_data/GSE304853"); os.makedirs(O, exist_ok=True)
F = {"K27_DMSO": ["input_GSM9158292_25_ChIP_K27Ac_DMSO_SKi_1", "input_GSM9158293_26_ChIP_K27Ac_DMSO_SKi_2"],
     "K27_20E":  ["input_GSM9158296_29_ChIP_K27Ac_20E1h_SKi_1", "input_GSM9158297_31_ChIP_K27Ac_20E1h_SKi_2"],
     "EcR_DMSO": ["input_GSM9158284_17_ChIP_EcR_DMSO_SKi_1", "input_GSM9158285_18_ChIP_EcR_DMSO_SKi_2"],
     "EcR_20E":  ["input_GSM9158288_21_ChIP_EcR_20E1h_SKi_1", "input_GSM9158289_22_ChIP_EcR_20E1h_SKi_2"]}
bw = {k: [pyBigWig.open(os.path.join(D, f + ".bw")) for f in v] for k, v in F.items()}
peaks = [l.rstrip("\n").split("\t") for l in open(os.path.join(P, "input_consensus_peaks.bed"))]
FL = 500
out = []
for c, s, e, name, *_ in peaks:
    s, e = int(s), int(e); a, b = max(0, s - FL), e + FL
    row = {"peak": name, "chr": c, "start": s, "end": e}
    for k, lst in bw.items():
        for i, h in enumerate(lst, 1):
            try: v = h.stats(c, a, min(b, h.chroms()[c]), type="mean")[0]
            except Exception: v = None
            row[f"{k}_rep{i}"] = 0.0 if v is None else v
    out.append(row)
with open(os.path.join(O, "output_peaks_H3K27ac_EcR_20E_vs_DMSO.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
print("peaks", len(out))

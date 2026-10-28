"""
=============================================================================
R3 Major 8: FLAG-EcR / FLAG-Usp ChIP-seq of the original Fig 6B (GSE139316, Mazina et al. 2020; 20E 0.3 uM 1 h) at E1, P, E2
  Only the 20E 1 h condition exists (no untreated FLAG ChIP), so binding is expressed relative to input:
  enrichment = log2((ChIP + 0.1) / (input + 0.1)), mean over the peak. Tag-only control ChIP (FLAG-BioID2) treated the same way.
  Reference distribution: the 9,229 distal consensus peaks (percentile)

Input:  input_consensus_peaks.bed, public_data/GSE139316/*.bigwig (GSE304853_00_download.sh),
        $ATAC_OUT/GSE304853/output_peaks_H3K27ac_EcR_with_ATAC.csv (GSE304853_02_H3K27ac_EcR_summary.R; peak class TSS/distal)
Output: $ATAC_OUT/GSE304853/output_FLAG_EcR_Usp_enrichment_all_peaks.csv; E1/P/E2 values printed
Source: 261016_revision/12_R3_Major8_H3K27ac/scripts/04_FLAG_EcR_Usp_E1E2.py (paths only changed)
=============================================================================
"""
import os, csv, numpy as np, pyBigWig
P = os.environ.get("ATAC_PROJECT", "."); B = os.path.join(os.environ.get("ATAC_OUT", "results"), "GSE304853")
D = os.path.join(P, "public_data/GSE139316")
o = lambda f: pyBigWig.open(os.path.join(D, f))
EcR, Usp, Tag = o("input_GSM4138004_ChIP-FLAG_S2_line_FLAG-EcR_20E-1h.bigwig"), o("input_GSM4138007_ChIP-FLAG_S2_line_FLAG-Usp_20E-1h.bigwig"), o("input_GSM4138010_ChIP-FLAG_S2_line_FLAG-BioID2_20E-1h.bigwig")
InE, InU = o("input_GSM4138017_Input_S2_line_FLAG-EcR_20E-1h.bigwig"), o("input_GSM4138018_Input_S2_line_FLAG-Usp_20E-1h.bigwig")
cls = {r["peak"]: r["class"] for r in csv.DictReader(open(os.path.join(B, "output_peaks_H3K27ac_EcR_with_ATAC.csv")))}
def mean(h, c, a, b):
    v = h.stats(c, a, b, type="mean")[0]; return 0.0 if v is None else v
rows = []
for c, s, e, name, *_ in (l.rstrip("\n").split("\t") for l in open(os.path.join(P, "input_consensus_peaks.bed"))):
    s, e = int(s), int(e)
    if c not in EcR.chroms(): continue
    ie, iu = mean(InE, c, s, e), mean(InU, c, s, e)
    rows.append(dict(peak=name, cls=cls.get(name), EcR=np.log2((mean(EcR, c, s, e) + .1) / (ie + .1)), Usp=np.log2((mean(Usp, c, s, e) + .1) / (iu + .1)),
                     Tag=np.log2((mean(Tag, c, s, e) + .1) / (ie + .1))))
dist = [r for r in rows if r["cls"] == "distal"]
with open(os.path.join(B, "output_FLAG_EcR_Usp_enrichment_all_peaks.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
pc = lambda k, v: round(100 * np.mean([r[k] <= v for r in dist]), 1)
print("distal peaks: median enrichment EcR %.2f  Usp %.2f  tag-only %.2f" % tuple(np.median([r[k] for r in dist]) for k in ("EcR", "Usp", "Tag")))
for n, p in (("E1", "peak6200"), ("P", "peak6201"), ("E2", "peak6202")):
    r = next(x for x in rows if x["peak"] == p)
    print(f"{n} {p}: FLAG-EcR {r['EcR']:+.2f} (pct {pc('EcR', r['EcR'])})  FLAG-Usp {r['Usp']:+.2f} (pct {pc('Usp', r['Usp'])})  tag-only {r['Tag']:+.2f} (pct {pc('Tag', r['Tag'])})")

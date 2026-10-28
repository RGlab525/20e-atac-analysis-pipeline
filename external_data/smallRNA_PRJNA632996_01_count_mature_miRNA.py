"""
=============================================================================
R3 Major 8: mature miRNA counts from the public S2 small RNA-seq (PRJNA632996)
  1. cut at the 3' adapter (TruSeq, AGATCGG...); keep inserts >= 18 nt
  2. compare with miRBase dme mature miRNAs from the 2nd base on (the 1st base is read as N by the sequencer);
     a read is assigned when it matches from base 2 to (shorter length - 2) and the length difference is <= 3 nt
     (isomiRs with shorter/longer 3' ends or 1-2 nt 3' additions); several hits -> smallest length difference, ties dropped
Usage:  python3 -I external_data/smallRNA_PRJNA632996_01_count_mature_miRNA.py <FASTQ folder with mature.fa> [run ...]
Output: $ATAC_OUT/smallRNA_PRJNA632996/output_mature_miRNA_counts.tsv (miRNA x sample), output_read_summary.tsv
Source: 261016_revision/14_R3_Major8_miR282/scripts/01_count_mature_miRNA.py (output folder only changed)
=============================================================================
"""
import sys, os, gzip, collections
RAW = sys.argv[1]
O = os.path.join(os.environ.get("ATAC_OUT", "results"), "smallRNA_PRJNA632996"); os.makedirs(O, exist_ok=True)
SAMPLES = [("SRR11794166", "DMSO_6h_1"), ("SRR11794165", "DMSO_6h_2"), ("SRR11794164", "DMSO_6h_3"),
           ("SRR11794169", "20E_6h_1"), ("SRR11794168", "20E_6h_2"), ("SRR11794167", "20E_6h_3"),
           ("SRR11794159", "DMSO_2d_1"), ("SRR11794158", "DMSO_2d_2"), ("SRR11794157", "DMSO_2d_3"),
           ("SRR11794163", "20E_2d_1"), ("SRR11794162", "20E_2d_2"), ("SRR11794160", "20E_2d_3")]
if len(sys.argv) > 2: SAMPLES = [x for x in SAMPLES if x[0] in sys.argv[2:]]   # testing: a subset of runs only
mat, name = {}, None
for l in open(os.path.join(RAW, "input_mature.fa")):
    l = l.strip()
    if l.startswith(">"): name = l[1:].split()[0]
    elif name and name.startswith("dme-"): mat[name] = mat.get(name, "") + l.replace("U", "T")
K = 15
idx = collections.defaultdict(list)
for n, s in mat.items(): idx[s[1:1 + K]].append(n)
def assign(r):
    c = idx.get(r[1:1 + K])
    if not c: return None
    best, bd = [], 99
    for n in c:
        s = mat[n]; m = min(len(s), len(r))
        if abs(len(s) - len(r)) > 3 or s[1:m - 2] != r[1:m - 2]: continue   # last 2 nt may differ (3' non-templated additions)
        d = abs(len(s) - len(r))
        if d < bd: best, bd = [n], d
        elif d == bd: best.append(n)
    return best[0] if len(best) == 1 else None
cnt = {lab: collections.Counter() for _, lab in SAMPLES}; summ = []
for srr, lab in SAMPLES:
    tot = ins = hit = 0; cache = {}
    with gzip.open(os.path.join(RAW, "input_" + srr + ".fastq.gz"), "rt") as f:
        for i, l in enumerate(f):
            if i % 4 != 1: continue
            tot += 1; p = l.find("AGATCGG", 17)
            if p < 18: continue
            r = l[:p]; ins += 1
            if r not in cache: cache[r] = assign(r)
            m = cache[r]
            if m: cnt[lab][m] += 1; hit += 1
    summ.append((lab, srr, tot, ins, hit)); print(lab, tot, ins, hit, cnt[lab]["dme-miR-282-5p"], cnt[lab]["dme-miR-282-3p"], flush=True)
with open(os.path.join(O, "output_mature_miRNA_counts.tsv"), "w") as f:
    f.write("miRNA\t" + "\t".join(l for _, l in SAMPLES) + "\n")
    for n in sorted(mat): f.write(n + "\t" + "\t".join(str(cnt[l][n]) for _, l in SAMPLES) + "\n")
with open(os.path.join(O, "output_read_summary.tsv"), "w") as f:
    f.write("sample\trun\ttotal_reads\tinsert_ge18nt\tmature_miRNA_reads\n")
    for s in summ: f.write("\t".join(map(str, s)) + "\n")

#!/usr/bin/env python3
"""Nearest same-chromosome transcript TSS for each REDfly enhancer.

The input REDfly table is tab-delimited.
"""
import argparse
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_gtf(path):
    tss_by_gene = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            match = re.search(r'gene_name "([^"]+)"', fields[8])
            if match:
                tss = int(fields[3]) if fields[6] == "+" else int(fields[4])
                tss_by_gene[match.group(1)].append((fields[0], tss))
    return tss_by_gene


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_redfly_tsv", type=Path)
    parser.add_argument("input_dm6_gtf", type=Path)
    parser.add_argument("output_summary_tsv", type=Path)
    parser.add_argument("output_histogram_pdf", type=Path)
    parser.add_argument("output_histogram_counts_tsv", type=Path)
    args = parser.parse_args()
    for path in (args.input_redfly_tsv, args.input_dm6_gtf):
        if not path.is_file():
            parser.error(f"Input missing: {path}")
    redfly = pd.read_csv(args.input_redfly_tsv, sep="\t")
    required = ("Gene Name", "Chromosome", "Start", "End", "Gene ID")
    missing = set(required) - set(redfly.columns)
    if missing:
        parser.error("Missing REDfly columns: " + ", ".join(sorted(missing)))
    tss_by_gene = parse_gtf(args.input_dm6_gtf)
    rows = []
    distances = []
    for _, row in redfly.iterrows():
        gene = row["Gene Name"]
        chrom = row["Chromosome"]
        start = int(row["Start"])
        end = int(row["End"])
        center = (start + end) / 2
        candidate_tss = [tss for tss_chrom, tss in tss_by_gene.get(gene, ())
                         if tss_chrom == chrom]
        distance = min((abs(center - tss) for tss in candidate_tss), default=np.nan)
        if np.isfinite(distance):
            distances.append(distance)
        rows.append({"FBgn_ID": row["Gene ID"], "Gene_Symbol": gene,
                     "Enhancer_Chrom": chrom, "Enhancer_Start": start,
                     "Enhancer_End": end, "Distance_to_TSS": distance})
    for path in (args.output_summary_tsv, args.output_histogram_pdf,
                 args.output_histogram_counts_tsv):
        path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(args.output_summary_tsv, sep="\t", index=False,
                              na_rep="NA")
    if not distances:
        raise ValueError("No enhancer could be mapped to a same-chromosome transcript TSS")
    values = np.asarray(distances)
    histogram_values = values / 1000 if values.max() > 1000 else values
    units = "kb" if values.max() > 1000 else "bp"
    fig, ax = plt.subplots(figsize=(10, 6))
    counts, edges, _ = ax.hist(histogram_values, bins=50, color="gold", edgecolor="black", alpha=.7)
    ax.axvline(10 if units == "kb" else 10000, color="red", linestyle="--", label="10 kb")
    ax.set(xlabel=f"Distance to TSS ({units})", ylabel="Frequency",
           title="REDfly enhancer to nearest transcript TSS")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.output_histogram_pdf)
    plt.close(fig)
    pd.DataFrame({"bin_start": edges[:-1], "bin_end": edges[1:],
                  "count": counts.astype(int)}).to_csv(args.output_histogram_counts_tsv,
                                                         sep="\t", index=False)
    print(f"Mapped distances: {len(values)}; within 10 kb: {int((values <= 10000).sum())}")


if __name__ == "__main__":
    main()

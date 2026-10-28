#!/usr/bin/env python3
"""Extract per-region BPM bigWig signal in 10 bp bins around TSS/enhancer centers.

The input BED is expected to have chromosome/start/end in columns 1-3,
region name in column 4, and optionally strand in column 7.
"""
import argparse
import csv
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pyBigWig


def load_regions(path, stranded):
    regions = []
    with path.open() as handle:
        for lineno, line in enumerate(handle, 1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or not fields[0].startswith("chr"):
                continue
            start, end = int(fields[1]), int(fields[2])
            if start >= end:
                raise ValueError(f"Invalid BED interval at line {lineno}")
            strand = fields[6] if stranded and len(fields) > 6 else "+"
            name = fields[3] if len(fields) > 3 else f"{fields[0]}:{start}-{end}"
            regions.append((fields[0], (start + end) // 2, strand, name))
    if not regions:
        raise ValueError("No BED regions found")
    return regions


def one_sample(task):
    sample, bigwig_path, regions, half_window, bin_size = task
    result = np.full((len(regions), 2 * half_window // bin_size), np.nan)
    with pyBigWig.open(str(bigwig_path)) as bigwig:
        chromosome_sizes = bigwig.chroms()
        for index, (chromosome, center, strand, _) in enumerate(regions):
            start, end = center - half_window, center + half_window
            if chromosome not in chromosome_sizes or start < 0 or end > chromosome_sizes[chromosome]:
                continue
            values = np.asarray(bigwig.stats(chromosome, start, end, type="mean",
                                             nBins=result.shape[1]), dtype=float)
            values[np.isnan(values)] = 0.0
            result[index] = values[::-1] if strand == "-" else values
    return sample, result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("regions_bed", type=Path)
    parser.add_argument("input_manifest_tsv", type=Path, help="TSV with sample and bigwig columns")
    parser.add_argument("output_npz", type=Path)
    parser.add_argument("--stranded", action="store_true", help="Reverse BED column 7 minus-strand regions")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    if not args.regions_bed.is_file():
        parser.error(f"BED not found: {args.regions_bed}")
    regions = load_regions(args.regions_bed, args.stranded)
    if not args.input_manifest_tsv.is_file():
        parser.error(f"Manifest not found: {args.input_manifest_tsv}")
    with args.input_manifest_tsv.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    if not manifest or not {"sample", "bigwig"}.issubset(manifest[0]):
        parser.error("Manifest needs sample and bigwig columns")
    samples = [row["sample"] for row in manifest]
    if len(samples) != len(set(samples)):
        parser.error("Duplicate sample names in manifest")
    paths = {}
    for row in manifest:
        path = Path(row["bigwig"]).expanduser()
        if not path.is_absolute():
            path = args.input_manifest_tsv.resolve().parent / path
        paths[row["sample"]] = path
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        parser.error("Missing bigWig files: " + ", ".join(missing))
    tasks = [(sample, paths[sample], regions, 1500, 10) for sample in samples]
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        matrices = [matrix for _, matrix in executor.map(one_sample, tasks)]
    signal = np.stack(matrices, axis=1)
    valid = ~np.isnan(signal).any(axis=(1, 2))
    if not valid.any():
        raise ValueError("All regions fall outside the bigWig chromosome bounds")
    args.output_npz.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output_npz, mat=signal[valid],
                        names=np.asarray([region[3] for region in regions])[valid],
                        samples=np.asarray(samples))
    profile = signal[valid].mean(axis=0)
    midpoint = profile.shape[1] // 2
    for index, sample in enumerate(samples):
        print(f"{sample}\tn={valid.sum()}\tcenter={profile[index, midpoint-1:midpoint+1].mean():.2f}\tedge={profile[index, :10].mean():.2f}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Assemble the pre/new count matrix from two six-time-point featureCounts tables.

Column order in each input: 0h, 1h, 3h, 6h, 12h, 24h after the six annotation columns.
Counts are copied unchanged; no normalization or peak filtering is performed.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

TIMES = ['0h', '1h', '3h', '6h', '12h', '24h']
ANNOTATION = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']

def read_counts(path, batch):
    table = pd.read_csv(path, sep='\t', comment='#')
    if list(table.columns[:6]) != ANNOTATION or table.shape[1] != 12:
        raise ValueError('Expected six annotation columns and six time-point count columns')
    if table.Geneid.isna().any() or table.Geneid.duplicated().any():
        raise ValueError('Missing or duplicate peak IDs')
    table = table.set_index('Geneid')
    counts = table.iloc[:, 5:].copy()
    values = counts.to_numpy(dtype=float)
    if not np.isfinite(values).all() or (values < 0).any() or (values != np.floor(values)).any():
        raise ValueError('Counts must be finite nonnegative integers')
    counts.columns = [f'{batch}_{time}' for time in TIMES]
    return table.iloc[:, :5], counts

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input_pre_counts', type=Path)
    parser.add_argument('input_new_counts', type=Path)
    parser.add_argument('output_dir', type=Path)
    args = parser.parse_args()
    ann_pre, pre = read_counts(args.input_pre_counts, 'pre')
    ann_new, new = read_counts(args.input_new_counts, 'new')
    if set(pre.index) != set(new.index):
        raise ValueError('Peak ID sets differ between batches')
    if not ann_pre.equals(ann_new.loc[pre.index]):
        raise ValueError('Peak annotations differ between batches')
    result = pd.concat([pre, new.loc[pre.index]], axis=1)
    result.index.name = 'peak_id'
    metadata = pd.DataFrame([
        dict(SampleID=f'{batch}_{time}', sample=f's{i + 1}', time=time, batch=batch)
        for batch in ['pre', 'new'] for i, time in enumerate(TIMES)
    ])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_dir / 'output_counts.csv')
    metadata.to_csv(args.output_dir / 'output_metadata.tsv', sep='\t', index=False)
    print(f'{result.shape[0]} peaks x {result.shape[1]} samples')

if __name__ == '__main__':
    main()

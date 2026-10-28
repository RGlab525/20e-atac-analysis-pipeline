#!/usr/bin/env python3
"""Select representative TSSs, annotate peaks and join saved Wald statistics.

Transcript TSV columns: gene_id, transcript_id, length, chr, start, end, strand.
Transcript start/end are 1-based inclusive. The source reference defines length
as genomic transcript span (end-start+1). Longest per gene is retained; equal
lengths retain input order. BED peaks are 0-based half-open.

TSS class: overlap of a peak with a representative 1-bp TSS (BED convention).
TSS gene names: original R raw-coordinate inclusive overlap, restricted to TSS
class peaks. This preserves the original annotation's boundary behaviour.
Distal gene links: any representative TSS overlapping any consensus peak and
within a BEDTools-style 10-kb window of the distal peak (not centre distance).
Saved Wald values are joined by peak ID, without statistical recalculation.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

STATS = ['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj']


def read_wald(path):
    d = pd.read_csv(path)
    key = next((k for k in ['peak_id', 'PeakID', 'peak', 'Unnamed: 0'] if k in d), None)
    if key is None or not set(STATS).issubset(d):
        raise ValueError('Wald table needs peak ID and standard DESeq2 result columns')
    d = d.rename(columns={key: 'peak_id'})
    if d.peak_id.duplicated().any(): raise ValueError('Duplicate Wald peak IDs')
    return d[['peak_id'] + STATS].set_index('peak_id')


def annotate(transcripts, peaks):
    required = {'gene_id', 'transcript_id', 'length', 'chr', 'start', 'end', 'strand'}
    if not required.issubset(transcripts): raise ValueError('Missing transcript columns')
    if not transcripts.strand.isin(['+', '-']).all(): raise ValueError('Invalid strand')
    if not (transcripts['length'] == transcripts.end - transcripts.start + 1).all():
        raise ValueError('Expected genomic span length = end-start+1')
    if (transcripts.start < 1).any() or (transcripts.end < transcripts.start).any():
        raise ValueError('Invalid transcript coordinates')
    if (transcripts.gene_id == '').any(): raise ValueError('Missing gene ID')
    t = transcripts.sort_values('length', ascending=False, kind='stable').drop_duplicates('gene_id').copy()
    t['tss_start'] = np.where(t.strand == '+', t.start, t.end) - 1
    t['tss_end'] = t.tss_start + 1
    t = t.sort_values(['chr', 'tss_start', 'tss_end'], kind='stable').reset_index(drop=True)
    pk = peaks.copy()
    pk['TSS'] = 'X'; pk['matched_gene'] = '-'; pk['enhancer'] = '-'
    active = set()
    for ch, g in pk.groupby('chr', sort=False):
        z = t[t.chr == ch]; pos = z.tss_start.to_numpy()
        for i, row in g.iterrows():
            hit = z.iloc[np.searchsorted(pos, row.start):np.searchsorted(pos, row.end)]
            if len(hit):
                pk.at[i, 'TSS'] = 'O'; active.update(hit.index)
                # Preserve original GRanges matching on unconverted BED coordinates.
                gene_hits = z.iloc[np.searchsorted(pos, row.start - 1):np.searchsorted(pos, row.end, side='right')]
                pk.at[i, 'matched_gene'] = '; '.join(dict.fromkeys(gene_hits.gene_id))
    at = t.loc[sorted(active)].copy()
    links = []
    for ch, g in pk[pk.TSS == 'X'].groupby('chr', sort=False):
        z = at[at.chr == ch]; pos = z.tss_start.to_numpy()
        for i, row in g.iterrows():
            hit = z.iloc[np.searchsorted(pos, max(0, row.start - 10000)):np.searchsorted(pos, row.end + 10000)]
            if len(hit): pk.at[i, 'enhancer'] = ';'.join(dict.fromkeys(hit.gene_id))
            for h in hit.itertuples():
                links.append([row.chr, row.start, row.end, row.peak_id, h.chr, h.tss_start, h.tss_end, h.gene_id, h.transcript_id, '.', h.strand])
    return t, at, pk, links


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('input_transcripts_tsv', type=Path)
    p.add_argument('input_consensus_bed', type=Path)
    p.add_argument('output_dir', type=Path)
    p.add_argument('--wald-1h', type=Path)
    p.add_argument('--wald-3h', type=Path)
    a = p.parse_args()
    tx = pd.read_csv(a.input_transcripts_tsv, sep='\t', keep_default_na=False)
    peaks = pd.read_csv(a.input_consensus_bed, sep='\t', header=None, usecols=[0,1,2,3], names=['chr','start','end','peak_id'], keep_default_na=False)
    if peaks.peak_id.duplicated().any(): p.error('Duplicate peak IDs')
    if (peaks.start < 0).any() or (peaks.end <= peaks.start).any(): p.error('Invalid BED coordinates')
    t, active, ann, links = annotate(tx, peaks)
    if bool(a.wald_1h) != bool(a.wald_3h): p.error('Supply both --wald-1h and --wald-3h, or neither')
    a.output_dir.mkdir(parents=True, exist_ok=True)
    def bed(d, name, transcript=False):
        cols = [d.chr, d.tss_start, d.tss_end, d.gene_id]
        if transcript: cols.append(d.transcript_id)
        cols.extend([pd.Series('.', index=d.index), d.strand])
        pd.concat(cols, axis=1).to_csv(a.output_dir/name, sep='\t', index=False, header=False)
    bed(t, 'output_tss_representative.bed')
    bed(t, 'output_tss_representative_with_transcript.bed', True)
    bed(active, 'output_active_tss.bed', True)
    t.to_csv(a.output_dir/'output_longest_transcripts.tsv', sep='\t', index=False)
    pd.DataFrame(links).to_csv(a.output_dir/'output_enhancer_gene_links.tsv', sep='\t', header=False, index=False)
    ann.to_csv(a.output_dir/'output_peak_gene_annotations.csv', index=False)
    ann.loc[ann.TSS == 'O', ['peak_id', 'matched_gene']].to_csv(a.output_dir/'output_tss_peaks.csv', index=False)
    if a.wald_1h:
        w1, w3 = read_wald(a.wald_1h), read_wald(a.wald_3h)
        for w in [w1,w3]:
            if set(peaks.peak_id) != set(w.index): p.error('Wald and BED peak ID sets differ')
        one = ann.merge(w3, left_on='peak_id', right_index=True, validate='one_to_one', sort=False)
        both = ann.merge(w1.add_suffix('_1h'), left_on='peak_id', right_index=True, validate='one_to_one', sort=False)
        both = both.merge(w3.add_suffix('_3h'), left_on='peak_id', right_index=True, validate='one_to_one', sort=False)
        one.to_csv(a.output_dir/'output_annotated_wald_3h.csv', index=False)
        both.to_csv(a.output_dir/'output_annotated_wald_1h_3h.csv', index=False)
        with pd.ExcelWriter(a.output_dir/'output_peak_annotations.xlsx') as writer:
            one.to_excel(writer, sheet_name='Wald_3h', index=False)
            both.to_excel(writer, sheet_name='260322_1h_3h', index=False)
    print(f'Representative TSSs: {len(t)}; active: {len(active)}; classes: {ann.TSS.value_counts().to_dict()}')

if __name__ == '__main__': main()

#!/usr/bin/env bash
# MACS2 paired-end peak calling for ATAC-seq.
set -euo pipefail
if (( $# != 2 )); then
  echo 'Usage: call_peaks_bampe.sh input_sorted.bam output_dir' >&2
  exit 2
fi
[[ -f "$1" ]] || { echo 'Input BAM missing' >&2; exit 2; }
mkdir -p "$2"
macs2 callpeak -t "$1" -f BAMPE -g dm --nomodel --extsize 150   --keep-dup all --call-summits -p 0.01 -n output --outdir "$2"
# MACS2 writes output_peaks.narrowPeak and output_summits.bed in output_dir.

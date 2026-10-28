#!/usr/bin/env bash
# =============================================================================
# Analysis code: S.C.; selected source version retained.
# Inputs and outputs: generic input_* / output_* paths supplied as command-line arguments.
# Method and command examples: see RUNNING.md.
# =============================================================================
# MACS2 paired-end peak calling for ATAC-seq.
set -euo pipefail
if (( $# != 2 )); then
  echo 'Usage: 01_call_peaks_bampe.sh input_sorted.bam output_dir' >&2
  exit 2
fi
[[ -f "$1" ]] || { echo 'Input BAM missing' >&2; exit 2; }
mkdir -p "$2"
macs2 callpeak -t "$1" -f BAMPE -g dm --nomodel --extsize 150   --keep-dup all --call-summits -p 0.01 -n output --outdir "$2"
# MACS2 writes output_peaks.narrowPeak and output_summits.bed in output_dir.

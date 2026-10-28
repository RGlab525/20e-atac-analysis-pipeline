#!/usr/bin/env bash
# Single-end fastp trimming with a minimum read length of 25 bp.
set -euo pipefail
if (( $# != 2 )); then
  echo 'Usage: trim_rna_fastp.sh input_reads.fastq output_dir' >&2
  exit 2
fi
[[ -f "$1" ]] || { echo 'Input FASTQ missing' >&2; exit 2; }
mkdir -p "$2"
fastp -i "$1" -o "$2/output_trimmed.fastq.gz" -l 25   -j "$2/output_fastp.json" -h "$2/output_fastp.html"

#!/usr/bin/env bash
# Retain full narrowPeak rows whose chromosome, start, and end occur in a BED file.
set -euo pipefail

if (( $# != 3 )); then
  printf 'Usage: %s <input_matched.bed> <input_peaks.narrowPeak> <output_selected.narrowPeak>\n' "$0" >&2
  exit 2
fi
matched_bed=$1
narrowpeak=$2
output=$3
for file in "$matched_bed" "$narrowpeak"; do
  if [[ ! -f "$file" ]]; then
    printf 'Input not found: %s\n' "$file" >&2
    exit 1
  fi
done
mkdir -p "$(dirname "$output")"
awk 'BEGIN { FS = OFS = "\t" }
  NR == FNR { if (NF >= 3) wanted[$1 FS $2 FS $3] = 1; next }
  NF >= 3 && (($1 FS $2 FS $3) in wanted) { print }
' "$matched_bed" "$narrowpeak" > "$output"

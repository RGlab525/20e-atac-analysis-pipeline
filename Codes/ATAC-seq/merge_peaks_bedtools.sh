#!/usr/bin/env bash
# BEDTools coordinate union across six ATAC-seq time points.
# Optional score filtering is applied to each input before coordinate merging.
set -euo pipefail

if (( $# < 7 || $# > 8 )); then
  echo 'Usage: merge_peaks_bedtools.sh input_0h.narrowPeak input_1h.narrowPeak input_3h.narrowPeak input_6h.narrowPeak input_12h.narrowPeak input_24h.narrowPeak output_merged.bed [minimum_score]' >&2
  exit 2
fi
for file in "${@:1:6}"; do
  [[ -f "$file" ]] || { echo "Input missing: $file" >&2; exit 2; }
done
output=$7
apply_filter=0
minimum_score=0
if (( $# == 8 )); then
  [[ "$8" =~ ^[0-9]+([.][0-9]+)?$ ]] || { echo 'minimum_score must be numeric' >&2; exit 2; }
  apply_filter=1
  minimum_score=$8
fi
mkdir -p "$(dirname "$output")"
tmp=$(mktemp "$(dirname "$output")/.output_merged.XXXXXX")
trap 'rm -f "$tmp"' EXIT
awk -v filter="$apply_filter" -v threshold="$minimum_score" '
  BEGIN { OFS="\t" }
  NF >= 3 {
    if (filter && (NF < 5 || !($5 + 0 > threshold + 0))) next
    print $1, $2, $3
  }
' "${@:1:6}" \
  | LC_ALL=C sort -k1,1 -k2,2n \
  | bedtools merge -i - > "$tmp"
mv "$tmp" "$output"

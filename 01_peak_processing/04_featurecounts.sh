#!/usr/bin/env bash
# Count six time-point BAMs against the supplied consensus-peak SAF.
# Input: SAF, output path, then BAMs in 0h/1h/3h/6h/12h/24h order.
# Method: original featureCounts command; -F SAF -T 14 -p -M.
# Use separately for each replicate/batch to retain replicate-level columns.
set -euo pipefail
if (( $# != 8 )); then
  echo 'Usage: 04_featurecounts.sh input_peaks.saf output_readcounts.txt input_0h.bam input_1h.bam input_3h.bam input_6h.bam input_12h.bam input_24h.bam' >&2
  exit 2
fi
annotation=$1
output=$2
shift 2
mkdir -p "$(dirname "$output")"
featureCounts -F SAF -a "$annotation" -o "$output" -T 14 -p -M "$@"

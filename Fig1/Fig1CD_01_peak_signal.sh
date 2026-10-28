#!/usr/bin/env bash
# =============================================================================
# Fig 1C/D (and Supp Fig S6): mean RPGC signal over every consensus peak at each time point
#   deepTools multiBigwigSummary BED-file over the 14,003 consensus peaks (merged-replicate RPGC bigWigs)
# Input : output_atac_{t}.tss.bw (00_preprocessing/03), input_consensus_peaks.bed
# Output: known_genes/output_all_peaks_mean_signal.tab (+ .npz)
# Source: 261016_revision/08_QC_peak_validity_fragment/scripts/known_genes_01_peak_signal_deepTools.sh (J.Y.H.)
# Changes from the original code: paths via config/config.sh.
# Tools: deepTools v3.5.4
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
K=$(out_dir known_genes)

BWS=(); for t in $TIMES; do BWS+=("$OUT_ROOT/00_preprocessing/bigwig/output_atac_${t}.tss.bw"); done
multiBigwigSummary BED-file -b "${BWS[@]}" --BED "$PROJECT_DIR/input_consensus_peaks.bed" --labels $TIMES \
   -o "$K/output_all_peaks_signal.npz" --outRawCounts "$K/output_all_peaks_mean_signal.tab" -p 4

#!/bin/bash
# =============================================================================
# Fig 2A/B, step 2: deepTools matrices; Fig 2G signal per window
#   computeMatrix reference-point, +/-1.5 kb, bin 10 bp (default)
#     A: TSS of the 4,774 TSS peaks (--referencePoint TSS)
#     B: centre of the 9,229 distal peaks (--referencePoint center)
#   multiBigwigSummary BED-file: mean signal per window for the Fig 2G box plots (windows from Fig2G_01_boxplot_windows.R)
#   The published heatmaps and profiles are drawn from these matrices (drawing code not included).
#
# Input:
#   output_atac_{t}.tss.bw                  RPGC bigWig of the merged replicates (00_preprocessing)
#   $ATAC_OUT/Fig2AB_regions/output_A_TSS_peaks_4774.bed, output_B_distal_peaks_9229.bed      (Fig2AB_01_region_beds.R)
#   $ATAC_OUT/Fig2G_windows/windows_TSS_peaks.bed, windows_ENH_distal9229.bed    (Fig2G_01_boxplot_windows.R)
# Output ($ATAC_OUT/Fig2AB/): matrices/{t}.A_new_TSSpeaks.gz, matrices/{t}.B_new_distal9229.gz,
#   boxplot/{set}_mean_signal.tab
#
# Source: 261016_revision/03_Fig2AB_heatmap_boxplot/scripts/02_run_deepTools_v3.5.4.sh (deepTools v3.5.4).
#   computeMatrix options are those of the original Fig 2A/B commands (2025-03-10_ATAC_TSS_enrichment).
#   Only the region sets of the revised figure are kept here; the original sets (all UCSC TSSs, EnhancerAtlas)
#   were drawn the same way for comparison.
# Requires: deepTools v3.5.4 (set DT to its bin folder if it is not on PATH)
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
DT="${DT:+$DT/}"
O="$(out_dir Fig2AB)"; R="$OUT_ROOT/Fig2AB_regions"; W="$OUT_ROOT/Fig2G_windows"
M="$O/matrices"; mkdir -p "$M" "$O/boxplot"
${DT}computeMatrix --version

for t in $TIMES; do
  BW="$OUT_ROOT/00_preprocessing/bigwig/output_atac_${t}.tss.bw"
  ${DT}computeMatrix reference-point -S "$BW" -R "$R/output_A_TSS_peaks_4774.bed"    --referencePoint TSS    -b 1500 -a 1500 -o "$M/output_${t}_TSS_matrix.gz"   -p "$THREADS"
  ${DT}computeMatrix reference-point -S "$BW" -R "$R/output_B_distal_peaks_9229.bed" --referencePoint center -b 1500 -a 1500 -o "$M/output_${t}_distal_matrix.gz" -p "$THREADS"
done

# box-plot signal: mean bigWig value per window, all six time points in one table per region set
BWS=(); for t in $TIMES; do BWS+=("$OUT_ROOT/00_preprocessing/bigwig/output_atac_${t}.tss.bw"); done
for s in TSS_peaks ENH_distal9229; do
  ${DT}multiBigwigSummary BED-file -b "${BWS[@]}" --BED "$W/output_windows_${s}.bed" --labels $TIMES \
     -o "$O/boxplot/output_${s}.npz" --outRawCounts "$O/boxplot/output_${s}_mean_signal.tab" -p "$THREADS"
done

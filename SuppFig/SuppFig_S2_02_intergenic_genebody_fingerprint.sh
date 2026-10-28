#!/bin/bash
# =============================================================================
# Supplementary Fig. S2, step 2: deepTools for the intergenic (S2A), gene-body (revision) and fingerprint (S2B) panels
#   S2A intergenic : computeMatrix scale-regions over the 3-kb regions (body 3,000 bp, bin 10 bp, --skipZeros),
#                    all six bigWigs in one matrix
#   Gene body      : the same command on TSS +1 kb to TES (revision, R2-4)
#   S2B fingerprint: plotFingerprint on the merged BAMs (default settings), cumulative coverage values only
#                    (--outRawCounts; the plot file is not written)
#
# Input:
#   output_atac_{t}.tss.bw, input_bam/output_atac_{t}.merge.sorted.bam
#   $ATAC_OUT/SuppFig_S2_regions/output_intergenic_down3000_unique.bed, output_genebody_TSSplus1kb_to_TES_unique.bed
# Output ($ATAC_OUT/SuppFig_S2/): output_Intergenic_S2A_down3000.mat.gz, output_GeneBody_TSSplus1kb_TES.mat.gz,
#   output_fingerprint_counts.tab
#
# Source: 2025-03-10_ATAC_TSS_enrichment (lines 908-909 plotFingerprint, 1286-1297 S2A computeMatrix);
#   261016_revision/03_Fig2AB_heatmap_boxplot/scripts/02_run_deepTools_v3.5.4.sh (gene body)
# Changes from the original code: --labels and --outRawCounts added to plotFingerprint (values unchanged); no plot file
# Requires: deepTools v3.5.4 (set DT to its bin folder if it is not on PATH)
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
DT="${DT:+$DT/}"
O="$(out_dir SuppFig_S2)"; R="$OUT_ROOT/SuppFig_S2_regions"
BWS=(); BAMS=(); for t in $TIMES; do BWS+=("$OUT_ROOT/00_preprocessing/bigwig/output_atac_${t}.tss.bw"); BAMS+=("$OUT_ROOT/00_preprocessing/bam/output_atac_${t}.merge.sorted.bam"); done

# S2A: intergenic
${DT}computeMatrix scale-regions -S "${BWS[@]}" -R "$R/output_intergenic_down3000_unique.bed" \
   --beforeRegionStartLength 0 --regionBodyLength 3000 --afterRegionStartLength 0 --skipZeros \
   -o "$O/output_Intergenic_S2A_down3000.mat.gz" -p "$THREADS"

# gene body (revision panel)
${DT}computeMatrix scale-regions -S "${BWS[@]}" -R "$R/output_genebody_TSSplus1kb_to_TES_unique.bed" \
   --beforeRegionStartLength 0 --regionBodyLength 3000 --afterRegionStartLength 0 --skipZeros \
   -o "$O/output_GeneBody_TSSplus1kb_TES.mat.gz" -p "$THREADS"

# S2B: fingerprint (cumulative coverage)
${DT}plotFingerprint -b "${BAMS[@]}" --labels $TIMES \
   --outRawCounts "$O/output_fingerprint_counts.tab" --numberOfProcessors "$THREADS"

#!/usr/bin/env bash
# =============================================================================
# Merge the two replicate libraries of each time point
#   Input : input_bam/input_atac_rep1_{t}.sorted.bam (rep1), input_bam/input_atac_rep2_{t}.sorted.bam (rep2)
#   Output: output_atac_{t}.merge.sorted.bam (+ .bai)
#   Used for peak calling (01_peak_processing), coverage tracks (03), Circos (Fig1B) and fragment-size plots (Supp Fig S1).
# Source: command recorded in the @PG header of input_bam/output_atac_{t}.merge.sorted.bam
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
OUT=$(out_dir 00_preprocessing/bam)

for t in $TIMES; do
    samtools merge "$OUT/output_atac_${t}.merge.sorted.bam" \
        "$PROJECT_DIR/input_bam/input_atac_rep1_${t}.sorted.bam" "$PROJECT_DIR/input_bam/input_atac_rep2_${t}.sorted.bam"
    samtools index "$OUT/output_atac_${t}.merge.sorted.bam"
done

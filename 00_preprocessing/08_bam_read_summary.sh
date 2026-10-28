#!/usr/bin/env bash
# =============================================================================
# Read summary of the filtered analysis BAMs (all *.sorted.bam in a folder)
#   Total reads, properly paired fragments, reads with MAPQ >= 30.
#   "Mapped fragments" in Supplementary Fig S1C = Total_Reads / 2 of these BAMs
#   (checked for all 12 libraries).
# Source: 2025-03-05_read_sammary_alltogether.sh (J.Y.H.); same commands, comments translated
# Changes from the original code: paths via config/config.sh; a zero multi-mapped count is written as 0
#   (the original wrote an empty field).
# Note: the BAMs are already filtered (trimming, Bowtie2 and filtering by S.C.), so the "unaligned" and
#       "overall alignment rate" columns do not describe the Bowtie2 alignment rate.
# Input : input_bam/*.sorted.bam (filtered BAMs of the 12 libraries)
# Output: $ATAC_OUT/00_preprocessing/qc/output_atac_combined_stats.tsv
# Usage: 08_bam_read_summary.sh [bam_folder]   (default: $PROJECT_DIR/input_bam)
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
BAM_DIR=${1:-$PROJECT_DIR/input_bam}
OUT=$(out_dir 00_preprocessing/qc)
STATS="$OUT/output_atac_combined_stats.tsv"

echo -e "Sample\tTotal_Reads\tConcordantly_Aligned_Reads\tConcordantly_Aligned_Percentage\tMulti_Mapped_Reads\tUnaligned_Reads\tOverall_Alignment_Rate\tProperly_Paired_Fragments\tProperly_Paired_Percentage\tQualified_Fragments\tFinal_Filtered_Reads" > "$STATS"

for BAM_FILE in "$BAM_DIR"/*.sorted.bam; do
    PREFIX=$(basename "${BAM_FILE%.sorted.bam}")
    echo "Processing $PREFIX..."
    TOTAL_READS=$(samtools view -c "$BAM_FILE")
    CONCORDANT=$(samtools view -c -F 4 -f 2 -q 1 "$BAM_FILE")
    CONCORDANT_PCT=$(awk "BEGIN {printf \"%.1f\", ($CONCORDANT / $TOTAL_READS) * 100}")
    MULTI=$(samtools view -F 4 "$BAM_FILE" | awk '$5 == 0 {n++} END {print n + 0}')
    UNALIGNED=$(samtools view -c -f 4 "$BAM_FILE")
    OVERALL=$(awk "BEGIN {printf \"%.1f\", (($TOTAL_READS - $UNALIGNED) / $TOTAL_READS) * 100}")
    PP_FRAG=$(( $(samtools view -c -f 2 "$BAM_FILE") / 2 ))
    PP_PCT=$(awk "BEGIN {printf \"%.1f\", ($PP_FRAG / ($TOTAL_READS / 2)) * 100}")
    # mitochondrial chromosome name, as in the original script (looks for "chrM" or "MT").
    # The analysis BAMs name it "M", so no chrM reads are excluded from the last two columns.
    CHR_NAMES=$(samtools view -H "$BAM_FILE" | grep "^@SQ" | cut -f 2 | cut -d':' -f2)
    if echo "$CHR_NAMES" | grep -q "chrM"; then MITO="chrM"
    elif echo "$CHR_NAMES" | grep -q "MT"; then MITO="MT"
    else MITO=""; fi
    if [ -n "$MITO" ]; then
        QUAL=$(samtools view -F 1024 -q 30 "$BAM_FILE" | grep -v "$MITO" | grep -c "^")
    else
        QUAL=$(samtools view -F 1024 -q 30 -c "$BAM_FILE")
    fi
    echo -e "${PREFIX}\t${TOTAL_READS}\t${CONCORDANT}\t${CONCORDANT_PCT}\t${MULTI}\t${UNALIGNED}\t${OVERALL}\t${PP_FRAG}\t${PP_PCT}\t$((QUAL / 2))\t${QUAL}" >> "$STATS"
done
column -t -s $'\t' "$STATS"

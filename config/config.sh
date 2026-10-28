# Shared settings for the shell scripts.
#   ATAC_PROJECT : folder that holds the input files (relative names in the scripts are resolved here)
#   ATAC_OUT     : output folder
#   THREADS      : number of threads
# Usage: source "$(dirname "$0")/../config/config.sh"

PROJECT_DIR="${ATAC_PROJECT:-.}"
OUT_ROOT="${ATAC_OUT:-results}"
THREADS="${THREADS:-8}"

TIMES="0h 1h 3h 6h 12h 24h"
GENOME_FA="${GENOME_FA:-$PROJECT_DIR/reference/input_dm6.fa}"
BOWTIE2_INDEX="${BOWTIE2_INDEX:-$PROJECT_DIR/reference/input_bowtie2_index}"
CHROM_SIZES="${CHROM_SIZES:-$PROJECT_DIR/input_dm6.chrom.sizes}"
EFFECTIVE_GENOME_SIZE=142573017   # dm6, used for RPGC normalization

out_dir() { mkdir -p "$OUT_ROOT/$1"; echo "$OUT_ROOT/$1"; }

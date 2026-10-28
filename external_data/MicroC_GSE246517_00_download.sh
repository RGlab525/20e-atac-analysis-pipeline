#!/bin/bash
# =============================================================================
# Download the pooled S2 Micro-C contact matrix (Brown et al. 2024, GSE246517; 200-bp bins, dm6)
# Output: $ATAC_PROJECT/public_data/GSE246517/GSE246517_S2.cool
# Source: file downloaded from GEO on 2026-10-28 for 261016_revision/19_MicroC_TAD
# =============================================================================
set -e
source "$(dirname "$0")/../config/config.sh"
mkdir -p "$PROJECT_DIR/public_data/GSE246517"; cd "$PROJECT_DIR/public_data/GSE246517"
curl -sSf -o input_microc_200bp.cool.gz https://ftp.ncbi.nlm.nih.gov/geo/series/GSE246nnn/GSE246517/suppl/GSE246517_S2.cool.gz
gunzip -k input_microc_200bp.cool.gz
ls -l

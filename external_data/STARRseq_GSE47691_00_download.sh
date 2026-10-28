#!/bin/bash
# =============================================================================
# Download Shlyueva et al. 2014 (GSE47691) S2 STARR-seq / DHS-seq peak files (dm3; untreated vs 20E 41 uM 24 h, 2 replicates)
#   and the UCSC dm3 -> dm6 chain (also used by the Hi-C TAD analysis)
# Output: $ATAC_PROJECT/public_data/GSE47691_STARR_DHS/
# Source: 261016_revision/15_R3_STARRseq_E1E2/scripts/00_download.sh (output folder only changed)
# =============================================================================
set -e
source "$(dirname "$0")/../config/config.sh"
mkdir -p "$PROJECT_DIR/public_data/GSE47691_STARR_DHS"; cd "$PROJECT_DIR/public_data/GSE47691_STARR_DHS"
base=https://ftp.ncbi.nlm.nih.gov/geo/samples
for f in GSM1357045_S2_STARRseq_no_ecd_rep1 GSM1357046_S2_STARRseq_no_ecd_rep2 \
         GSM1357047_S2_STARRseq_with_ecd_rep1 GSM1357048_S2_STARRseq_with_ecd_rep2 \
         GSM1154837_S2_DHSseq_no_ecd_rep1 GSM1154838_S2_DHSseq_no_ecd_rep2 \
         GSM1154839_S2_DHSseq_with_ecd_rep1 GSM1154840_S2_DHSseq_with_ecd_rep2; do
  g=${f%%_*}
  curl -sSf -o input_${f}.peaks.txt.gz "$base/${g:0:7}nnn/$g/suppl/${f}.peaks.txt.gz"
done
curl -sSf -o input_dm3ToDm6.over.chain.gz https://hgdownload.soe.ucsc.edu/goldenPath/dm3/liftOver/dm3ToDm6.over.chain.gz
ls -l

#!/bin/bash
# =============================================================================
# Download the authors' bigWig files used in the revision analyses (R1 Major 2, R3 Major 8)
#   GSE304853 (Evdokimova et al. 2025; S2 cells, DMSO vs 20E 0.3 uM for 1 h, control dsRNA "SKi", 2 replicates)
#     RNA-seq, EcR ChIP-seq, H3K27ac ChIP-seq
#   GSE139316 (Mazina et al. 2020; S2 FLAG-EcR / FLAG-Usp / FLAG-BioID2 lines, 20E 1 h): ChIP and input bigWigs
# Output: $ATAC_PROJECT/public_data/GSE304853/, $ATAC_PROJECT/public_data/GSE139316/
# Source: file names as used in 261016_revision/12_R3_Major8_H3K27ac/public_data (downloaded from GEO 2026-10).
#   The GEO supplementary-file URL pattern below is the standard one; check the GEO sample pages if a link fails.
# =============================================================================
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
geo() {  # geo <GSM id> <file name> <output folder>
    local g=$1 f=$2 d=$3; mkdir -p "$d"
    curl -sSf -o "$d/input_$f.gz" "https://ftp.ncbi.nlm.nih.gov/geo/samples/${g:0:7}nnn/$g/suppl/$f.gz" && gunzip -f "$d/input_$f.gz" ||
    curl -sSf -o "$d/input_$f" "https://ftp.ncbi.nlm.nih.gov/geo/samples/${g:0:7}nnn/$g/suppl/$f"
}
D1="$PROJECT_DIR/public_data/GSE304853"
for f in GSM9158268_1_RNA_DMSO_SKi_1 GSM9158269_2_RNA_DMSO_SKi_2 GSM9158272_5_RNA_20E1h_SKi_1 GSM9158273_6_RNA_20E1h_SKi_2 \
         GSM9158284_17_ChIP_EcR_DMSO_SKi_1 GSM9158285_18_ChIP_EcR_DMSO_SKi_2 GSM9158288_21_ChIP_EcR_20E1h_SKi_1 GSM9158289_22_ChIP_EcR_20E1h_SKi_2 \
         GSM9158292_25_ChIP_K27Ac_DMSO_SKi_1 GSM9158293_26_ChIP_K27Ac_DMSO_SKi_2 GSM9158296_29_ChIP_K27Ac_20E1h_SKi_1 GSM9158297_31_ChIP_K27Ac_20E1h_SKi_2; do
    geo "${f%%_*}" "$f.bw" "$D1"
done
D2="$PROJECT_DIR/public_data/GSE139316"
for f in GSM4138004_ChIP-FLAG_S2_line_FLAG-EcR_20E-1h GSM4138007_ChIP-FLAG_S2_line_FLAG-Usp_20E-1h GSM4138010_ChIP-FLAG_S2_line_FLAG-BioID2_20E-1h \
         GSM4138017_Input_S2_line_FLAG-EcR_20E-1h GSM4138018_Input_S2_line_FLAG-Usp_20E-1h; do
    geo "${f%%_*}" "$f.bigwig" "$D2"
done
ls -l "$D1" "$D2"

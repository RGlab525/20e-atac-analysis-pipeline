#!/bin/bash
# =============================================================================
# R3 Major 8: download the public S2 small RNA-seq (Jin et al. 2020 RNA Biol, PRJNA632996; DMSO vs 20E 5 uM, 6 h and 2 d,
#   3 replicates; 12 runs, about 6 GB) from ENA, and miRBase mature.fa
# Usage: smallRNA_PRJNA632996_00_download.sh <download folder>
# Source: 261016_revision/14_R3_Major8_miR282/scripts/00_download.sh
# =============================================================================
set -euo pipefail
mkdir -p "${1:?input_download_dir}"; cd "$1"
for r in SRR11794157 SRR11794158 SRR11794159 SRR11794160 SRR11794162 SRR11794163 SRR11794164 SRR11794165 SRR11794166 SRR11794167 SRR11794168 SRR11794169; do
  echo "$r"
done | xargs -P 12 -n 1 bash -c 'r=$1; curl -s -C - --retry 5 -o "input_$r.fastq.gz" "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR117/0${r: -2}/$r/$r.fastq.gz"' _
curl -s -L -o input_mature.fa https://www.mirbase.org/download/mature.fa
ls -la

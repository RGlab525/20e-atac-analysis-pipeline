#!/usr/bin/env bash
# HOMER known and de novo motif enrichment in TSS / enhancer DARs.
#
# Manuscript: Fig 3E (TSS), Fig 4E (enhancer); Methods "Motif enrichment analysis" (HOMER v5.1)
# Input : $ATAC_OUT/motif_beds/{TSS,enhancer}_3h_{up,down}.bed (Fig3E_4E_motif_input_beds.R)
# Output: $ATAC_OUT/HOMER/{TSS,enhancer}_{up,down}/ (homerResults.html, knownResults.txt, ...)
# Source: the command recorded by HOMER in motifFindingParameters.txt of the runs used for the figures, e.g.
#         "260610_enh_3h_down.bed dm6 .../enhancer_down -size given -mask -p 30"
#         (default background: GC- and length-matched random genomic regions).
set -euo pipefail
source "$(dirname "$0")/../config/config.sh"
IN="$OUT_ROOT/motif_beds"
OUT=$(out_dir HOMER)

for cls in TSS enhancer; do
    for dir in up down; do
        findMotifsGenome.pl "$IN/output_${cls}_3h_${dir}.bed" dm6 "$OUT/${cls}_${dir}" -size given -mask -p "$THREADS"
    done
done

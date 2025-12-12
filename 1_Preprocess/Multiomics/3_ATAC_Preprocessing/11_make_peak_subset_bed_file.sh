#!/usr/bin/env bash
PATH_IN=/jukebox/krienen/victor/marm_hmba_atac_50k/data/marm_subcortex_filtered_peaks.txt
PATH_OUT=/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/data/atac/atac_462k_peaks.bed
awk -F '[:\\-]' 'BEGIN{OFS="\t"} {print $1, $2, $3}' "$PATH_IN" > "$PATH_OUT"

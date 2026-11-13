#!/usr/bin/env bash

PDIR="/jukebox/krienen/lab_members/mike/251007_HMBA_Spatial_ATAC/DESeq2_AllRegionNT_RNA"
BATCH="${PDIR}/pipeline/batchscript_region_nt_deseq2.sh"
REGION_NT_SET="${PDIR}/data/all_region_nt.txt"
LOGDIR="${PDIR}/DESeq2_logs"

mkdir -p "$LOGDIR"

while IFS=_ read -r region nt; do
  [[ -z "$region" || -z "$nt" ]] && continue
  sbatch \
    --job-name="deseq_${region}_${nt}" \
    --output="${LOGDIR}/${region}_${nt}_%j.out" \
    "$BATCH" "$region" "$nt"
done < "$REGION_NT_SET"

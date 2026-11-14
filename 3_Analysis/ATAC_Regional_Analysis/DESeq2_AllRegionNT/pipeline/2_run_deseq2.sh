#!/usr/bin/env bash

REPO="$(git rev-parse --show-toplevel)"  # root of this repository (run from anywhere inside it)
PDIR="${REPO}/3_Analysis/ATAC_Regional_Analysis/DESeq2_AllRegionNT"
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

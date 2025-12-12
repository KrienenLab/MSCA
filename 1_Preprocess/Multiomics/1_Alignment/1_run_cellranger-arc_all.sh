#!/usr/bin/env bash
cd $(dirname $0)
cat ./sample_names.txt | while read SAMPLE; do
    cd "$SAMPLE"
    if [ ! -e "${SAMPLE}/outs" ]; then
        [ ! -e "${SAMPLE}/_lock" ] || rm "${SAMPLE}/_lock"
        sbatch -J "cellranger-arc_${SAMPLE}" ../wrapper_slurm_cellranger-arc.sh "$SAMPLE"
    fi
    cd ..
done

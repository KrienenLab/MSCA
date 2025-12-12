#!/usr/bin/env bash

# Run this from the parent directory containing all sample subfolders
cd $(dirname $0)

module purge
module load anacondapy/2024.02
conda activate cellbender

cat ./sample_names.txt | while read SAMPLE; do
    [ -e "$SAMPLE/QC_FAIL" ] && continue
    cd "$SAMPLE"
    if [ ! -e "cellbender_out/cellbender_output_filtered.h5" ]; then
        sbatch -J "cellbender_${SAMPLE}" ../wrapper_slurm_cellbender.sh "$SAMPLE"
    fi
    cd ..
done
exit 0

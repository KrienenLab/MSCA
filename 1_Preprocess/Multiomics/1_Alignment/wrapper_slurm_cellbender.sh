#!/usr/bin/env bash
#SBATCH -o cellb_%j.out
#SBATCH -p all
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --gres=gpu:1

# Usage:
# Run this script from within the sample directory, e.g.:
# cd .../alignments/.../THE_SAMPLE_NAME
# sbatch ../wrapper_slurm_cellbender.sh THE_SAMPLE_NAME

module purge
module load anacondapy/2024.02
conda activate cellbender

SAMPLE=$1

[ -e cellbender_out ] || mkdir cellbender_out
cd cellbender_out

cellbender remove-background \
--input ../${SAMPLE}/outs/raw_feature_bc_matrix.h5 \
--output cellbender_output.h5 \
--exclude-feature-types Peaks \
--cuda

# Optional repack for Seurat
PATH_H5="cellbender_output_filtered.h5"
if [ -e "$PATH_H5" ]; then
    PATH_SEUR_H5="cellbender_output_filtered_seurat.h5"
    [ -e "$PATH_SEUR_H5" ] || ptrepack --complevel 5 "${PATH_H5}:/matrix" "${PATH_SEUR_H5}:/matrix"
fi

exit 0
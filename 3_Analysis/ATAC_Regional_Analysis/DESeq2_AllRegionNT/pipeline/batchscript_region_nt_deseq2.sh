#!/usr/bin/env bash
#SBATCH -o deseq_%j.out
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --cpus-per-task=1        # 
#SBATCH --mem=300G               # 
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends

source $(conda info --base)/etc/profile.d/conda.sh
conda activate /jukebox/scratch/md6347/envs/r4.4

REPO="$(git rev-parse --show-toplevel)"  # root of this repository (run from anywhere inside it)
PATH_SCRIPT="${REPO}/3_Analysis/Nonneuron/FigS5H_microglia_DEG/RunDESeq2apeglm_bulked_h5ad.R"
PATH_H5AD="${REPO}/3_Analysis/ATAC_Regional_Analysis/common_data/Cluster_v4_by_462k_peaks_modcopy.h5ad"
PATH_OUT="${REPO}/3_Analysis/ATAC_Regional_Analysis/DESeq2_AllRegionNT/DESeq2_results"
REGION=$1
NT=$2

# Rscript RunDESeq2apeglm_bulked_h5ad.R <pseudobulk.h5ad> 
# <group_column_in_obs> 
# <numerator_group> 
# <denominator_group> 
# <out_dir> 
# [counts_layer]

Rscript "$PATH_SCRIPT" "$PATH_H5AD" \
"${REGION}_vs_OtherNeur" \
"${REGION}_${NT}" \
"OtherNeur" \
"${PATH_OUT}/${REGION}_${NT}"

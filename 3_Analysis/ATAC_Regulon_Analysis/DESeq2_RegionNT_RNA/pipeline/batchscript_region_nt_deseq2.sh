#!/usr/bin/env bash
#SBATCH -o deseq_%j.out
#SBATCH --time=24:00:00          # total run time limit (HH:MM:SS)
#SBATCH --cpus-per-task=1        # 
#SBATCH --mem=30G                # 
#SBATCH --mail-type=begin        # send email when job begins
#SBATCH --mail-type=end          # send email when job ends

source $(conda info --base)/etc/profile.d/conda.sh
conda activate /jukebox/scratch/md6347/envs/r4.4

PATH_SCRIPT=/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/analysis/code/util/RunDESeq2apeglm_bulked_h5ad.R
PATH_H5AD=/jukebox/krienen/lab_members/mike/251022_HMBA_Spatial_ATAC/DESeq2_RegionNT_RNA/data/Cluster_v4_by_genes_modcopy.h5ad
PATH_OUT=/jukebox/krienen/lab_members/mike/251022_HMBA_Spatial_ATAC/DESeq2_RegionNT_RNA/DESeq2_results
REGION=$1
NT=$2

# Rscript RunDESeq2apeglm_bulked_h5ad.R <pseudobulk.h5ad> 
# <group_column_in_obs> 
# <numerator_group> 
# <denominator_group> 
# <out_dir> 
# [counts_layer]

Rscript $PATH_SCRIPT $PATH_H5AD \
"${REGION}_vs_OtherNeur" \
"${REGION}_${NT}" \
"OtherNeur" \
"${PATH_OUT}/${REGION}_${NT}"

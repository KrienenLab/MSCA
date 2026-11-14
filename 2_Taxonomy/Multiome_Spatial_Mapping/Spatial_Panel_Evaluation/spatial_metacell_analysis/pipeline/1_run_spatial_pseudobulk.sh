#!/usr/bin/env bash

# run from anywhere inside this repository
cd "$(git rev-parse --show-toplevel)/2_Taxonomy/Multiome_Spatial_Mapping/Spatial_Panel_Evaluation/spatial_metacell_analysis"
conda activate /jukebox/scratch/md6347/envs/scvi_daisy

python pipeline/pseudobulk_spatial_by_cluster.py \
--h5ad /jukebox/krienen/aibs_fileshare/hmba-marmoset-wg-802451596237-us-west-2/Subcortex/spatial/20251022_preprint_freeze/Xenium-CJ23.56.004-marmoset-subcortex-cluster-v4_PREPRINT_FREEZE_20251022_subcortex.h5ad \
--out-counts data/spatial_pseudobulk_Cluster_v4.csv \
--out-cluster-stats data/spatial_Cluster_v4_stats.csv

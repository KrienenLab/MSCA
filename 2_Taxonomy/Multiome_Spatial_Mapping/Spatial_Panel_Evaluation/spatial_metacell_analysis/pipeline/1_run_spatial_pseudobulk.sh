#!/usr/bin/env bash

cd /jukebox/krienen/lab_members/mike/260216_Spatial_Panel_Eval/spatial_metacell_analysis
conda activate /jukebox/scratch/md6347/envs/scvi_daisy

python pipeline/pseudobulk_spatial_by_cluster.py \
--h5ad /jukebox/krienen/aibs_fileshare/hmba-marmoset-wg-802451596237-us-west-2/Subcortex/spatial/20251022_preprint_freeze/Xenium-CJ23.56.004-marmoset-subcortex-cluster-v4_PREPRINT_FREEZE_20251022_subcortex.h5ad \
--out-counts /jukebox/krienen/lab_members/mike/260216_Spatial_Panel_Eval/spatial_metacell_analysis/data/spatial_pseudobulk_Cluster_v4.csv \
--out-cluster-stats /jukebox/krienen/lab_members/mike/260216_Spatial_Panel_Eval/spatial_metacell_analysis/data/spatial_Cluster_v4_stats.csv

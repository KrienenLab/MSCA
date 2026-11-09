#!/usr/bin/env bash
# Fig. S5H: genes differentially expressed in the BST-localized microglia clusters vs. the other microglia clusters
# (DESeq2 on the cluster metacells, apeglm-shrunken log2 fold changes; positive = higher in BST clusters).
#
# Run from this folder, in an environment with R (DESeq2, apeglm, reticulate) and a Python with anndata:
#     bash FigS5H_microglia_DEG.sh
#
# Input   data/Microglia_cluster_metacell.h5ad   raw-count metacells of the 14 microglia clusters (obs: Cluster_v4)
# Output  output/DESeq2_apeglm_results.csv, DESeqDataSet.rds, sample_group_mapping.csv

set -euo pipefail

pseudo_bulk="data/Microglia_cluster_metacell.h5ad"
group_column_in_obs="Cluster_v4"
out_dir="output"

# BST-localized clusters vs. all other microglia clusters
numerator_group="cluster_854,cluster_855,cluster_857,cluster_858"
denominator_group="cluster_846,cluster_847,cluster_848,cluster_849,cluster_850,cluster_851,cluster_852,cluster_853,cluster_856,cluster_859"

Rscript RunDESeq2apeglm_bulked_h5ad.R "${pseudo_bulk}" "${group_column_in_obs}" \
    "${numerator_group}" "${denominator_group}" "${out_dir}"

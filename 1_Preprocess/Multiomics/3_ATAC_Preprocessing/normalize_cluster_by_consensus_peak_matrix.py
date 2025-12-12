#!/usr/bin/env python3
from pathlib import Path
import os
import anndata as ad
import numpy as np
from crested.pp import normalize_peaks

pdir           = '/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster'
path_peaks     = Path(pdir, 'data/atac/atac_consensus_peaks.bed')
dir_count_mat  = Path(pdir, "data/atac/atac_count_matrix")
path_in        = Path(dir_count_mat, 'Cluster_v4_by_consensus_peaks.h5ad')
path_out       = Path(dir_count_mat, 'Cluster_v4_by_consensus_peaks_normalized.h5ad')

print("Importing cluster-by-peak metacell count matrix into memory...")
adata = ad.read_h5ad(path_in)

# --------------------------------------------------
# Renormalize peaks
# --------------------------------------------------


print("Normalizing peaks using CREsted constitutive peak method")
#
# top_k_percent (float (default: 0.01)) – The percentage (expressed as a fraction) of top values to consider for Gini score calculation.
# 
# The top_k_percent parameters can be tuned based on potential bias towards cell types. 
# If some weights are overcompensating too much, consider increasing the top_k_percent. 
# 
# -> got all NaN when tried running directly on count matrix
adata.X = adata.X.astype(np.float32)
normalize_peaks(adata, top_k_percent=0.01)

# print("Exporting barplot of the normalization weights")
# p = crested.pl.bar.normalization_weights(
#     adata, title="Normalization Weights by Cluster", x_label_rotation=90, show=False
# ).get_figure()
# p.savefig(Path(plotdir, 'norm_weights_per_cluster.pdf'))
# plt.close(p)

print("Writing out...")
adata.write_h5ad(path_out)

print("Script finished.")

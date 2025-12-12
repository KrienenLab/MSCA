#!/usr/bin/env python3
from pathlib import Path
import anndata as ad

pdir           = '/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster'
path_peak_set  = '/jukebox/krienen/victor/marm_hmba_atac_50k/data/marm_subcortex_filtered_peaks.txt'
dir_count_mat  = Path(pdir, "data/atac/atac_count_matrix")
path_in        = Path(dir_count_mat, 'Cluster_v4_by_consensus_peaks_normalized.h5ad')
path_out       = Path(dir_count_mat, 'Cluster_v4_by_462k_peaks_normalized.h5ad')

def subset_peaks(path_in, path_out, path_peaks):
    # peaks is a txt file with feature names
    print("Importing anndata into memory...")
    adata = ad.read_h5ad(path_in)

    print("Importing peak names from txt file...")
    with open(path_peaks, "r") as f:
        peak_names = [line.strip() for line in f if line.strip()]
    
    assert all([p in adata.var_names for p in peak_names]), "Some peaks not found in adata"

    print(f"Subsetting anndata from {adata.shape[1]} to {len(peak_names)} peaks")
    adata = adata[:, peak_names].copy()

    print("Writing out...")
    adata.write_h5ad(path_out)
    return None

# Raw metacells
print("Subsetting 'Cluster_v4_by_consensus_peaks.h5ad'")
subset_peaks(
    Path(dir_count_mat, 'Cluster_v4_by_consensus_peaks.h5ad'),
    Path(dir_count_mat, 'Cluster_v4_by_462k_peaks.h5ad'),
    path_peak_set
)
print("Done.")

# Normalized metacells
print("Subsetting 'Cluster_v4_by_consensus_peaks_normalized.h5ad'")
subset_peaks(
    Path(dir_count_mat, 'Cluster_v4_by_consensus_peaks_normalized.h5ad'),
    Path(dir_count_mat, 'Cluster_v4_by_462k_peaks_normalized.h5ad'),
    path_peak_set
)
print("Done.")

# Raw cells ----> return w/ more memory...
# print("Subsetting 'cell_by_consensus_peaks.h5ad'")
# subset_peaks(
#     Path(dir_count_mat, 'cell_by_consensus_peaks.h5ad'),
#     Path(dir_count_mat, 'cell_by_462k_peaks.h5ad'),
#     path_peak_set
# )
# print("Done.")

print("Script finished.")

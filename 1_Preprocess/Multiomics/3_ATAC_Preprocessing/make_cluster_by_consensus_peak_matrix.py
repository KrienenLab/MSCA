#!/usr/bin/env python
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp
from pathlib import Path

pdir          = Path("/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster")
dir_count_mat = Path(pdir, "data/atac/atac_count_matrix")
path_in       = Path(dir_count_mat, "cell_by_consensus_peaks.h5ad")
path_out      = Path(dir_count_mat, "Cluster_v4_by_consensus_peaks.h5ad")
cluster_col   = "Cluster_v4"

print("Importing cell-by-peak count matrix into memory...")
adata = ad.read_h5ad(path_in)

print("Getting one-hot encoding for clusters...")
codes      = adata.obs[cluster_col].astype("category").cat.codes.to_numpy()
cats       = adata.obs[cluster_col].astype("category").cat.categories
n_cells    = adata.shape[0]
n_clusters = cats.size

# sparse one-hot matrix (n_clusters × n_cells)
# (and snapatac2 also uses int32 by default)
indicator = sp.csr_matrix(
    (np.ones(n_cells, dtype=np.int32),  # fill 1s for
     (codes, np.arange(n_cells))),      # (row, col) pairs
    shape=(n_clusters, n_cells),
    dtype=np.int32,
)

print("Matrix multiplication to generate sparse cluster-by-peak matrix...")
# sparse matrix–matrix multiplication automatically sums over cells that belong to the same cluster
cluster_peak = indicator @ adata.X

print("Making & saving AnnData/h5ad...")
cluster_adata = ad.AnnData(
    X   = cluster_peak.astype(np.int32),
    obs = pd.DataFrame(index=cats),
    var = adata.var.copy()
)
cluster_adata.write_h5ad(path_out)

adata.file.close() 
print("Script finished.")


######################################
# Code tests
######################################

if False:
    import numpy as np
    import scipy.sparse as sp

    # Parameters for the toy dataset
    n_cells = 100
    n_peaks = 50
    n_clusters = 7
    rng = np.random.default_rng(0)

    # Assign each cell to a random cluster
    codes = rng.integers(0, n_clusters, n_cells)

    # Create a sparse cell-by-peak matrix with integer counts (1–9)
    X = sp.random(
        n_cells,
        n_peaks,
        density=0.1,
        format="csr",
        data_rvs=lambda s: rng.integers(1, 10, size=s, dtype=np.int32),
    )
    X.data = X.data.astype(np.int32)  # ensure int32 dtype

    # Build indicator matrix (n_clusters × n_cells)
    indicator = sp.csr_matrix(
        (
            np.ones(n_cells, dtype=np.int32),  # data: all 1s
            (codes, np.arange(n_cells)),       # (row, col) coordinates
        ),
        shape=(n_clusters, n_cells),
        dtype=np.int32,
    )

    # Collapse cell counts to cluster counts
    cluster_peak = indicator @ X

    # Section 3: dtype sanity & random spot‑checks
    print("input X dtype      :", X.dtype)
    print("cluster_peak dtype :", cluster_peak.dtype)

    for _ in range(10):
        cl = rng.integers(0, n_clusters)
        pk = rng.integers(0, n_peaks)
        rows = np.where(codes == cl)[0]
        direct_sum = X[rows, pk].sum()
        assert cluster_peak[cl, pk] == direct_sum, f"Mismatch at cluster {cl}, peak {pk}"

    print("✓ All random spot‑checks passed without errors.")

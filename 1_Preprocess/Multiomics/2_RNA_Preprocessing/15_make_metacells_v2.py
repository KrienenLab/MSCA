#!/usr/bin/env python
import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp
from pathlib import Path

pdir      = Path("/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster")
path_in   = Path(pdir, "data/rna/rna_raw_v4_freeze.h5ad")
path_anno = Path(pdir, "metadata/preprint_freeze/subcortex_anno_table_preprint.csv")
dir_out   = Path(pdir, "data/rna/metacell_v2/subcortex")
dir_out.mkdir(parents=True, exist_ok=True)

adata = ad.read_h5ad(path_in)
anno  = pd.read_csv(path_anno).set_index('marm_cluster_id')

# keep nothing but cluster IDs from h5ad
adata.obs = adata.obs.loc[:, ['Cluster_v4']]

# keep only the annotation columns being used for metacells (not including clusters)
anno = anno.loc[:, ['Subcortex_Group_v4']]
assert all(adata.obs.Cluster_v4.isin(anno.index)), "Not all clusters in annotation table"
adata.obs = adata.obs.join(anno, on='Cluster_v4') 

def make_metacell(adata, cluster_col, normalize = True):
    print(f"Getting one-hot encoding for `{cluster_col}`...")
    codes      = adata.obs[cluster_col].astype("category").cat.codes.to_numpy()
    cats       = adata.obs[cluster_col].astype("category").cat.categories
    n_cells    = adata.shape[0]
    n_clusters = cats.size

    # sparse one-hot matrix (n_clusters x n_cells)
    indicator = sp.csr_matrix(
        (np.ones(n_cells, dtype=np.int32),  # fill 1s for (row, col) pairs
         (codes, np.arange(n_cells))),      
        shape=(n_clusters, n_cells),
        dtype=np.int32,
    )

    print("Generating sparse metacell count matrix...")
    bulked = indicator @ adata.X

    print("Creating new AnnData...")
    cluster_adata = ad.AnnData(
        X   = bulked.astype(np.int32),
        obs = pd.DataFrame(index=cats),
        var = adata.var.copy()
    )

    if normalize:
        print("Creating layer 'counts' for raw and log1p(CPM) normalizing metacells in X...")
        import scanpy as sc
        cluster_adata.layers['counts'] = cluster_adata.X
        sc.pp.normalize_total(cluster_adata, target_sum=1e6)
        sc.pp.log1p(cluster_adata)
    
    return(cluster_adata)


make_metacell(adata, "Cluster_v4").write_h5ad(Path(dir_out, "Cluster_v4_metacell.h5ad"))
make_metacell(adata, "Subcortex_Group_v4").write_h5ad(Path(dir_out, "Subcortex_Group_v4_metacell.h5ad"))


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

    # Build indicator matrix (n_clusters x n_cells)
    indicator = sp.csr_matrix(
        (
            np.ones(n_cells, dtype=np.int32),  # data: all 1s
            (codes, np.arange(n_cells)),       # (row, col) coordinates
        ),
        shape=(n_clusters, n_cells),
        dtype=np.int32,
    )

    # Collapse cell counts to cluster counts
    bulked = indicator @ X

    # Section 3: dtype sanity & random spot‑checks
    print("input X dtype      :", X.dtype)
    print("bulked dtype :", bulked.dtype)

    for _ in range(10):
        cl = rng.integers(0, n_clusters)
        pk = rng.integers(0, n_peaks)
        rows = np.where(codes == cl)[0]
        direct_sum = X[rows, pk].sum()
        assert bulked[cl, pk] == direct_sum, f"Mismatch at cluster {cl}, peak {pk}"

    print("All tests passed without errors.")

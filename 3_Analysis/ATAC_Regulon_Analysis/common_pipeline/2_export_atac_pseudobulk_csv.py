#!/usr/bin/env python
# Exporting csv file of the matrix to facilitate some analysis in R.
import anndata as ad
import pandas as pd
from pathlib import Path

pdir          = Path("/jukebox/krienen/lab_members/mike/251022_HMBA_Spatial_ATAC")
pdir_dat      = Path("/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster")
dir_count_mat = Path(pdir_dat, "data/atac/atac_count_matrix")
path_cl_by_pk = Path(dir_count_mat, "Cluster_v4_by_462k_peaks.h5ad")
outdir        = Path(pdir, "common_data")


# --------------------------------------------------
# Imports
# --------------------------------------------------

# Import cell-by-peak matrix (raw)
cl_by_pk = ad.read_h5ad(path_cl_by_pk)


# --------------------------------------------------
# Export csv files for convenience in R
# --------------------------------------------------

from scipy import sparse

# transpose matrix (faster r/w); and densify if sparse (which it is)
X = cl_by_pk.X.T
if sparse.issparse(X):
   import numpy as np
   X = X.astype(np.int64).toarray()  

pd.DataFrame(
   X, columns=cl_by_pk.obs_names, index=cl_by_pk.var_names
).to_csv(Path(outdir, "Cluster_v4_by_462k_peaks_matrix.csv"))

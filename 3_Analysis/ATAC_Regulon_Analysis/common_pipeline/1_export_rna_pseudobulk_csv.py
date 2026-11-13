#!/usr/bin/env python
# Exporting csv file of the matrix to facilitate some analysis in R.
import anndata as ad
import pandas as pd
from pathlib import Path

pdir          = Path("/jukebox/krienen/lab_members/mike/251022_HMBA_Spatial_ATAC")
pdir_dat      = Path("/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster")
path_cl_by_gn = Path(pdir_dat, "data/rna/metacell_v2/subcortex/Cluster_v4_metacell.h5ad")
outdir        = Path(pdir, "common_data")


# --------------------------------------------------
# Imports
# --------------------------------------------------

# Import cell-by-peak matrix (normalized in X; so take the raw)
cl_by_gn = ad.read_h5ad(path_cl_by_gn)
cl_by_gn.X = cl_by_gn.layers['counts']
del cl_by_gn.layers['counts'], cl_by_gn.uns['log1p']

# --------------------------------------------------
# Export csv files for convenience in R
# --------------------------------------------------

from scipy import sparse

# transpose matrix (faster r/w); and densify if sparse (which it is)
X = cl_by_gn.X.T
if sparse.issparse(X):
   import numpy as np
   X = X.astype(np.int64).toarray()  

pd.DataFrame(
   X, columns=cl_by_gn.obs_names, index=cl_by_gn.var_names
).to_csv(Path(outdir, "Cluster_v4_by_genes_matrix.csv"))

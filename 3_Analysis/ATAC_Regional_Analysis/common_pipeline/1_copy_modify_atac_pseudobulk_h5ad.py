#!/usr/bin/env python
#
# For my convenience, will copy the h5ad pseudobulk ATAC file into
# my analysis directory. Here, the only modification will be joining
# in the complete annotation table (as of 20251008) and the 'Dominant_Region' 
# annotations, but more modifications will be made in other analysis scripts.
#
# Using the Coarse (group) & Fine (cluster) region annotations, I will generate
# a single set of "Dominant_Region" annotations that apply at the cluster-level. 
# These will use those cluster-level refinements where available, and default
# to the group-level annotation everywhere else.
#
# Will also be exporting csv files of the matrix and the obs to
# facilitate some analysis in R.
import anndata as ad
import pandas as pd
from pathlib import Path

repo          = Path(__file__).resolve().parents[3]  # root of this repository
pdir          = Path(repo, "3_Analysis/ATAC_Regional_Analysis")
pdir_dat      = Path("/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster")
dir_count_mat = Path(pdir_dat, "data/atac/atac_count_matrix")
path_cl_by_pk = Path(dir_count_mat, "Cluster_v4_by_462k_peaks.h5ad")
path_anno_tbl = Path(pdir_dat, "metadata/subcortex_anno_table_archive/subcortex_anno_table_20251008.csv")
path_reg_anno = Path(pdir_dat, "metadata/subcortex_anno_table_archive/subcortex_anno_table_cluster_region_20251008.csv")
outdir        = Path(pdir, "common_data")


# --------------------------------------------------
# Imports
# --------------------------------------------------

# Import cell-by-peak matrix (raw)
cl_by_pk = ad.read_h5ad(path_cl_by_pk)

# Import main annotation table
anno_tbl = pd.read_csv(path_anno_tbl, index_col='marm_cluster_id')

# Import current region annotations [note junk clusters not included in that file]
reg_anno = pd.read_csv(path_reg_anno, index_col="marm_cluster_id")
assert all(reg_anno.index.isin(anno_tbl.index))
assert all(reg_anno.Subcortex_Group_v4.isin(anno_tbl.Subcortex_Group_v4))


# --------------------------------------------------
# Make a merged 'Dominant_Region' column
# --------------------------------------------------

# >>> reg_anno.Dominant_Region_Fine.drop_duplicates().to_list()
# [nan, 'STR', 'Multi', 'AMY', 'CTX', 'HiF']
# >>> reg_anno.Dominant_Region_Coarse.drop_duplicates().to_list()
# [nan, 'THM', 'MB', 'STR', 'AMY', 'HTH', 'HB', 'Multi', 'CTX', 'HiF', 'GP', 'BF', 'OlfC', 'CLEn', 'EXA', 'Pin']

assert all(reg_anno.Dominant_Region_Fine.isin(reg_anno.Dominant_Region_Coarse))

reg_anno['Dominant_Region'] = reg_anno['Dominant_Region_Coarse']
# >>> sum(reg_anno.Dominant_Region == reg_anno.Dominant_Region_Coarse)
# 1633
idx = reg_anno.index[~reg_anno['Dominant_Region_Fine'].isna()]
# >>> len(idx)
# 277
reg_anno.loc[idx, 'Dominant_Region'] = reg_anno.loc[idx, 'Dominant_Region_Fine']
# >>> sum(reg_anno.Dominant_Region == reg_anno.Dominant_Region_Coarse)
# 1533
# >>> sum(reg_anno.Dominant_Region == reg_anno.Dominant_Region_Fine)
# 277


# --------------------------------------------------
# Join Dominant_Region annotations into larger annotation table
# --------------------------------------------------

reg_anno = reg_anno.loc[:, ['Dominant_Region_Fine', 'Dominant_Region_Coarse', 'Dominant_Region']]
assert not any(reg_anno.columns.isin(anno_tbl.columns))
anno_tbl = anno_tbl.join(reg_anno)


# --------------------------------------------------
# Join into anndata and export
# --------------------------------------------------

assert all(cl_by_pk.obs_names.isin(anno_tbl.index))
cl_by_pk.obs = cl_by_pk.obs.join(anno_tbl)

## drop problem columns (although I should probably just select the columns I want...)
## -> actually, so annoying to have all those columns, will make a selection instead
# cl_by_pk.obs.drop(['Need Help', 'Preprint_include'], axis=1, inplace=True)

cols_keep = [
   'Subcortex_Class_v4', 'Subcortex_Group_v4', 'NT', 'Subcortex_Subclass_hcpc_v4', 'N_cells_in_spatial_mapping', 
   'CDM_Cluster_v4_bootstrapping_probability_median', 'CDM_Cluster_v4_aggregate_probability_median', 'CDM_Cluster_v4_avg_correlation_median', 
   'nsforest_markers', 'Dominant_Region_Fine', 'Dominant_Region_Coarse', 'Dominant_Region'
]
cl_by_pk.obs = cl_by_pk.obs.loc[:, cols_keep]

# Export
cl_by_pk.write_h5ad(Path(outdir, "Cluster_v4_by_462k_peaks_modcopy.h5ad"))


# --------------------------------------------------
# Export csv files for convenience in R
# --------------------------------------------------

cl_by_pk.obs.to_csv(Path(outdir, "Cluster_v4_by_462k_peaks_obs.csv"))
cl_by_pk.var.to_csv(Path(outdir, "Cluster_v4_by_462k_peaks_var.csv"))

from scipy import sparse

# transpose matrix (faster r/w); and densify if sparse (which it is)
X = cl_by_pk.X.T
if sparse.issparse(X):
   import numpy as np
   X = X.astype(np.int64).toarray()  

pd.DataFrame(
   X, columns=cl_by_pk.obs_names, index=cl_by_pk.var_names
).to_csv(Path(outdir, "Cluster_v4_by_462k_peaks_matrix.csv"))

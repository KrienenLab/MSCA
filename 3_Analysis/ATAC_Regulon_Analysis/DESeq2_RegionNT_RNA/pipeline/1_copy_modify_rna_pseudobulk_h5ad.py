#!/usr/bin/env python
#
# Copying and modify RNA pseudobulk h5ad file to have 
# additional columns to code for various DESeq2
# comparisons. Each column will be f'{region}_vs_OtherNeur',
# and will, for each region, code its f'{Region}_{NT}',
# code all other neurons as 'OtherNeur', and code non-neurons
# and junk differently (to not be involved in the comparisons).
#
# This script will also separate MSNs from GABA (so treat like 
# another neurotransmitter).
# 
# The text file output is for a bash script to run each DESeq2
# test in a slurm job
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path

pdir            = Path("/jukebox/krienen/lab_members/mike/251022_HMBA_Spatial_ATAC/DESeq2_RegionNT_RNA")
pdir_dat        = Path("/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster")
path_cl_by_gn   = Path(pdir_dat, "data/rna/metacell_v2/subcortex/Cluster_v4_metacell.h5ad")
path_anno_tbl   = Path(pdir_dat, "metadata/preprint_freeze/subcortex_anno_table_preprint.csv")
path_h5ad_out   = Path(pdir, "data/Cluster_v4_by_genes_modcopy.h5ad")
out_list_reg_nt = Path(pdir, "data/all_region_nt.txt")


# --------------------------------------------------
# Imports
# --------------------------------------------------

# Import cell-by-peak matrix (normalized in X; so take the raw)
adata = ad.read_h5ad(path_cl_by_gn)
adata.X = adata.layers['counts']
del adata.layers['counts'], adata.uns['log1p']

# Import main annotation table
anno_tbl = pd.read_csv(path_anno_tbl, index_col='marm_cluster_id')


# --------------------------------------------------
# Join obs
# --------------------------------------------------

assert all(adata.obs_names.isin(anno_tbl.index)), "Not all clusters in anno table"
adata.obs = adata.obs.join(anno_tbl)


# --------------------------------------------------
# Subset anndata to preprint_analysis_include
# --------------------------------------------------

# >>> adata.obs.loc[:, ['Subcortex_Class_v4', 'preprint_analysis_include']].value_counts()
# Subcortex_Class_v4  preprint_analysis_include
# BG                  True                         403
# TEL Glut            True                         354
# THM Glut            True                         178
# MB-HB Glut          True                         115
# Nonneuron           True                         113
# MB-HB GABA          True                         112
# CNU-HTH GABA        True                          87
# HTH Glut            True                          76
# Junk                False                         74 [check]
# HTH GABA            True                          57
# CNU-HTH Glut        True                          47
# CNU-HTH GABA        False                         43 [?]
# SEP GABA            True                          38
# THM GABA            True                          38
# HTH GABA            False                         22 [?]
# CTX-CGE GABA        True                          14
# CNU-HTH Glut        False                         11 [?]
# SEP GABA            False                          8 [?]
# CTX-MGE GABA        True                           5
# MB-HTH-HB Sero      True                           4
# TEL Glut            False                          4
# THM Glut            False                          4
# MB-HB Glut          False                          2
# HTH Hist            True                           1
# Pineal Gland        False                          1
# Name: count, dtype: int64

###
### note that's 22% of clusters w/ dominant region HTH being dropped (45/203); total of 76 clusters w/ an HTH
###
# >>> adata.obs.loc[adata.obs.Cluster_Dominant_Region == "HTH", ['preprint_analysis_include']].value_counts()
# preprint_analysis_include
# True                         158
# False                         45


# >>> adata.shape
# (1811, 35787)
adata = adata[adata.obs.preprint_analysis_include, ].copy()
# >>> adata.shape
# (1642, 35787)


# --------------------------------------------------
# Create obs column for each region's tests
# --------------------------------------------------

### the MSNs:
# adata.obs.Subcortex_Group_v4[
#     adata.obs.Subcortex_Group_v4.str.contains("MSN")
# ].drop_duplicates().to_list()
# ['STRd D2 Striosome MSN', 'STRd D2 Matrix MSN', 'STRd D1 Matrix MSN', 'STRv D2 MSN', 'STRv D1 MSN', 'STRd D1 Striosome MSN', 'STRd D1D2 Hybrid MSN', 'AMY-STR ITC-NUDAP MSN', 'STRd D2 StrioMat Hybrid MSN']

# For each region, will compare each NT to all other neurons outside the region [although dropped 'Pineal Gland']
nn_classes = ['Nonneuron', 'Pineal Gland']

all_regions = sorted(
    [x for x in set(adata.obs['Cluster_Dominant_Region']) if pd.notna(x)]
)

# >>> all_regions
# ['AMY', 'BF', 'CLEn', 'CTX', 'EXA', 'GP', 'HB', 'HTH', 'HiF', 'MB', 'Multi', 'OlfC', 'STR', 'THM']

for region in all_regions:
    k = f'{region}_vs_OtherNeur'
    msn_clusts = adata.obs.index[
        (adata.obs['Subcortex_Group_v4'].str.contains("MSN")) & 
        (adata.obs['Cluster_Dominant_Region'] == region)
    ].drop_duplicates().to_list()
    adata.obs[k] = adata.obs['Cluster_Dominant_Region'].astype(str) + "_" + adata.obs['Cluster_NT'].astype(str)
    adata.obs.loc[adata.obs.index.isin(msn_clusts), k] = region + "_MSN"
    adata.obs.loc[adata.obs['Cluster_Dominant_Region'] != region, k] = 'OtherNeur'
    adata.obs.loc[adata.obs['Subcortex_Class_v4'].isin(nn_classes), k] = 'NN'


# --------------------------------------------------
# Write out modified h5ad file [for use with DESeq2]
# --------------------------------------------------

adata.write_h5ad(path_h5ad_out)


# --------------------------------------------------
# Also write out the set of all region_nt pairs
# --------------------------------------------------

## Each Region+NT where there's >1 cluster
valid_series = []
for region in all_regions:
    if region == "Multi":
        continue
    s = adata.obs[f'{region}_vs_OtherNeur']
    s = s[~s.isin(['OtherNeur', 'NN', 'Junk'])]  # keep only Region_NT labels
    valid_series.append(s.astype(str))

all_pairs = pd.concat(valid_series, ignore_index=True)
counts = all_pairs.value_counts()
pairs = sorted(counts[counts > 1].index.tolist())

with open(out_list_reg_nt, "w") as f:
    f.write("\n".join(map(str, pairs)))
    f.write("\n")

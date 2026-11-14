#!/usr/bin/env python
# Will modify the same h5ad already copied/modified to have
# the 'Dominant_Region' annotation.
#
# This script will add columns to code for various DESeq2
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

repo            = Path(__file__).resolve().parents[4]  # root of this repository
pdir            = Path(repo, "3_Analysis/ATAC_Regional_Analysis")
path_h5ad       = Path(pdir, "common_data/Cluster_v4_by_462k_peaks_modcopy.h5ad")
out_list_reg_nt = Path(pdir, "DESeq2_AllRegionNT/data/all_region_nt.txt")


# --------------------------------------------------
# Imports
# --------------------------------------------------

adata = ad.read_h5ad(path_h5ad, backed="r+")


# --------------------------------------------------
# Create obs column for each region's tests
# --------------------------------------------------

## the MSNs:
# >>> adata.obs.Subcortex_Group_v4[
# ...     adata.obs.Subcortex_Group_v4.str.contains("MSN")
# ... ].drop_duplicates().to_list()
# ['STRd D2 Striosome MSN', 'STRd D2 Matrix MSN', 'STRd D1 Matrix MSN', 'STRv D2 MSN', 'STRv D1 MSN', 'STRd D1 Striosome MSN', 'STRd D1D2 Hybrid MSN', 'AMY-STR NUDAP-ITC MSN', 'STRd D2 StrioMat Hybrid MSN']

# For each region, will compare each NT to all other neurons outside the region
nn_classes = ['Nonneuron', 'Pineal Gland']

all_regions = sorted(
    [x for x in set(adata.obs['Dominant_Region']) if pd.notna(x)]
)

# >>> all_regions
# ['AMY', 'BF', 'CLEn', 'CTX', 'EXA', 'GP', 'HB', 'HTH', 'HiF', 'MB', 'Multi', 'OlfC', 'Pin', 'STR', 'THM']

for region in all_regions:
    k = f'{region}_vs_OtherNeur'
    msn_clusts = adata.obs.index[
        (adata.obs.Subcortex_Group_v4.str.contains("MSN")) & 
        (adata.obs.Dominant_Region == region)
    ].drop_duplicates().to_list()
    adata.obs[k] = adata.obs['Dominant_Region'].astype(str) + "_" + adata.obs['NT'].astype(str)
    adata.obs.loc[adata.obs.index.isin(msn_clusts), k] = region + "_MSN"
    adata.obs.loc[adata.obs['Dominant_Region'] != region, k] = 'OtherNeur'
    adata.obs.loc[adata.obs['Subcortex_Class_v4'].isin(nn_classes), k] = 'NN'
    adata.obs.loc[adata.obs['Subcortex_Class_v4'] == "Junk", k] = 'Junk'


# --------------------------------------------------
# Write out modified h5ad file [for use with DESeq2]
# --------------------------------------------------

adata.write_h5ad(path_h5ad)


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

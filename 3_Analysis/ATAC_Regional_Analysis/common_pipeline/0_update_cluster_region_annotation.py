#!/usr/bin/env python
#
# Regional annotations ("Dominant_Region") were first made at the Group level
# (Dominant_Region_Coarse) and later refined for individual clusters
# (Dominant_Region_Fine). For example, a number of hippocampal clusters belong
# to Groups that otherwise localize to cortex.
#
# The cluster-level region file was made before the 2025-10-08 annotation
# freeze. This script updates its Class, Group and neurotransmitter columns to
# match that freeze, and overwrites the file. The region columns are untouched.
import pandas as pd
from pathlib import Path

pdir_dat      = Path("/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster")
path_anno_tbl = Path(pdir_dat, "metadata/subcortex_anno_table_archive/subcortex_anno_table_20251008.csv")
path_reg_anno = Path(pdir_dat, "metadata/subcortex_anno_table_archive/subcortex_anno_table_cluster_region_20251008.csv")


# --------------------------------------------------
# Imports
# --------------------------------------------------

# Import main annotation table (to verify groups still match)
anno_tbl = pd.read_csv(path_anno_tbl)
anno_tbl.set_index('marm_cluster_id', inplace=True)

# Import target file for having region annotations
reg_anno = pd.read_csv(path_reg_anno, index_col="marm_cluster_id")
# assert all(reg_anno.Subcortex_Group_v4.isin(anno_tbl.Subcortex_Group_v4))

# also this appeared (numbered rows, will get rid of those)
if "Unnamed: 0" in reg_anno.columns:
    reg_anno.drop("Unnamed: 0", axis=1, inplace=True)

# --------------------------------------------------
# Update the cluster_region file to match anno table
# --------------------------------------------------

for k in ["Subcortex_Class_v4", "Subcortex_Group_v4", "NT"]:
    reg_anno[k] = anno_tbl.loc[reg_anno.index, k]


# --------------------------------------------------
# Overwrite file
# --------------------------------------------------

reg_anno.to_csv(path_reg_anno)

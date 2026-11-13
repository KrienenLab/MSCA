#!/usr/bin/env python
#
# We had a set of "Dominant_Region" annotations made at the Group level.
# Now, we have an additional set of refinements for some clusters that 
# result in some regions have much better representation. (For example,
# a number of hippocampal clusters are in groups that otherwise localize
# generally to cortex).
# 
# Daisy had previously made a cluster-level file that contained the current
# group names (as of 251006) and had the Dominant_Region_Coarse (group-level)
# included, but Fenna made some updates on night of 251006, so now (251007)
# I will update that file to contain the updated Dominant_Region_Fine.
#
# Then that version had a mistake in the input file and was missing some
# region assignments (I don't know...), so I have to redo.
#
# Ok, 2025-10-08, we made a new freeze, I'll align the regional annotations 
# with that, and keep everything in the archive folder. This way at least 
# my annotations will match something (and not be in the middle
# of ongoing updates).
#
# SO NOW all this script will do is update the groups and classes in the
# regional annotation file, and leave everything else to be done later.
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


# ######## SCRAP (FOR THE H5AD MOD-COPY SCRIPT) ###########


# --------------------------------------------------
# Update groups for 251007 update [might as well!]
# --------------------------------------------------

# From the 251007 update, below are the changes I gather in the group names.
# Not sure, will see about this comment...:
# ---
#   Group NDB-SI LHX6-LHX8-GBX1 GABA (aka NDB-SI LHX8-STK32A GABA in bg annotation) 
#   is no longer included in the BG preprint. This group is now merged with 
#   NDB-SI-MPA-LPA LHX6-LHX8-GBX1 GABA, and also changed from the BG class to CNU-HTH GABA class. 
# ---

# group_update = {
#     "Oligo RBFOX1"         : "Oligo PLEKHG1",
#     "STR PVALB-RSPO2 GABA" : "STR LYPD6-RSPO2 GABA",
#     "STR SST GABA"         : "STR LYPD6-RSPO2 GABA",
#     "OB Dopa-Gaba"         : "OB Dopa-GABA",
#     "OB-in FRMD7 GABA"     : "OB FRMD7 GABA"
# }

# upd_mask = reg_anno.Subcortex_Group_v4.isin(group_update)
# reg_anno.loc[upd_mask, "Subcortex_Group_v4"] = [
#     group_update[x] for x in reg_anno.Subcortex_Group_v4[upd_mask]
# ]
## assert all(reg_anno.Subcortex_Group_v4.isin(anno_tbl.Subcortex_Group_v4))

## still assertion error
# >>> reg_anno.Subcortex_Group_v4[~reg_anno.Subcortex_Group_v4.isin(anno_tbl.Subcortex_Group_v4)].drop_duplicates().to_list()
# ['Oligo PLEKHG1', 'NDB-SI LHX6-LHX8-GBX1 GABA', 'STR LYPD6-RSPO2 GABA', 'OB Dopa-GABA', 'OB FRMD7 GABA', 'L5 PrS Glut']

## ok i see that went a little wrong

###
### so... at this point, I actually don't need to figure out the group updates, because I have cluster-indexed region labels.
### i just need to verify after-the-fact the Dominant_Region_Coarse labeling logic is still intact (so check for groups that
### changed if the Dominant_Region_Coarse would become ambiguous)
###

# ok so start over, re-import the original tables

# >>> anno_tbl.shape
# (1811, 144)
# >>> reg_anno.shape
# (1746, 7)

tst = pd.DataFrame(
    {
        "grp_old" : reg_anno.Subcortex_Group_v4, 
        "grp_new" : anno_tbl.loc[reg_anno.index, 'Subcortex_Group_v4']
    },
    index = reg_anno.index
)
# >>> tst.loc[tst.grp_old != tst.grp_new, :]
#                                     grp_old                             grp_new
# marm_cluster_id
# cluster_911      NDB-SI LHX6-LHX8-GBX1 GABA  NDB-SI-MPA-LPA LHX6-LHX8-GBX1 GABA
# cluster_912      NDB-SI LHX6-LHX8-GBX1 GABA  NDB-SI-MPA-LPA LHX6-LHX8-GBX1 GABA
# cluster_914      NDB-SI LHX6-LHX8-GBX1 GABA  NDB-SI-MPA-LPA LHX6-LHX8-GBX1 GABA
# cluster_916      NDB-SI LHX6-LHX8-GBX1 GABA  NDB-SI-MPA-LPA LHX6-LHX8-GBX1 GABA
# cluster_917      NDB-SI LHX6-LHX8-GBX1 GABA  NDB-SI-MPA-LPA LHX6-LHX8-GBX1 GABA
# cluster_1430                    L5 PrS Glut                                 NaN



# --------------------------------------------------
# Check entries
# --------------------------------------------------

# [!] Should match 'Dominant_Region_Coarse'

# >>> list(reg_upd.for_mike.drop_duplicates())
# [nan, 'MSN', 'Junk', 'Multi', 'STR', 'AMY', 'CTX', 'HPF']

# >>> list(reg_anno.Dominant_Region_Fine.drop_duplicates())
# [nan, 'MSN', 'STR', 'AMY', 'CTX', 'Multi', 'HPF']

# >>> list(reg_anno.Dominant_Region_Coarse.drop_duplicates())
# [nan, 'THM', 'MB', 'STR', 'AMY', 'HTH', 'HB', 'Multi', 'CTX', 'HiF', 'GP', 'BF', 'OlfC', 'EXA', 'Pin']


mapping = {
    'MSN'   : 'STR',
    'Junk'  : 'Multi',
    'Multi' : 'Multi',
    'STR'   : 'STR',
    'AMY'   : 'AMY',
    'CTX'   : 'CTX',
    'HPF'   : 'HiF',

}
reg_upd['for_mike'] = reg_upd['for_mike'].map(mapping)


# --------------------------------------------------
# Overwrite Dominant_Region_Fine
# --------------------------------------------------

reg_anno['Dominant_Region_Fine'] = reg_upd['for_mike']


# --------------------------------------------------
# Make a unified Dominant_Region column
# --------------------------------------------------

reg_anno['Dominant_Region'] = reg_anno['Dominant_Region_Coarse']
# >>> sum(reg_anno.Dominant_Region == reg_anno.Dominant_Region_Coarse)
# 1453
idx = reg_anno.index[~reg_anno['Dominant_Region_Fine'].isna()]
# >>> len(idx)
# 277
reg_anno.loc[idx, 'Dominant_Region'] = reg_anno.loc[idx, 'Dominant_Region_Fine']
# >>> sum(reg_anno.Dominant_Region == reg_anno.Dominant_Region_Coarse)
# 1354
# >>> sum(reg_anno.Dominant_Region == reg_anno.Dominant_Region_Fine)
# 277


# --------------------------------------------------
# Overwrite target file
# --------------------------------------------------

reg_anno.to_csv(path_reg_anno)

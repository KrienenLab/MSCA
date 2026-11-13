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
# Also, using these two sets of annotations, I will generate a single set of 
# "Dominant_Region" annotations that apply at the cluster-level. These
# will use those cluster-level refinements where available, and default
# to the group-level annotation everywhere else.
# 
import pandas as pd
from pathlib import Path

pdir_dat      = Path("/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster")
path_anno_tbl = Path(pdir_dat, "metadata/subcortex_anno_table.csv")
path_reg_anno = Path(pdir_dat, "metadata/subcortex_anno_table_cluster_region_1007.csv")
path_reg_upd  = Path("/jukebox/krienen/lab_members/mike/251007_HMBA_Spatial_ATAC/common_data/251007_Marmoset_subcortex_annot_table.csv")


# --------------------------------------------------
# Imports
# --------------------------------------------------

# Import main annotation table
anno_tbl = pd.read_csv(path_anno_tbl)
anno_tbl.set_index('marm_cluster_id', inplace=True)

# Import target file for having region annotations
reg_anno = pd.read_csv(path_reg_anno, index_col="marm_cluster_id")
assert all(reg_anno.Subcortex_Group_v4.isin(anno_tbl.Subcortex_Group_v4))

# Import updated anno table that contains updated Dominant_Region_Fine (in column "for_mike")
reg_upd = pd.read_csv(path_reg_upd).set_index('marm_cluster_id')

# >>> anno_tbl.shape
# (1811, 144)

# >>> reg_upd.shape
# (1811, 144)

# >>> reg_anno.shape
# (1746, 6)

# >>> anno_tbl.loc[~anno_tbl.index.isin(reg_anno.index), "Subcortex_Class_v4"].drop_duplicates().to_list()
# ['Junk']

### I actually don't need anno_tbl for anything

# Subset to non-junk (to match the dominant region table...)
reg_upd = reg_upd.loc[reg_anno.index, ]
assert all(reg_anno.index == reg_upd.index)


# --------------------------------------------------
# Update entries
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
# Overwrite target file
# --------------------------------------------------

reg_anno.to_csv(path_reg_anno)

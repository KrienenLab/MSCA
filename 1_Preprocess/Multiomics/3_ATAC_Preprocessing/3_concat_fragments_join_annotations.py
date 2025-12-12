#!/usr/bin/env python3
import snapatac2 as snap
from pathlib import Path
import re
import pandas as pd
import polars as pl

pdir = '/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster'
indir = Path(pdir, 'data/atac/atac_h5ads/anndatas')
path_anno = Path(pdir, 'metadata/subcortex_anno_table.csv')
outfile = Path(pdir, 'data/atac/atac_h5ads/_dataset.h5ads')

# --------------------------------------------------
# Import annotation look-up table
# --------------------------------------------------

# 
# Have 'Cluster_v4' in adata for joining.
# In the annotation file, it's 'marm_cluster_id'. Will rename
# to keep as 'Cluster_v4' and take the columns highlighted by Daisy
# to join into the data.
# 

anno = pd.read_csv(path_anno).rename(columns={'marm_cluster_id': 'Cluster_v4'})
anno.set_index('Cluster_v4', inplace=True)
cols_add = ['Subcortex_Class_v4', 'Subcortex_Group_v4']


# --------------------------------------------------
# Concatenate ATAC h5ad files & their metadata
# --------------------------------------------------

#
# This is tedious, complete lack of methods in snapatac's API for this.
# I need to extract all the metadata from obs of each input file, concatenate
# and save that for each column, and then repopulate obs in the new concatenated 
# _dataset.h5ad file
# 

h5ad_files = {
    re.sub(".h5ad$", "", p.name): snap.read(p, backed='r+') 
    for p in indir.iterdir() if p.suffix == '.h5ad'
}

# Since I also CANNOT figure out a way to access the column/element names in .obs from each file,
# so here are all of the columns that were added in the initialization script,
# plus the computed metrics
#
# All I could figure out is: 
#   ad_mem = h5ad_files["P0023_1"].to_memory() 
# which works but whatever
cols_added = [
    'barcoded_cell_sample_name', 'tissue_name', 'dissociated_cell_sample_name', 
    'donor_name', 'mit_name', 'sample_name', 'donor_slab_tile', 'hemisphere', 
    'alignment_job_id', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.Rpl', 
    'percent.Hbb', 'scDblFinder.class', 'scDblFinder.score', 'scDblFinder.weighted', 
    'scDblFinder.cxds_score', 'S.Score', 'G2M.Score', 'Phase',
    'n_fragment', 'frac_dup', 'frac_mito', 'tsse', 'Cluster_v4'
]

print("Concatenating cell metadata & joining additional annotations...")
cell_meta = {}
for k in cols_added:
    cell_meta[k] = pl.concat([x.obs[k] for x in h5ad_files.values()])

for k in cols_add:
    cell_meta[k] = pl.concat([pl.Series(anno.loc[x.obs['Cluster_v4'], k]) for x in h5ad_files.values()])

print("Concatenating the AnnDatas...")
adata = snap.AnnDataSet(
    adatas=list(h5ad_files.items()), 
    filename=outfile,
    add_key='barcoded_cell_sample_name'
)
adata.uns["species_name"] = "Marmoset"

print("Repopulating cell metadata...")
assert all(adata.obs['barcoded_cell_sample_name'] == cell_meta['barcoded_cell_sample_name'])
for k, v in cell_meta.items():
    if k == 'barcoded_cell_sample_name':
        continue
    adata.obs[k] = v

# >>> adata
# AnnDataSet object with n_obs x n_vars = 651082 x 0 backed at '/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/data/atac/atac_h5ads/_dataset.h5ads'
# contains 133 AnnData objects with keys: 'P0023_1', 'P0010_2', 'P0021_4', 'P0025_2', 'P0074_3', 'P0014_1', 'P0023_3', 'P0080_2', 'P0018_1', 'P0053_2', 'P0077_2', 'P0061_2', 'P0068_1', 'P0076_4', 'P0072_1', 'P0050_1', 'P0073_4', 'P0080_1', 'P0079_1', 'P0067_4', 'P0069_4', 'P0011_2', 'P0048_1', 'P0075_2', 'P0071_2', 'P0052_1', 'P0041_1', 'P0062_2', 'P0067_2', 'P0074_4', 'P0034_1', 'P0023_2', 'P0057_1', 'P0066_2', 'P0064_1', 'P0076_3', 'P0070_4', 'P0059_3', 'P0065_2', 'P0030_1', 'P0068_2', 'P0010_1', 'P0077_1', 'P0032_1', 'P0069_2', 'P0072_2', 'P0068_4', 'P0070_2', 'P0070_3', 'P0052_2', 'P0008_1', 'P0037_1', 'P0038_1', 'P0078_1', 'P0078_2', 'P0040_1', 'P0068_3', 'P0014_2', 'P0060_2', 'P0022_3', 'P0042_2', 'P0044_1', 'P0018_2', 'P0059_1', 'P0021_2', 'P0009_1', 'P0049_1', 'P0056_1', 'P0066_1', 'P0065_1', 'P0056_2', 'P0022_2', 'P0053_1', 'P0009_2', 'P0066_3', 'P0036_1', 'P0047_1', 'P0031_1', 'P0011_1', 'P0067_1', 'P0025_1', 'P0062_1', 'P0022_1', 'P0045_1', 'P0008_2', 'P0015_2', 'P0073_1', 'P0060_1', 'P0060_4', 'P0022_4', 'P0046_1', 'P0064_2', 'P0075_1', 'P0079_2', 'P0060_3', 'P0071_4', 'P0058_2', 'P0071_1', 'P0079_4', 'P0061_1', 'P0069_1', 'P0080_3', 'P0071_3', 'P0024_2', 'P0080_4', 'P0075_3', 'P0074_2', 'P0079_3', 'P0039_1', 'P0051_1', 'P0059_2', 'P0058_1', 'P0067_5', 'P0045_2', 'P0023_4', 'P0068_5', 'P0024_1', 'P0081_2', 'P0057_2', 'P0059_4', 'P0076_1', 'P0074_1', 'P0005_2', 'P0028_1', 'P0063_2', 'P0043_1', 'P0067_3', 'P0076_2', 'P0043_2', 'P0073_2', 'P0073_3', 'P0015_1', 'P0042_1'
#     obs: 'barcoded_cell_sample_name', 'tissue_name', 'dissociated_cell_sample_name', 'donor_name', 'mit_name', 'sample_name', 'donor_slab_tile', 'hemisphere', 'alignment_job_id', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.Rpl', 'percent.Hbb', 'scDblFinder.class', 'scDblFinder.score', 'scDblFinder.weighted', 'scDblFinder.cxds_score', 'S.Score', 'G2M.Score', 'Phase', 'n_fragment', 'frac_dup', 'frac_mito', 'tsse', 'Cluster_v4', 'Subcortex_Class_v4', 'Subcortex_Group_v4'
#     uns: 'reference_sequences', 'genome', 'AnnDataSet', 'species', 'species_name'

print("Done. Closing file...")
adata.close()

print("Script finished.")

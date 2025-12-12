## Sync data with the latest annotation
## Daisy,last updated May 6TH, 2025
print("starting script")

## load libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import sys

import anndata as ad
import scanpy as sc

from pathlib import Path
import os

import warnings
import matplotlib
import re

warnings.filterwarnings('ignore')
warnings.simplefilter('ignore')

## change plotting settings
pd.options.mode.chained_assignment = None  # default='warn'

## change settings for plotting
sc.settings.set_figure_params(frameon=False, dpi_save = 500)

## ----------------------------------------------------------- ##
outdir = sys.argv[1]
data_path = sys.argv[2]
annotation_path = sys.argv[3]
metadata_path = sys.argv[4]

outdir = "/jukebox/krienen/marm_hmba_integration/250204_reprocess_and_recluster/BG_preprint/20250506"
# data_path = "/jukebox/krienen/marm_hmba_integration/250204_reprocess_and_recluster/BG_preprint/rna_bg_marm_20250430_preprint_BG.h5ad"
data_path = "/jukebox/krienen/marm_hmba_integration/250204_reprocess_and_recluster/BG_preprint/rna_bg_marm_20250506_preprint_BG.h5ad"
annotation_path = "/jukebox/krienen/marm_hmba_integration/250204_reprocess_and_recluster/BG_preprint/Marmoset_BG_annot_table_v3_preprint_20250506.csv"
metadata_path = "/jukebox/krienen/marm_hmba_integration/250204_reprocess_and_recluster/metadata/processed_subcortex_metadata.csv"

## check if outdir exists
if not os.path.exists(outdir):
    os.makedirs(outdir)

## load data
adata = ad.read_h5ad(data_path, backed='r')

## load annotation and metadata
annot_table = pd.read_csv(annotation_path)
metadata = pd.read_csv(metadata_path, index_col=0)

## double check if all the clusters in the annotation are in the data
cells_in_annot = metadata[metadata['cluster_id'].isin(annot_table['cluster_id'])].index
if cells_in_annot.isin(adata.obs['old_cell_barcode']).all():
    print("All cells in the annotation are in the data")
else:
    print("Some cells in the annotation are not in the data")
    missing_cells = cells_in_annot[~cells_in_annot.isin(adata.obs.index)]
    print(f"Number of missing cells: {len(missing_cells)}")
    sys.exit(1)

## columns in the adata to remove
remove_col = ['BG_Neighborhood','BG_Class', 'BG_Subclass', 'BG_Group','AIT117_scANVI_Cluster',
       'AIT117_scANVI_Group', 'AIT117_scANVI_Subclass', 'AIT117_scANVI_Class',
       'AIT117_scVI_Study_entropy', 'AIT117_scVI_toRemove',
       'AIT119_scANVI_Group', 'AIT119_scVI_toRemove',
       'AIT119_scVI_Study_entropy', 'Marm_Group_v2', 'Marm_scANVI_Group_v2','cluster_id', 'Group', 'Marm_BG_Cluster_v3',]

## remove columns from adata
for col in remove_col:
    if col in adata.obs.columns:
        adata.obs.drop(col, axis=1, inplace=True)
        print(f"Removed {col} from adata")
    else:
        print(f"{col} not found in adata")

## use bg_cluster_id, merge the annotation table into the adata.obs
annot_table = annot_table.set_index('cluster_id')
adata.obs = adata.obs.merge(annot_table, left_on='bg_cluster_id', right_index=True, how='left')

## rename the columns in the adata.obs
rename_dict = {
    'bg_cluster_id': 'cluster_id',
    'Group_y': 'Group',
    "Cluster_MapMyCells_Cluster_label": "Siletti_MapMyCells_Cluster_label",}
adata.obs.rename(columns=rename_dict, inplace=True)

## rearrange the columns in the adata.obs
col_order = ['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'sample_name', 'percent.mt',
       'percent.Rpl', 'percent.Hbb', 'scDblFinder.class', 'scDblFinder.score',
       'scDblFinder.weighted', 'scDblFinder.cxds_score', 'S.Score',
       'G2M.Score', 'Phase', 'donor_name','experiment_start_date',
       'donor_slab_tile', 'hemisphere', 'mit_name', 'tissue_name',
       'dissociated_cell_sample_name', 'facs_population_plan', 'facs_sorted',
       'cell_prep_type', 'study', 'expc_cell_capture',
       'enriched_cell_sample_quantity_count', 'port_well',
       'barcoded_cell_sample_label', 'RNA_library_creation_date',
       'RNA_library_prep_set', 'RNA_library_name',
       'RNA_tapestation_avg_size_bp', 'RNA_library_num_cycles',
       'RNA_lib_quantification_ng', 'RNA_library_pool_name',
       'ATAC_library_creation_date', 'ATAC_library_prep_set',
       'ATAC_library_name', 'ATAC_tapestation_avg_size_bp',
       'ATAC_library_num_cycles', 'ATAC_lib_quantification_ng',
       'ATAC_library_pool_name', 'cDNA_amplification_method',
       'cDNA_amplification_date', 'amplified_cdna_name', 'cDNA_pcr_cycles',
       'percent_cdna_longer_than_400bp', 'cdna_amplified_quantity_ng',
       'cDNA_library_input_ng', 'RNA_r1_index', 'RNA_r2_index', 'ATAC_index',
       'alignment_job_id', 'cell_barcode', 'alignment_job_database',
       'cell_label', 'old_cell_barcode', 'cluster_id', 'accession_cluster',
       'Neighborhood', 'Class', 'Subclass', 'Group',
       'Siletti_MapMyCells_Supercluster_label',
       'Siletti_MapMyCells_Cluster_label',
       'Siletti-reannot_MapMyCells_Class_label',
       'Siletti-reannot_MapMyCells_Subclass_label',
       'Siletti-reannot_MapMyCells_Cluster_label',
       'ABCmouse_MapMyCells_CLAS_label', 'ABCmouse_MapMyCells_SUBC_label',
       'ABCmouse_MapMyCells_CLUS_label']

adata.obs = adata.obs[col_order]

## add umaps plots to visually check the annotation
cols_plot = ['cluster_id', 'Group', 'Subclass', 'Class', 'Neighborhood']
for col in cols_plot:
    adata.obs[col] = adata.obs[col].astype('str')

os.chdir(outdir)
embeddings = ['X_umap_BG_integrated','X_umap_donor']
for each_embedding in embeddings:
    if each_embedding in adata.obsm.keys():
        sc.pl.embedding(adata, basis = each_embedding, color='cluster_id', show=False, palette = list(matplotlib.colors.CSS4_COLORS.keys()), save = "_cluster_id.pdf") ## 594 BG clusters
        sc.pl.embedding(adata, basis = each_embedding, color='Group', show=False, save = "_Group.pdf")
        sc.pl.embedding(adata, basis = each_embedding, color='Subclass',  show=False, save = "_Subclass.pdf")
        sc.pl.embedding(adata, basis = each_embedding, color='Class', show=False, save = "_Class.pdf")
        sc.pl.embedding(adata, basis = each_embedding, color='Neighborhood',  show=False, save = "_Neighborhood.pdf")
    else:
        print(f"{each_embedding} not found in adata.obsm")
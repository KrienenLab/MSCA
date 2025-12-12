#!/usr/bin/env python

import pandas as pd
from pathlib import Path

pdir                = Path("/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster")
path_umap           = Path(pdir, "metadata/rna_umap_cluster_0715_clean.csv")
path_anno           = Path(pdir, "metadata/subcortex_anno_table.csv")
# path_multiome_cmeta = Path(pdir, "data/rna/rna_cell_metadata.csv")           # cell metadata
path_multiome_lmeta = Path(pdir, "metadata/merged_datalog_metadata.csv")     # library metadata
path_out            = Path(pdir, "metadata/rna_umap_and_metadata_for_cytosplore.csv")


# --------------------------------------------------
# Imports
# --------------------------------------------------

umap  = pd.read_csv(path_umap, index_col=0)
anno  = pd.read_csv(path_anno, index_col=0) # cluster ID
# cmeta = pd.read_csv(path_multiome_cmeta, index_col=0)
lmeta = pd.read_csv(path_multiome_lmeta)
lmeta.index = lmeta['barcoded_cell_sample_name']


# --------------------------------------------------
# Verify all CBCs accounted for
# --------------------------------------------------

umap['barcoded_cell_sample_name'] = umap.index.str.replace(".*-", "", regex=True)
assert all(umap.barcoded_cell_sample_name.isin(lmeta.barcoded_cell_sample_name)), "missing/unaligned library metadata"
# assert all(umap.barcoded_cell_sample_name.isin(cmeta.barcoded_cell_sample_name)), "missing/unaligned cell metadata"


# --------------------------------------------------
# [SKIP] Get cell-level metadata
# --------------------------------------------------

# >>> cmeta.columns
# Index(['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'barcoded_cell_sample_name',
#        'sample_name', 'percent.mt', 'percent.Rpl', 'percent.Hbb',
#        'scDblFinder.class', 'scDblFinder.score', 'scDblFinder.weighted',
#        'scDblFinder.cxds_score', 'S.Score', 'G2M.Score', 'Phase'],
#       dtype='object')


## --> nevermind, i actually won't keep any


# --------------------------------------------------
# Get library-level metadata
# --------------------------------------------------

# >>> lmeta.columns
# Index(['sample_name', 'experiment_start_date', 'donor_slab_tile', 'hemisphere',
#        'mit_name', 'donor_name', 'tissue_name', 'dissociated_cell_sample_name',
#        'facs_population_plan', 'facs_sorted', 'cell_prep_type', 'study',
#        'expc_cell_capture', 'enriched_cell_sample_quantity_count', 'port_well',
#        'barcoded_cell_sample_name', 'RNA_library_creation_date',
#        'RNA_library_prep_set', 'RNA_library_name',
#        'RNA_tapestation_avg_size_bp', 'RNA_library_num_cycles',
#        'RNA_lib_quantification_ng', 'RNA_library_pool_name',
#        'ATAC_library_creation_date', 'ATAC_library_prep_set',
#        'ATAC_library_name', 'ATAC_tapestation_avg_size_bp',
#        'ATAC_library_num_cycles', 'ATAC_lib_quantification_ng',
#        'ATAC_library_pool_name', 'cDNA_amplification_method',
#        'cDNA_amplification_date', 'amplified_cdna_name', 'cDNA_pcr_cycles',
#        'percent_cdna_longer_than_400bp', 'cdna_amplified_quantity_ng',
#        'cDNA_library_input_ng', 'RNA_r1_index', 'RNA_r2_index', 'ATAC_index',
#        'alignment_job_id'],
#       dtype='object')

cols_keep_lmeta = [
    "sample_name",
    "donor_slab_tile",
    "hemisphere",
    "mit_name",
    "donor_name",
    "tissue_name"
]


# --------------------------------------------------
# Get annotations
# --------------------------------------------------

cols_keep_anno = ['Subcortex_Class_v4', 'Subcortex_Group_v4', 'NT']


# --------------------------------------------------
# Join and export
# --------------------------------------------------

out = umap.join(anno.loc[:, cols_keep_anno], on='Cluster_v4')
out = out.join(lmeta.loc[:, cols_keep_lmeta], on='barcoded_cell_sample_name')
out.to_csv(path_out)

#!/usr/bin/env python3
import anndata as ad
import pandas as pd
from pathlib import Path
import sys
from scipy.io import mmread
from scipy.sparse import csr_matrix, vstack

# paths
pdir = Path('/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster')
path_curated_rna = Path(pdir, 'data/rna')
path_cell_meta = Path(path_curated_rna, 'rna_cell_metadata.csv')
path_gene_meta = Path(path_curated_rna, 'rna_gene_metadata.csv')
path_library_meta = Path(pdir, "metadata/merged_datalog_metadata.csv")
out_path_h5ad = Path(path_curated_rna, 'rna_raw.h5ad')

# metadata
cell_meta = pd.read_csv(path_cell_meta, index_col = 0)
gene_meta = pd.read_csv(path_gene_meta, index_col = 0)
library_meta = pd.read_csv(path_library_meta, index_col = 0)

if all(gene_meta.iloc[:, 0].isna()):
    gene_meta = gene_meta.drop(gene_meta.columns[0], axis=1)

# ------------------------------
# modify cell_meta to have columns of interest
# ------------------------------

# >>> cell_meta.columns
# Index(['orig.ident', 'nCount_RNA', 'nFeature_RNA', 'barcoded_cell_sample_name',
#        'sample_name', 'percent.mt', 'percent.Rpl', 'percent.Hbb',
#        'scDblFinder.class', 'scDblFinder.score', 'scDblFinder.weighted',
#        'scDblFinder.cxds_score', 'S.Score', 'G2M.Score', 'Phase'],
#       dtype='object')

cols_add = [
    'barcoded_cell_sample_name',
    'donor_name', 
    'donor_slab_tile', 
    'hemisphere',
    'mit_name', 
    'tissue_name', 
    'dissociated_cell_sample_name',
    'alignment_job_id'
]
library_meta = library_meta.loc[:, cols_add]
library_meta.set_index('barcoded_cell_sample_name', inplace=True)
cell_meta = cell_meta.join(library_meta, on = 'barcoded_cell_sample_name')

col_order = [
    'barcoded_cell_sample_name', 
    'tissue_name', 
    'dissociated_cell_sample_name', 
    'donor_name', 
    'mit_name', 
    'sample_name', 
    'donor_slab_tile', 
    'hemisphere', 
    'alignment_job_id',
    'nCount_RNA', 
    'nFeature_RNA', 
    'percent.mt', 
    'percent.Rpl', 
    'percent.Hbb', 
    'scDblFinder.class', 
    'scDblFinder.score', 
    'scDblFinder.weighted', 
    'scDblFinder.cxds_score', 
    'S.Score', 
    'G2M.Score', 
    'Phase'
]
cell_meta = cell_meta.loc[:, col_order]

# ------------------------------
# make, export h5ad
# ------------------------------

# count matrices
n_chunks = len(list(Path.glob(path_curated_rna, "rna_counts_*.mtx")))
counts = csr_matrix(mmread(Path(path_curated_rna, "rna_counts_1.mtx")))

for i in range(2, n_chunks+1):
    if cell_meta.sample_name.iloc[counts.shape[0]-1] == cell_meta.sample_name.iloc[counts.shape[0]]:
        sys.exit("Chunked RNA count matrix does not appear to match the cell metadata!")
    pth = Path(path_curated_rna, f"rna_counts_{i}.mtx")
    counts = vstack([counts, csr_matrix(mmread(pth))])

adata = ad.AnnData(counts, obs = cell_meta, var = gene_meta)
adata.write_h5ad(out_path_h5ad)

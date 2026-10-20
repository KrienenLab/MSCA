## Run NSForest algorithm to find minimal marker genes for each of the scrattch clusters
## Daisy Dan
## Last updated June 16th, 2025

# import libraries
import sys
import os
sys.path.insert(0, os.path.abspath("../"))
sys.path.insert(0, os.path.abspath("../nsforest/nsforesting"))
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import nsforest as ns
from nsforest import utils
from nsforest import preprocessing as pp
from nsforest import nsforesting
from nsforest import evaluating as ev
from nsforest import plotting as pl
from pathlib import Path

## load data
pdir = "/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/"
data_path = Path(pdir, 'data', 'rna','rna_raw_v4_freeze.h5ad')
scvi_path = Path(pdir, 'data', 'rna','scvi_model_donor','X_scVI.csv')
adata = sc.read(data_path)
X_scVI = pd.read_csv(scvi_path, index_col=0)
X_scVI = X_scVI.loc[adata.obs.index]
adata.obsm['X_scVI'] = X_scVI.values

## define which column for cell type annotation
cluster_header = 'Cluster_v4'
output_folder = Path(pdir, 'data', 'rna','nsforest')

# prepare data for NSForest
adata.obs[cluster_header] = pd.Categorical(adata.obs[cluster_header].astype(str))

## preprocess data
adata.layers['UMIs'] = adata.X.copy()  # store raw counts in a layer

ns.pp.dendrogram(adata, cluster_header, use_rep="X_scVI", save = True, output_folder = output_folder, outputfilename_suffix = cluster_header)

## calculate cluster medians per gene
adata = ns.pp.prep_medians(adata, cluster_header)

## calculate binary scores per gene per cluster
adata = ns.pp.prep_binary_scores(adata, cluster_header)

print(f"Saving new anndata object as...\n{filename}")
adata.write_h5ad(Path(output_folder, 'rna_raw_v4_nsforest.h5ad'), compression = "gzip")

print("finished preprocessing data for NSForest")

# run NSForest
outputfilename_prefix = "Cluster"
results = nsforesting.NSForest(adata, cluster_header, output_folder = output_folder, outputfilename_prefix = outputfilename_prefix)

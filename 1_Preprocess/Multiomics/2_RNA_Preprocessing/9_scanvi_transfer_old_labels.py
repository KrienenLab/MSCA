import os, sys
import scanpy as sc
import numpy as np
import pandas as pd
import scvi
from scipy.sparse import csr_matrix
from scvi.model import SCANVI
# from scvi.inference import UnsupervisedTrainer, SemiSupervisedTrainer
from sklearn.preprocessing import LabelEncoder
import click
import scipy
import re
import matplotlib
from pathlib import Path

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

##
pdir = "/jukebox/krienen/marm_hmba_integration/250204_reprocess_and_recluster"
path_h5ad = Path(pdir, 'data', 'transcriptomic_clustering','rna_clustered.h5ad')
tmp_dir = '/scratch/sdan/temp'
scvi_dir = Path(pdir, 'data', 'scvi_model_donor','scvi_model_donor')
old_annot_path = '/jukebox/krienen/marm_hmba_integration/240905_reprocess_and_recluster/metadata/250218_subcortex_metadata.csv'

outdir = Path(pdir, 'data','scanvi_transfer')

if not Path.is_dir(outdir):
    os.mkdir(outdir)
os.chdir(outdir)


max_epochs_scanvi = 30

## ----------------------------------------
## Load data
print("Loading data...")
adata = sc.read_h5ad(path_h5ad)

## load labels from the old data (241213)
old_annot = pd.read_csv(old_annot_path)
old_annot['cell_id'] = old_annot['Unnamed: 0'].str.split('_').str[0]
old_annot['cell_index'] = old_annot['sample_name'] + '_' + old_annot['cell_id']
old_annot = old_annot[old_annot['cell_index'].isin(adata.obs.index)]
old_annot = old_annot.set_index('cell_index')

## put the old annotation onto the new data
adata.obs.loc[old_annot.index, 'Group_2024'] = old_annot['Group_v2']

## clean up the annotation
adata.obs['Group_2024'] = adata.obs['Group_2024'].astype(str)
adata.obs['Group_2024'][adata.obs['Group_2024'].isin(['nan','Split','low_quality','Unassigned','<NA>'])] = pd.NA


labels_key = "labels_scanvi"
adata.obs[labels_key] = adata.obs["Group_2024"].values
adata.obs[labels_key] = adata.obs[labels_key].fillna("Unknown")

adata_hvg = adata[:, adata.var.highly_variable].copy()

click.echo("Loading trained scVI donor model...")
scvi_model = scvi.model.SCVI.load(scvi_dir, adata_hvg) 

click.echo("Training scANVI...")
scanvi_model = SCANVI.from_scvi_model(scvi_model, labels_key=labels_key, unlabeled_category = "Unknown")
scanvi_model.train(max_epochs=max_epochs_scanvi, n_samples_per_label=100)

## save the scanvi model
click.echo("Saving trained scANVI...")
scanvi_model.save(dir_path="transfer_old_annot", save_anndata=False, overwrite=True)

SCANVI_LATENT_KEY = "X_scANVI"
adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation()

## save the scanvi data
X_scANVI_df = pd.DataFrame(adata.obsm[SCANVI_LATENT_KEY], index=adata.obs_names)
X_scANVI_df.to_csv("X_scANVI_df_cluster.csv")

sc.pp.neighbors(adata, use_rep=SCANVI_LATENT_KEY)
X_umap_integrated_scanvi = sc.tl.umap(adata, min_dist=0.3, copy=True)
adata.obsm["X_umap_scanvi_old_annot"] = X_umap_integrated_scanvi.obsm["X_umap"]

## save umap
UMAP_df = pd.DataFrame(adata.obsm['X_umap_scanvi_old_annot'], index=adata.obs_names)
UMAP_df.to_csv('X_umap_scanvi_old_annot.csv')

# ## MDE
# SCANVI_MDE_KEY = "X_scANVI_MDE"
# adata.obsm[SCANVI_MDE_KEY] = scvi.model.utils.mde(adata.obsm[SCANVI_LATENT_KEY], accelerator="gpu")

# ## save mde
# MDE_df = pd.DataFrame(adata.obsm[SCANVI_MDE_KEY], index=adata.obs_names)
# MDE_df.to_csv('X_scANVI_MDE_cluster.csv')

SCANVI_PREDICTIONS_KEY = "predictions_scanvi"
adata.obs[SCANVI_PREDICTIONS_KEY] = scanvi_model.predict(adata_hvg)
adata.obs[SCANVI_PREDICTIONS_KEY].to_csv("predictions_scanvi_cluster.csv", index = True)

click.echo("Done!")
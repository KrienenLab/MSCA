print("Job starting...")

import anndata as ad
import scanpy as sc
import scvi
import pandas as pd
import numpy as np
from pathlib import Path
import os
import matplotlib.pyplot as plt
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

# paths
pdir = Path('/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster')
path_h5ad = Path(pdir, 'data/rna', 'rna_raw.h5ad')
outdir = Path(pdir, 'data/rna', 'scvi_model_donor')

if not Path.is_dir(outdir):
    os.mkdir(outdir)

print("Importing h5ad...")
adata = ad.read_h5ad(path_h5ad)

## use the sample name to get the donor name
adata.obs['donor_name'] = adata.obs['sample_name'].str.split('_').str[2]

# adata.X = adata.X.astype(np.int32) ## Data is float64 from NeMO which is excessive for counts.

## Normalize the count matrix after storing the raw counts in the raw slot
adata.raw = adata.copy()
print("Normalizing data")
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)

print("Finding highly variable genes")
adata.layers["UMIs"] = adata.raw.X ## This is a shallow copy, no extra space is used!
sc.pp.highly_variable_genes(
    adata, 
    n_top_genes=4000, 
    layer="UMIs", 
    subset=False, 
    flavor="seurat_v3", 
    batch_key="donor_name")

adata_hvg = adata[:, adata.var.highly_variable].copy()

# print("Saving merged h5ad with HVGs selected...")

# adata.write_h5ad(Path(outdir, "adata.h5ad"))
## re-import
# adata = ad.read_h5ad(Path(outdir, "adata.h5ad"))


# --------------------------------------------------
# scVI model
# --------------------------------------------------

## Run scVI with known confounders
scvi.model.SCVI.setup_anndata(
    adata_hvg,
    layer="UMIs",
    batch_key="donor_name", 
)

print("Training scVI")
scvi_model_donor = scvi.model.SCVI(
    adata_hvg, 
    dispersion="gene-batch", 
    n_hidden=256, 
    n_latent=64, 
    n_layers=3)

scvi_model_donor.train(accelerator = "gpu", max_epochs = 300)
scvi_model_donor.save(Path(outdir, "scvi_model_donor"))

## re-import
# scvi_model_donor = scvi.model.SCVI.load(Path(outdir, "scvi_model_donor"), adata_hvg) 


# --------------------------------------------------
# Get, save scVI embeddings
# --------------------------------------------------

adata.obsm["X_scVI"] = scvi_model_donor.get_latent_representation()

## also save the latent space as a csv file
scVI = pd.DataFrame(adata.obsm["X_scVI"], index=adata.obs.index)
scVI.to_csv(Path(outdir,"X_scVI.csv"))

print("getting neighbors")
sc.pp.neighbors(adata, use_rep = "X_scVI")

# setting min_dist=0.25 and keeping spread=1; will see, we generally prefer 
# more clumping which this would encourage
print("getting umap")
adata.obsm["X_umap_donor"] = sc.tl.umap(adata, min_dist = 0.3, copy = True).obsm["X_umap"]
umap = pd.DataFrame(adata.obsm["X_umap_donor"], index=adata.obs.index)
umap.to_csv(Path(outdir,"X_umap_donor.csv"))

# --------------------------------------------------
# Export h5ad
# --------------------------------------------------
print("exporting h5ad")
adata.write_h5ad(Path(outdir, "rna_scvi_model_donor.h5ad"), compression = "gzip")

print("Done")
## re-import
# adata = ad.read_h5ad(Path(outdir, "rna_scvi_model_donor.h5ad"))

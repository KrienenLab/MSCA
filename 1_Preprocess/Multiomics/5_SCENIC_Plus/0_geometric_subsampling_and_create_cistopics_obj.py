#!/usr/bin/env python3
import os
import scanpy as sc
import pandas as pd
import numpy as np
from geosketch import gs
from fbpca import pca 
import pycisTopic
from pathlib import Path
import pickle


# ------------------------------------------------------------
# Geosketch-based down-sampling, modality splits (Neuron /
# Nonneuron) and cisTopic object creation in uint16.
# ------------------------------------------------------------

# ------------------------------------------------------------------
# 1‒ INPUT FILES & CONSTANTS
# ------------------------------------------------------------------
PATH_RNA_FULL   = (
    "/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/"
    "data/rna/rna_raw_v4_freeze.h5ad"
)
PATH_ATAC_FULL  = (
    "/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/"
    "data/atac/atac_count_matrix/cell_by_consensus_peaks.h5ad"
)

PATH_SCVI_EMB   = (
    "/usr/people/nw8333/jukebox/scratch/subcortical_pseudobulks/raw_data/"
    "X_scVI.csv"
)
PATH_META       = (
    "/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/"
    "metadata/subcortex_anno_table.csv"
)

OUT_DIR         = Path(
    "/jukebox/krienen/victor/marm_hmba_atac_50k/data"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

LABEL_COL   = "Subcortex_Class_v4"
N_TARGET    = 50000                    # desired n per group
SHRINK_UINT16 = True                     # ATAC → uint16

# ------------------------------------------------------------ #
# 1) load data                                                  #
# ------------------------------------------------------------ #
print(" Loading AnnData objects (backed mode) …")
adata_rna  = sc.read_h5ad(PATH_RNA_FULL,  backed="r")
adata_atac = sc.read_h5ad(PATH_ATAC_FULL, backed="r")

print(" Loading scVI embedding + metadata …")
X_scvi = pd.read_csv(PATH_SCVI_EMB, index_col=0)
meta   = pd.read_csv(PATH_META,    index_col=0)

# ------------------------------------------------------------ #
# 2) attach embedding & labels to the RNA object                #
# ------------------------------------------------------------ #
adata_rna.obsm["X_scVI"] = X_scvi.loc[adata_rna.obs_names]
adata_rna.obs = adata_rna.obs.join(
    meta[[LABEL_COL, "Subcortex_Group_v4"]],
    on="Cluster_v4"
)
'''
# boolean mask: True for cells we want to keep
keep = adata_rna.obs["Subcortex_Group_v4"] != "Junk"

# subset and make a copy to avoid view warnings
adata_rna = adata_rna[keep].copy()
'''
# ------------------------------------------------------------ #
# 3) Filter an ATAC object so it only keeps peaks with >0.5 non-zero for pseudobulks and tsse > 5
# ------------------------------------------------------------ #
adata_bitmap = sc.read_h5ad(
    "/usr/people/nw8333/jukebox/scratch/subcortical_pseudobulks/intermediate_data/"
    "peak_bitmaps/peak_pb_group_level_bitmap_threshold=0.5.h5ad",
    backed='r'
)
adata_bitmap = adata_bitmap.to_memory()
# Filter peaks
peaks_mask = np.array(adata_bitmap.X.sum(axis=0)).flatten() > 0

# Match those peaks by name in the ATAC object
mask_common = adata_atac.var_names.isin(adata_bitmap.var_names[peaks_mask])

# Filter ATAC object 
adata_atac = adata_atac[:, mask_common].to_memory()

# Keep only cells with tsse > 5
adata_atac = adata_atac[adata_atac.obs["tsse"] > 5].copy()

# ------------------------------------------------------------ #
# 4) SVD on scVI space                                          #
# ------------------------------------------------------------ #
print(" Computing SVD of scVI space …")
X_lat   = adata_rna.obsm["X_scVI"].to_numpy(dtype=np.float32, copy=False)
Xc      = X_lat - X_lat.mean(axis=0)
U, s, _ = pca(Xc, k=X_lat.shape[1], raw=True)
X_svd   = U * s            # same shape as U  (n_cells × n_latent)
adata_rna.obsm["X_svd"] = X_svd

# helper ------------------------------------------------------------------ #
def geosketch_group(mask, tag):
    """Return RNA + ATAC subsets after geosketch, save files, return ATAC AnnData."""
    n_avail  = mask.sum()
    n_pick   = min(N_TARGET, n_avail)
    if n_pick < N_TARGET:
        print(f"  {tag}: only {n_avail:,} cells available – sampling {n_pick:,}")

    shared_barcodes = (
    adata_rna.obs_names[mask]                   # your biological filter
    .intersection(adata_atac.obs_names)         # in both modalities
    )

    #Geosketch needs numeric rows; get them in RNA, then convert back to barcodes:
    idx     = adata_rna.obs_names.get_indexer(shared_barcodes)
    chosen  = gs(X_svd[idx], n_pick, replace=False, seed=42)
    keep_bc = shared_barcodes[chosen]              # final ordered list

    # Slice both modalities with the very same list
    rna_sub  = adata_rna[keep_bc].to_memory()
    atac_sub = adata_atac[keep_bc].to_memory()

    # Guarantee identical order
    assert (rna_sub.obs_names == atac_sub.obs_names).all()

    # Convert string columns to categoricals before writing
    atac_sub.strings_to_categoricals()

    # optional: shrink ATAC counts to uint16
    if SHRINK_UINT16:
        atac_sub.X = atac_sub.X.astype(np.uint16)

    # write AnnData files
    fn_rna  = f"rna_{tag}_geosketch_{N_TARGET//1000}kcells_marm_subcortex.h5ad"
    fn_atac = f"atac_{tag}_geosketch_{N_TARGET//1000}kcells_marm_subcortex.h5ad"
    rna_sub.write( OUT_DIR / fn_rna )
    atac_sub.write(OUT_DIR / fn_atac)
    print(f"    – {fn_rna:55s}  n={rna_sub.n_obs:,}")
    print(f"    – {fn_atac:55s}  n={atac_sub.n_obs:,}")

    return atac_sub, tag

# ------------------------------------------------------------ #
# 5) run group-specific geosketch                               #
# ------------------------------------------------------------ #
print(f" Geosketching {N_TARGET:,} cells per group …")
nonneuron_mask = adata_rna.obs[LABEL_COL] == "Nonneuron"

neuron_labels = [
    "BG", "CNU-HTH GABA", "CNU-HTH Glut", "CTX GABA", "HTH GABA", "HTH Glut",
    "HTH Hist", "MB-HB GABA", "MB-HB Glut", "MB-HTH-HB Sero", "Pineal Gland",
    "SEP GABA", "TEL Glut", "THM GABA", "THM Glut"
]
neuron_mask    = adata_rna.obs[LABEL_COL].isin(neuron_labels)

atac_neu,  tag_neu  = geosketch_group(neuron_mask, "neuron")
atac_non,  tag_non  = geosketch_group(nonneuron_mask, "nonneuron")

# ------------------------------------------------------------ #
# 6) build cisTopic objects (uint16)          #
# ------------------------------------------------------------ #
from scipy.sparse import csr_matrix
from pycisTopic.cistopic_class import *

print(" Building cisTopic objects …")
def build_cisto(ad_atac, tag):
    X_csr = csr_matrix(ad_atac.X.T, dtype=np.uint16)
    cisto = create_cistopic_object(
        X_csr,
        cell_names   = ad_atac.obs_names.tolist(),
        region_names = ad_atac.var_names.tolist()
    )

    fname = f"cistopic_obj_{tag}_geosketch_{N_TARGET//1000}kcells_marm_subcortex.pkl"
    with open(OUT_DIR / fname, "wb") as fh:
        pickle.dump(cisto, fh)
    print(f"    – {fname:55s}  cells={len(cisto.cell_names)}")

build_cisto(atac_neu, tag_neu)
build_cisto(atac_non, tag_non)

print(" Finished – separate geosketch subsets and cisTopic pickles created")


'''
# Read full RNA and ATAC objects (backed mode)
adata_rna = sc.read_h5ad("/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/data/rna/rna_raw_v4_freeze.h5ad", backed='r')
adata_peaks = sc.read_h5ad("/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/data/atac/atac_count_matrix/cell_by_consensus_peaks.h5ad", backed='r')

# Load latent space and metadata
scvi_df = pd.read_csv('/usr/people/nw8333/jukebox/scratch/subcortical_pseudobulks/raw_data/X_scVI.csv', index_col=0)
metadata_df = pd.read_csv('/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/metadata/subcortex_anno_table.csv', index_col=0)

# Subset metadata of interest
meta = metadata_df[['Subcortex_Class_v4', 'Subcortex_Group_v4']].copy()

# Add metadata and latent space
adata_rna.obsm['X_scVI'] = scvi_df.loc[adata_rna.obs_names]
adata_rna.obs = adata_rna.obs.join(meta, on='Cluster_v4')

# Compute SVD on the scVI embedding 
# Centre each latent dimension  (important for SVD - PCA equivalence)
X_latent = adata_rna.obsm['X_scVI'].to_numpy(dtype=np.float32, copy=False)
Xc = X_latent - X_latent.mean(axis=0)

#Run a thin SVD with fbpca
k = X_latent.shape[1]                            # number of components
U, s, Vt = pca(Xc, k=k, raw=True)  # raw=True to skip an internal re-centering

#    U : (n_cells × k)        left singular vectors (orthonormal)
#    s : (k,)                 singular values
#    Vt: (k × n_latent_dims)  right singular vectors

# Build the reduced representation (scores) exactly like PCA does
X_dimred = U * s                 # broadcasting multiplies each column by s

# Store in the AnnData object 
adata_rna.obsm["X_svd"]  = X_dimred

# Subsample representative nuclei by geometric sketching of transcriptional latent space
#N = int(len(adata_rna)*0.20) # Target 20% of the data, in this case 139,562 nuclei
N = 100_000
sketch_index = gs(X_dimred, N, replace=False, seed=42)  
adata_rna_subset = adata_rna[sketch_index].to_memory()
adata_rna_subset

# Get RNA barcodes that actually exist in ATAC data
valid_barcodes = adata_peaks.obs_names.isin(adata_rna_subset.obs_names)

# Subset and write ATAC 
out_path = (
    "/jukebox/krienen/victor/marm_hmba_atac_downsampled/data/"
    "atac_geosketch_100kcells_marm_subcortex.h5ad"
)

adata_peaks_subset = adata_peaks[valid_barcodes, :].to_memory()

# Convert string columns to categoricals before writing
adata_peaks_subset.strings_to_categoricals()
adata_peaks_subset

# Save subset objects
adata_peaks_subset.write(out_path)
adata_rna_subset.write("/jukebox/krienen/victor/marm_hmba_atac_downsampled/data/rna_geosketch_100kcells_marm_subcortex.h5ad")

# Create cistopicsobject
from scipy.sparse import csr_matrix
from pycisTopic.cistopic_class import *

counts_peaks_subset = adata_peaks_subset.X

X = csr_matrix(counts_peaks_subset.T)
cell_names = adata_peaks_subset.obs_names.tolist()
region_names = adata_peaks_subset.var_names.tolist()

cistopic_obj = create_cistopic_object(X, cell_names=cell_names, region_names=region_names)
out_cistopic_path = Path("/jukebox/krienen/victor/marm_hmba_atac_downsampled/data") / \
           "cistopic_obj_geosketch_100kcells_marm_subcortex.pkl"

# Save cistopicsobject
with open(out_cistopic_path, "wb") as fh:
    pickle.dump(cistopic_obj, fh)

cistopic_obj
'''
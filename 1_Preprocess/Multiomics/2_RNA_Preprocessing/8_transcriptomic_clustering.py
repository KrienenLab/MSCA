#!/usr/bin/env python3
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
from pathlib import Path
import os
import shutil
import matplotlib.pyplot as plt
import json
import pickle

import transcriptomic_clustering as tc
from transcriptomic_clustering.iterative_clustering import (
    build_cluster_dict, iter_clust, OnestepKwargs, summarize_final_clusters
)

from transcriptomic_clustering.final_merging import (
    final_merge, FinalMergeKwargs,
)

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

# paths
pdir = Path('/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster')
path_h5ad = Path(pdir, 'data/rna/scvi_model_donor','rna_scvi_model_donor.h5ad')
outdir = Path(pdir, 'data/rna/transcriptomic_clustering')
tmp_dir = Path(outdir, ".tmp")

if not Path.is_dir(outdir):
    os.mkdir(outdir)

if not tmp_dir.is_dir():
    os.mkdir(tmp_dir)

os.chdir(outdir)

print("Importing h5ad...")
adata = ad.read_h5ad(path_h5ad)

## set up transcriptomic clustering parameters
def setup_transcriptomic_clustering():
    means_vars_kwargs = {
        'low_thresh': 0.6931472, # lowest value required for a gene to pass filtering. set to 1 originally, 0.6931472 to match to bigcat
        'min_cells': 4 # minimum number of cells expressed required for a gene to pass filtering
    }
    highly_variable_kwargs = {
        'max_genes': 4000 # originally 3000, 4000 to match to bigcat
    }
    pca_kwargs = {
        'cell_select': 30000, # originally 500000 cells
        'n_comps': 50,
        'svd_solver': 'randomized'
    }
    filter_pcs_kwargs = {
        'known_components': None,
        'similarity_threshold': 0.7,
        'method': 'zscore', # or elbow
        'zth': 2,
        'max_pcs': None,
    }
    ## Leave empty if you don't want to use known_modes
    filter_known_modes_kwargs = {
        # 'known_modes': known_modes_df, # a pd dataframe. index is obs (cell) names, columns are known modes. Originally commented out
        'similarity_threshold': 0.7
    }
    ## !!NEW!! Original method: "PCA", allows the user to select any obsm latent space such as "X_scVI" for leiden clustering.
    latent_kwargs = {
        # 'latent_component': "X_pca"
        'latent_component': "X_scVI"
    }
 
    cluster_louvain_kwargs = {
        'k': 15, # number of nn, originally 150, change to 15
        'nn_measure': 'euclidean',
        'knn_method': 'annoy',
        'louvain_method': 'taynaud', #'vtraag',
        'weighting_method': 'jaccard',
        'n_jobs': 30, # cpus # originally 8`
        'resolution': 1.0 # resolution of louvain for taynaud method
    }
    merge_clusters_kwargs = {
        'thresholds': {
            'q1_thresh': 0.5,
            'q2_thresh': None,
            'cluster_size_thresh': 10, ## originally uses 50, 10 to match to bigcat
            'qdiff_thresh': 0.7,
            'padj_thresh': 0.05,
            'lfc_thresh': 0.6931472, # log2 fold change threshold for DE genes
            'score_thresh': 100, # originally uses 200, 100 to match to bigcat
            'low_thresh': 0.6931472, # originally uses 1 # applied to log2(cpm+1) to determine if a gene is expressed or not, 0.6931472 to match to bigcat
            'min_genes': 5
        },
        'k': 4, # number of nn for de merge, originaly 2, 4 to match to bigcat
        'de_method': 'ebayes'
    }
 
    onestep_kwargs = OnestepKwargs(
        means_vars_kwargs = means_vars_kwargs,
        highly_variable_kwargs = highly_variable_kwargs,
        pca_kwargs = pca_kwargs,
        filter_pcs_kwargs = filter_pcs_kwargs,
        filter_known_modes_kwargs = filter_known_modes_kwargs,
        latent_kwargs = latent_kwargs,
        cluster_louvain_kwargs = cluster_louvain_kwargs,
        merge_clusters_kwargs = merge_clusters_kwargs
    )
    return onestep_kwargs


## =================================================================================================== ##
onestep_kwargs = setup_transcriptomic_clustering()

## Run clustering
print("Running clustering")
clusters, markers = iter_clust(
    adata,
    min_samples=4,
    onestep_kwargs=onestep_kwargs,
    random_seed=123,
    tmp_dir=tmp_dir
)
with open(os.path.join(outdir, 'clustering_results_0305.pkl'), 'wb') as f:
    pickle.dump(clusters, f)

##
with open(os.path.join(outdir,'clustering_markers_0305.pkl'), 'wb') as f:
    pickle.dump(markers, f)

## load the .pickle files
with open(os.path.join(outdir, 'clustering_results_0305.pkl'), 'rb') as f:
    clusters = pickle.load(f)

with open(os.path.join(outdir,'clustering_markers_0305.pkl'), 'rb') as f:
    markers = pickle.load(f)


## Convert to a list of lists
clusters_as_lists = [array.tolist() for array in clusters]

def setup_merging(): 
    merge_clusters_kwargs = {
        'thresholds': {
            'q1_thresh': 0.5,
            'q2_thresh': None,
            'cluster_size_thresh': 10, 
            'qdiff_thresh': 0.7, 
            'padj_thresh': 0.05, 
            'lfc_thresh': 0.6931472, 
            'score_thresh': 100, 
            'low_thresh': 0.6931472, 
            'min_genes': 5
        },
        'k': 4,
        'de_method': 'ebayes',
    }
    latent_kwargs = {  
        "latent_component": "X_scVI"
    }
    merge_kwargs = FinalMergeKwargs(
        merge_clusters_kwargs=merge_clusters_kwargs,
        latent_kwargs=latent_kwargs  
    )
    return merge_kwargs

merge_kwargs = setup_merging()

## Run the final merging
clusters_after_merging, markers = final_merge(
    adata, 
    clusters_as_lists, 
    # markers, # required for PCA, but optional if using a pre-computed latent space
    n_samples_per_clust=20, 
    random_seed=123, 
    # n_jobs = 30, # modify this to the number of cores you want to use
    n_markers = None, ## skip calculating markers
    # return_markers_df = False, # return the pair-wise DE results for each cluster pair. If False (default), only return a set of markers (top 20 of up and down regulated genes in each pair comparison)
    final_merge_kwargs=merge_kwargs
)

##
print("finished clustering")

with open(os.path.join(outdir, 'clustering_results_final_merging_0305.pkl'), 'wb') as f:
    pickle.dump(clusters_after_merging, f)

##
# with open(os.path.join(outdir,'clustering_markers_final_merging.pkl'), 'wb') as f:
#     pickle.dump(markers, f)

## Gather clustering results
print("Summarizing final clusters")
n_cells = sum(len(i) for i in clusters_after_merging)
print(f"Number of cells: {n_cells}")
cl = ['unknown']*n_cells

for i in range(len(clusters_after_merging)):
    for j in clusters_after_merging[i]:
        cl[j] = i+1

##
print("Saving clusters")
clusters = pd.DataFrame({'cl': cl}, index=adata.obs_names)
np.all(adata.obs_names == clusters.index)


# ##
# cluster_dict = build_cluster_dict(clusters)

# ##
# print("Saving clusters")
# adata.obs["cluster"] = ""
# for cluster in cluster_dict.keys():
#     adata.obs.cluster[cluster_dict[cluster]] = cluster

# ##
adata.obs['cluster'] = clusters['cl']
adata.obs['cluster'] = adata.obs['cluster'].astype("category")

## save the clusters
clusters = pd.DataFrame(adata.obs["cluster"], index = adata.obs.index)
clusters.to_csv("ts_cluster.csv")

# with open('clusters.json', 'w') as f:
#         json.dump(cluster_dict, f)

##
print("Saving clustered data")
adata.write(Path(outdir,"rna_clustered.h5ad"), compression="gzip")

shutil.rmtree(tmp_dir)
print("Script finished.")

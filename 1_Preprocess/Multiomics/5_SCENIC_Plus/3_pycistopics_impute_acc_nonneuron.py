#!/usr/bin/env python3
import os
import scanpy as sc
import numpy as np
import anndata as ad
import pandas as pd
from pathlib import Path
import pickle
import pycisTopic
pycisTopic.__version__

os.chdir("/jukebox/krienen/victor/marm_hmba_atac_50k/pycisTopic/nonneuron")

# Set dir to store the output of pycisTopic
out_dir = "outs"
os.makedirs(out_dir, exist_ok = True)

with open("/jukebox/krienen/victor/marm_hmba_atac_50k/data/cistopic_obj_nonneuron_geosketch_50kcells_marm_subcortex.pkl", "rb") as fp:
    cistopic_obj = pickle.load(fp)

from pycisTopic.topic_binarization import binarize_topics

region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, method='ntop', ntop = 3_000,
    plot=True, num_columns=5
)

region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=True, num_columns=5
)

binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=True,
    num_columns=5, nbins=100)


from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)

'''# Get a boolean peaks that are used in over half of the cells of at least one group
adata_bitmap = sc.read_h5ad("/usr/people/nw8333/jukebox/scratch/subcortical_pseudobulks/intermediate_data/peak_bitmaps/peak_pb_group_level_bitmap_threshold=0.5.h5ad", backed='r')
adata_bitmap = adata_bitmap.to_memory()
peaks_mask = np.array(adata_bitmap.X.sum(axis=0)).flatten() > 0
peak_names = adata_bitmap.var_names[peaks_mask]'''

imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    scale_factor=10**6,
    chunk_size=5000
)


# Update pycistopic object after processing 
pickle.dump(
    imputed_acc_obj,
    open(os.path.join(out_dir, "imputed_acc_obj_geosketch_50kcells_marm_subcortex_nonneuron.pkl"),"wb")
)
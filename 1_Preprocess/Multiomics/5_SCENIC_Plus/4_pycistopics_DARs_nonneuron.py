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

with open("/jukebox/krienen/victor/marm_hmba_atac_50k/pycisTopic/nonneuron/outs/imputed_acc_obj_geosketch_50kcells_marm_subcortex_nonneuron.pkl", "rb") as fp:
    imputed_acc_obj = pickle.load(fp)

from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)

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

normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

'''adata_atac = sc.read_h5ad("/jukebox/krienen/victor/marm_hmba_atac_downsampled/data/atac_nonneuron_geosketch_50kcells_marm_subcortex.h5ad", backed='r')
metadata = adata_atac.obs[['Subcortex_Class_v4', 'Subcortex_Group_v4']].copy()
cistopic_obj.add_cell_data(metadata, split_pattern='___')'''

variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True
)

print(len(variable_regions))

markers_dict= find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='Subcortex_Group_v4',
    var_features=variable_regions,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=5,
    _temp_dir='/jukebox/scratch/vn0027/',
    split_pattern = '___'
)

print("Number of DARs found:")
print("---------------------")
for x in markers_dict:
    print(f"  {x}: {len(markers_dict[x])}")

# Save region sets
os.makedirs(os.path.join(out_dir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "region_sets", "DARs_cell_type"), exist_ok = True)

from pycisTopic.utils import region_names_to_coordinates

for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )

for cell_type in markers_dict:
    region_names_to_coordinates(
        markers_dict[cell_type].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(out_dir, "region_sets", "DARs_cell_type", f"{cell_type}.bed"),
        sep = "\t",
        header = False, index = False
    )

# Save DAR dictionary and imputed accessibility object
pickle.dump(
    markers_dict,
    open(os.path.join(out_dir, "markers_dict.pkl"), "wb")
)

pickle.dump(
    variable_regions,
    open(os.path.join(out_dir, "variable_regions.pkl"), "wb")
)

# Update pycistopic object after processing 
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj_geosketch_50kcells_marm_subcortex_nonneuron.pkl"),"wb")
)
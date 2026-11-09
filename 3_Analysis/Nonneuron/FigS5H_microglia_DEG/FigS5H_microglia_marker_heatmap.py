#!/usr/bin/env python
"""Fig. S5H: mean expression of microglia marker genes in each microglia cluster (heatmap).

Usage: python FigS5H_microglia_marker_heatmap.py [path/to/nonneuron_cells.h5ad]

Input   cell-level h5ad of the non-neuronal cells (default data/rna_nonneurons.h5ad); .X = normalized
        expression, obs columns 'Subcortex_Group_v4' and 'Cluster_v4'. Only the Microglia group is used.
Output  output/microglia_marker_genes_heatmap.svg
"""

import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

script_dir = Path(__file__).parent.resolve()
adata_path = Path(sys.argv[1]) if len(sys.argv) > 1 else script_dir / "data" / "rna_nonneurons.h5ad"
out_dir = script_dir / "output"
out_dir.mkdir(exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 400,
    "font.family": "sans-serif",
    "lines.linewidth": 0.3,
    "font.size": 10,
    "font.weight": "normal",
    "svg.fonttype": "none",
    "axes.labelsize": 10,
    "axes.linewidth": 0.3,
    "axes.titlesize": 10,
    "xtick.major.size": 2,
    "xtick.major.width": 1,
    "ytick.major.size": 2,
    "ytick.major.width": 1,
    "savefig.format": "svg",
})

MARKER_GENES = [
    "TREM2", "CX3CR1",                   # homeostatic microglia
    "CD83", "HIF1A",                     # transitional microglia
    "ITGAX", "SPP1",                     # disease-associated microglia
    "PCSK2", "JAZF1", "HCRTR2", "KCNQ3", "TEAD1", "VCL", "CNRIP1", "PNPLA7", "ST6GALNAC3", "NCKAP5",  # BST clusters
]
BST_CLUSTERS = ["cluster_857", "cluster_854", "cluster_855", "cluster_858"]  # shown first
COLOR_MAX = 5

# ----------------------------------------
# Mean expression of the marker genes per microglia cluster
# ----------------------------------------
adata = ad.read_h5ad(adata_path, backed="r")
microglia = adata[adata.obs["Subcortex_Group_v4"].isin(["Microglia"])].to_memory()
microglia.obs["Cluster_v4"] = microglia.obs["Cluster_v4"].astype("category").cat.remove_unused_categories()

x = microglia[:, MARKER_GENES].X
expression = pd.DataFrame(x.toarray() if hasattr(x, "toarray") else np.asarray(x),
                          index=microglia.obs_names, columns=MARKER_GENES)
mean_expression = expression.groupby(microglia.obs["Cluster_v4"], observed=True).mean().T  # genes x clusters

cluster_order = BST_CLUSTERS + [c for c in microglia.obs["Cluster_v4"].cat.categories if c not in BST_CLUSTERS]
mean_expression = mean_expression[cluster_order]

# ----------------------------------------
# Heatmap
# ----------------------------------------
fig, ax = plt.subplots(figsize=(3.5, 3))
image = ax.imshow(mean_expression, aspect="auto", cmap="magma", rasterized=True, vmax=COLOR_MAX)
ax.set_xticks(np.arange(mean_expression.shape[1]))
ax.set_xticklabels([c.replace("cluster_", "") for c in mean_expression.columns], rotation=90, ha="center")
ax.set_yticks(np.arange(mean_expression.shape[0]))
ax.set_yticklabels(mean_expression.index)
ax.set_xlabel("Microglia Cluster")
ax.set_title("Microglia Marker Genes Mean Expression")
fig.colorbar(image, ax=ax).set_label("Mean Expression")
fig.tight_layout()
fig.savefig(out_dir / "microglia_marker_genes_heatmap.svg", dpi=300, transparent=True)
plt.close(fig)

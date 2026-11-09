#!/usr/bin/env python
"""Fig. 4B: PCA of the astrocyte cluster metacells and the genes that drive PC1.

Usage: python Fig4B_Astrocyte_PCA.py [path/to/nonneuron_cells.h5ad]

Inputs (./data)
    Astro_cluster_metacell_log1p_hvg.csv   log1p expression of the highly variable genes in the metacells of
                                           the astrocyte clusters (rows = clusters)
    S4 - Cluster Info.xlsx                 cluster -> group (sheet 'ST4A Cluster Info')
Optional cell-level h5ad of the non-neuronal cells (first argument; .X = normalized expression, obs column
'Subcortex_Group_v4'), only for the dot plot of the top PC1 genes; the dot plot is skipped without it.

Outputs (./output)
    Astro_metacell_pca.svg                 PC1 vs PC2, colored by group
    Astro_PC1_loadings.csv                 PC1 loading of every gene, highest first
    astro_top_pc1_genes_dotplot.svg        dot plot of TOP_PC1_GENES by group (needs the cell-level h5ad)
The top PC1 genes are printed.
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA

script_dir = Path(__file__).parent.resolve()
data_dir = script_dir / "data"
out_dir = script_dir / "output"
out_dir.mkdir(exist_ok=True)
cells_path = Path(sys.argv[1]) if len(sys.argv) > 1 else None

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

N_COMPONENTS = 10
TOP_N = 30  # number of top PC1 genes to print
GROUP_COLORS = {"Astro-TE NN": "#6BC9B2", "Astro-WM NN": "#E09E50", "Astro-NT NN": "#A260C6"}
# genes with the highest PC1 loadings, shown in the dot plot
TOP_PC1_GENES = ["MOXD1", "GRM3", "ADGRV1", "FOXG1", "MID1"]
DOTPLOT_GROUP_ORDER = ["Astro-TE NN", "Astro-WM NN", "Astro-NT NN"]

# ----------------------------------------
# Metacells of the astrocyte clusters and their groups
# ----------------------------------------
cluster_info = pd.read_excel(data_dir / "S4 - Cluster Info.xlsx", sheet_name="ST4A Cluster Info", skiprows=35)
cluster_to_group = cluster_info.set_index("Cluster_id")["Subcortex_Group_v4"]

expression = pd.read_csv(data_dir / "Astro_cluster_metacell_log1p_hvg.csv", index_col=0)
groups = expression.index.to_series().map(cluster_to_group)
if groups.isna().any() or not groups.isin(GROUP_COLORS).all():
    raise ValueError(f"Clusters without an astrocyte group in S4: {groups[~groups.isin(GROUP_COLORS)].index.tolist()}")

# ----------------------------------------
# PCA
# ----------------------------------------
pca = PCA(n_components=N_COMPONENTS)
scores = pca.fit_transform(expression.values)
print(f"Explained variance ratio of the {N_COMPONENTS} PCA components: {pca.explained_variance_ratio_}")

fig, ax = plt.subplots(figsize=(5, 3))
ax.scatter(scores[:, 0], scores[:, 1], s=20, c=groups.map(GROUP_COLORS))
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_title("Astro Metacell PCA")
sns.despine(ax=ax)
ax.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=GROUP_COLORS[g], markersize=8, label=g)
                   for g in groups.unique()], frameon=False, loc=(1, 0.5))
fig.tight_layout()
fig.savefig(out_dir / "Astro_metacell_pca.svg")
plt.close(fig)

# ----------------------------------------
# Genes driving PC1
# ----------------------------------------
pc1_loadings = (pd.DataFrame({"gene": expression.columns, "pc1_loading": pca.components_[0]})
                .sort_values("pc1_loading", ascending=False))
pc1_loadings.to_csv(out_dir / "Astro_PC1_loadings.csv", index=False)
print(f"Top {TOP_N} genes driving PC1:")
print(pc1_loadings["gene"].head(TOP_N).to_list())

# ----------------------------------------
# Dot plot of the top PC1 genes (cell-level data)
# ----------------------------------------
if cells_path is None or not cells_path.exists():
    print("No cell-level h5ad given: skipping the dot plot")
else:
    import scanpy as sc

    adata = sc.read_h5ad(cells_path)
    adata_astro = adata[adata.obs["Subcortex_Group_v4"].str.contains("Astro"), :].copy()
    dotplot = sc.pl.dotplot(adata_astro, var_names=TOP_PC1_GENES, groupby="Subcortex_Group_v4",
                            categories_order=DOTPLOT_GROUP_ORDER, return_fig=True)
    dotplot.savefig(out_dir / "astro_top_pc1_genes_dotplot.svg")
    plt.close("all")

#!/usr/bin/env python
"""Fig. 2A-B: mapping of marmoset subcortical groups onto human (Siletti) superclusters.

1. Barplot of the number of marmoset groups / clusters mapped to each Siletti supercluster.
2. Dotplot of the MapMyCells co-clustering: dot size = number of cells of a marmoset group assigned to a
   Siletti supercluster, dot color = Spearman correlation between the marmoset group metacell and the
   Siletti supercluster metacell (homologous genes).

Inputs (./data)
    S4 - Cluster Info.xlsx, S5 - Group Info.xlsx     supplementary tables (sheet 'ST4A Cluster Info' of S4)
    supercluster_color_hex.csv                       Siletti supercluster colors
    siletti_genes.csv                                human Ensembl ID -> gene symbol
    marm_human_homologous_gene.csv                   marmoset-human homologous gene symbols
    rna_clean_0715_cell_metadata.csv                 cell-level metadata (Cluster_v4, Siletti_MapMyCells_Supercluster_label)
    marm_Subcortex_Group_v4_metacell.h5ad            marmoset group metacells
    human_supercluster_term_metacell.h5ad            human Siletti supercluster metacells
Outputs (./output): svg figures
"""

from pathlib import Path

import anndata as ad
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

try:
    from anndata.io import read_elem
except ImportError:  # anndata < 0.11
    from anndata.experimental import read_elem

script_dir = Path(__file__).parent.resolve()
data_dir = script_dir / "data"
out_dir = script_dir / "output"
out_dir.mkdir(exist_ok=True)

NONNEURON_CLASS = "Nonneuron"
# Siletti non-neuronal superclusters that cells of neuronal groups can still be assigned to
SILETTI_NONNEURONS = ["Astrocyte", "Committed oligodendrocyte precursor", "Oligodendrocyte precursor"]
LEGEND_CELL_COUNTS = [100, 1000, 10000]  # dot sizes shown in the size legend
SUPERCLUSTER_ORDER = [  # y axis of the log2 barplot, top to bottom
    "Splatter", "Thalamic excitatory", "Amygdala excitatory", "Medium spiny neuron", "Miscellaneous",
    "CGE interneuron", "Deep-layer corticothalamic and 6b", "Midbrain-derived inhibitory",
    "Deep-layer intratelencephalic", "Hippocampal CA4", "MGE interneuron", "Mammillary body",
    "Upper-layer intratelencephalic", "LAMP5-LHX6 and Chandelier", "Hippocampal dentate gyrus",
    "Hippocampal CA1-3", "Eccentric medium spiny neuron", "Deep-layer near-projecting",
    "Upper rhombic lip", "Lower rhombic lip",
]

plt.rcParams.update({
    "figure.dpi": 100,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
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


# --------------------------------------------------
# Helpers
# --------------------------------------------------
def top_label(labels):
    """'Splatter(0.98)| Thalamic excitatory(0.02)' -> 'Splatter' (the first, highest-proportion label)."""
    return labels.str.split("|").str[0].str.extract(r"^\s*([^(]*)\(")[0].str.strip()


def read_expression(path, gene_names=None):
    """Metacell expression (.X) as a DataFrame (metacells x genes); optionally rename genes by `gene_names`."""
    # X, obs and var are read directly: .uns (not needed) can be written by a newer anndata than the one installed
    with h5py.File(path, "r") as f:
        adata = ad.AnnData(X=read_elem(f["X"]), obs=read_elem(f["obs"]), var=read_elem(f["var"]))
    if gene_names is not None:
        original = pd.Series(adata.var_names)
        adata.var_names = original.map(gene_names).fillna(original).to_numpy()
        adata.var_names_make_unique()
    x = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    return pd.DataFrame(x, index=adata.obs_names, columns=adata.var_names)


def spearman_matrix(a, b):
    """Spearman correlation between every row of `a` and every row of `b` (rows x rows, shared columns)."""
    ranks_a = a.rank(axis=1).to_numpy()
    ranks_b = b.rank(axis=1).to_numpy()
    corr = np.corrcoef(ranks_a, ranks_b)[: len(a), len(a):]
    return pd.DataFrame(corr, index=a.index, columns=b.index)


def cocluster_table(cell_metadata, cluster_info):
    """Number of cells per (Siletti supercluster, marmoset group) for the neuronal groups."""
    cells = cell_metadata.merge(cluster_info[["Cluster_id", "Subcortex_Group_v4", "Subcortex_Class_v4"]],
                                left_on="Cluster_v4", right_on="Cluster_id")
    cells = cells[cells["Subcortex_Class_v4"] != NONNEURON_CLASS]
    table = (cells.groupby(["Siletti_MapMyCells_Supercluster_label", "Subcortex_Group_v4"]).size()
             .unstack(fill_value=0))
    return table.loc[~table.index.isin(SILETTI_NONNEURONS)]


# --------------------------------------------------
# Load tables
# --------------------------------------------------
# tables start after the block describing the columns (28 rows in S5, 35 rows in S4)
group_info = pd.read_excel(data_dir / "S5 - Group Info.xlsx", skiprows=28)
cluster_info = pd.read_excel(data_dir / "S4 - Cluster Info.xlsx", sheet_name="ST4A Cluster Info", skiprows=35)

group_info = group_info[group_info["Subcortex_Class_v4"] != NONNEURON_CLASS].copy()
group_info["siletti_supercluster"] = top_label(group_info["Siletti_MapMyCells_Supercluster_label"])
supercluster_colors = pd.read_csv(data_dir / "supercluster_color_hex.csv", index_col=0)["color_hex"]

# --------------------------------------------------
# Barplots: marmoset groups / clusters per Siletti supercluster
# --------------------------------------------------
counts = (group_info.groupby("siletti_supercluster")
          .agg(n_groups=("Subcortex_Group_v4", "count"), n_clusters=("N_clusters_in_group", "sum")))
counts["supercluster"] = counts.index
palette = supercluster_colors.to_dict()

# groups and clusters side by side, log2 axes, superclusters in a fixed order (top to bottom)
if set(SUPERCLUSTER_ORDER) != set(counts.index):
    raise ValueError("SUPERCLUSTER_ORDER does not match the superclusters in S5: "
                     f"{sorted(set(SUPERCLUSTER_ORDER) ^ set(counts.index))}")
counts_ordered = counts.loc[SUPERCLUSTER_ORDER]
fig, axes = plt.subplots(1, 2, figsize=(3, 2.5), sharex=True)
for ax, metric, title in zip(axes, ["n_groups", "n_clusters"],
                             ["Number of Mapped\nMarmoset Groups", "Number of Mapped\nMarmoset Clusters"]):
    sns.barplot(data=counts_ordered, y="supercluster", x=metric, hue="supercluster", palette=palette,
                order=SUPERCLUSTER_ORDER, legend=False, ax=ax, width=0.5)
    ax.set_xscale("log", base=2)
    powers = np.unique(np.linspace(np.floor(np.log2(counts[metric].min())),
                                   np.ceil(np.log2(counts[metric].max())), 3).round().astype(int))
    superscripts = str.maketrans("0123456789-+", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺")
    ax.set_xticks(2.0 ** powers)
    ax.set_xticklabels([f"2{str(p).translate(superscripts)}" for p in powers])
    ax.tick_params(axis="x", labelsize=6)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title(title, fontsize=6)
axes[0].tick_params(axis="y", labelsize=6)
axes[1].tick_params(axis="y", left=False, labelleft=False)
plt.subplots_adjust(wspace=0.01)
fig.tight_layout()
fig.savefig(out_dir / "supercluster_group_cluster_bars_log2.svg", bbox_inches="tight")
plt.close(fig)

# --------------------------------------------------
# Dotplot: co-clustering (size) and metacell correlation (color)
# --------------------------------------------------
cell_metadata = pd.read_csv(data_dir / "rna_clean_0715_cell_metadata.csv",
                            usecols=["Cluster_v4", "Siletti_MapMyCells_Supercluster_label"])
cocluster = cocluster_table(cell_metadata, cluster_info)

siletti_genes = pd.read_csv(data_dir / "siletti_genes.csv", index_col=0)["Gene"]
homologous_genes = pd.read_csv(data_dir / "marm_human_homologous_gene.csv", index_col=0).iloc[:, 0]
marm = read_expression(data_dir / "marm_Subcortex_Group_v4_metacell.h5ad")
human = read_expression(data_dir / "human_supercluster_term_metacell.h5ad", gene_names=siletti_genes)

shared_genes = sorted(set(marm.columns) & set(human.columns) & set(homologous_genes))
correlation = spearman_matrix(human[shared_genes], marm[shared_genes])  # superclusters x groups

# groups ordered by their Siletti supercluster; every neuronal group needs metacells and cells
group_order = group_info.sort_values("siletti_supercluster", kind="stable")["Subcortex_Group_v4"]
missing = sorted(set(group_order) - set(correlation.columns))
if missing:
    raise ValueError(f"Groups without a metacell (group names must match S5): {missing}")
group_order = group_order[group_order.isin(cocluster.columns)].to_numpy()
superclusters = cocluster.index.intersection(correlation.index)

cocluster = cocluster.loc[superclusters, group_order]
correlation = correlation.loc[superclusters, group_order]
vmin, vmax = correlation.to_numpy().min(), correlation.to_numpy().max()

# the wide dotplot is split into two stacked panels
fig, axes = plt.subplots(2, 1, figsize=(9, 6.5), sharey=True)
for ax, cols in zip(axes, np.array_split(np.arange(len(group_order)), 2)):
    x, y = np.meshgrid(np.arange(len(cols)), np.arange(len(superclusters)))
    ax.scatter(x.flatten(), y.flatten(), s=np.sqrt(cocluster.to_numpy()[:, cols].flatten()),
               c=correlation.to_numpy()[:, cols].flatten(), cmap="Reds", vmin=vmin, vmax=vmax,
               alpha=0.8, edgecolors="k")
    ax.set_xlim(-0.5, len(cols) - 0.5)
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(group_order[cols], rotation=90, ha="center", fontsize=6)
    ax.set_yticks(np.arange(len(superclusters)))
    ax.set_yticklabels(superclusters)
    ax.tick_params(axis="y", labelsize=6)
    sns.despine(ax=ax)
axes[0].set_ylabel("Siletti Supercluster", fontsize=6)
fig.tight_layout()
fig.savefig(out_dir / "co_clust_dotplot_siletti_supercluster_marmoset_subcortex_group_v4_split_noleg.svg",
            dpi=500, transparent=True, bbox_inches="tight")
plt.close(fig)

# legends as separate figures
fig, ax = plt.subplots(figsize=(2.2, 1.2))
handles = [plt.scatter([], [], s=np.sqrt(c), color="gray", edgecolors="k", alpha=0.8) for c in LEGEND_CELL_COUNTS]
ax.legend(handles, [f"{c} cells" for c in LEGEND_CELL_COUNTS], title="Cell count", loc="center",
          frameon=False, fontsize=8, title_fontsize=9)
ax.axis("off")
fig.tight_layout()
fig.savefig(out_dir / "co_clust_dotplot_legend_sizes.svg", dpi=300, transparent=True, bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(2.0, 3.0))
colorbar = fig.colorbar(mpl.cm.ScalarMappable(cmap="Reds", norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax)),
                        ax=ax, orientation="vertical")
colorbar.ax.tick_params(labelsize=8)
colorbar.set_label("Spearman Correlation", fontsize=8)
ax.axis("off")
fig.tight_layout()
fig.savefig(out_dir / "co_clust_dotplot_legend_colorbar.svg", dpi=300, transparent=True, bbox_inches="tight")
plt.close(fig)

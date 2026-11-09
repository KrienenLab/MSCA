#!/usr/bin/env python
"""Homology clusters between human (Siletti), marmoset and mouse; input of Fig. 2C (expressolog scores).

Marmoset clusters and mouse supertypes were mapped to Siletti clusters with MapMyCells. Those mapped with a
proportion > MAPPING_THRESHOLD to a Siletti cluster that both species map to are averaged (raw metacell
expression) into one 'pseudo' Siletti cluster per species. The three tables (Siletti clusters shared by all
species x highly variable homologous genes) are written for Fig2C_expressolog_human_marm_mouse.R.

Inputs (./data)
    S4 - Cluster Info.xlsx                       marmoset cluster -> Siletti cluster (sheet 'ST4A Cluster Info')
    siletti_on_mouse_abc_anno_table.csv          mouse -> Siletti cluster (MapMyCells)
    siletti_genes.csv                            human Ensembl ID -> gene symbol
    human_marm_mouse_genes.csv                   human Ensembl ID -> gene symbol of the homologous gene set
    hvg_genes.txt                                highly variable genes (one per line)
    marm_Cluster_v4_metacell.h5ad, human_cluster_id_metacell.h5ad, mouse_abc_supertype_metacell.h5ad
                                                 metacells (raw expression in .raw)
Outputs: homology_clusters_{marmoset,human,mouse}_all.csv in ../Fig2C_expressolog_human_marm_mouse/data
"""

from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

try:
    from anndata.io import read_elem
except ImportError:  # anndata < 0.11
    from anndata.experimental import read_elem

script_dir = Path(__file__).parent.resolve()
data_dir = script_dir / "data"
out_dir = script_dir.parent / "Fig2C_expressolog_human_marm_mouse" / "data"
out_dir.mkdir(parents=True, exist_ok=True)

MAPPING_THRESHOLD = 0.7  # fraction of a cluster's cells mapped to a Siletti cluster
MARMOSET_EXCLUDED_CLASSES = ["Nonneuron", "Junk"]
MOUSE_EXCLUDED_CLASSES = ["31 OPC-Oligo", "34 Immune", "30 Astro-Epen"]


# --------------------------------------------------
# Helpers
# --------------------------------------------------
def homologous_siletti_cluster(labels, keys):
    """Siletti cluster ('siletti_382') of every row whose top MapMyCells label has proportion > threshold.

    `labels` look like 'Splat_382(0.85)| Splat_383(0.1)'; the result is indexed by `keys`.
    """
    parsed = labels.str.split("|").str[0].str.extract(r"^\s*([^(]*)\(([0-9.]+)\)")
    keep = parsed[1].astype(float) > MAPPING_THRESHOLD
    siletti = "siletti_" + parsed.loc[keep, 0].str.strip().str.split("_").str[1]
    return pd.Series(siletti.to_numpy(), index=keys[keep.to_numpy()])


def read_raw_expression(path, gene_names=None, drop_unmapped=False):
    """Raw metacell expression (metacells x genes).

    Genes are renamed with the `gene_names` Series (Ensembl ID -> symbol) if given. Unmapped genes are dropped
    (drop_unmapped) or keep their ID; duplicated symbols are then made unique ('GENE', 'GENE-1', ...).
    """
    # obs and raw are read directly: .uns (not needed) can be written by a newer anndata than the one installed
    with h5py.File(path, "r") as f:
        obs_names = read_elem(f["obs"]).index
        raw_x = read_elem(f["raw/X"])
        raw_var_names = read_elem(f["raw/var"]).index
    x = raw_x.toarray() if hasattr(raw_x, "toarray") else np.asarray(raw_x)
    expr = pd.DataFrame(x, index=obs_names, columns=raw_var_names)
    if gene_names is None:
        return expr
    symbols = pd.Series(raw_var_names).map(gene_names)
    if drop_unmapped:
        expr.columns = symbols.to_numpy()
        return expr.loc[:, symbols.notna().to_numpy()]
    expr.columns = ad.utils.make_index_unique(pd.Index(symbols.fillna(pd.Series(raw_var_names))))
    return expr


def mean_by_homolog(expr, homolog):
    """Average the metacells that map to the same Siletti cluster."""
    mapped = homolog.reindex(expr.index).dropna()
    return expr.loc[mapped.index].groupby(mapped).mean()


# --------------------------------------------------
# Mapping of marmoset clusters / mouse supertypes to Siletti clusters
# --------------------------------------------------
marm_anno = pd.read_excel(data_dir / "S4 - Cluster Info.xlsx", sheet_name="ST4A Cluster Info", skiprows=35)
marm_anno = marm_anno[~marm_anno["Subcortex_Class_v4"].isin(MARMOSET_EXCLUDED_CLASSES)]
marm_homolog = homologous_siletti_cluster(marm_anno["Siletti_MapMyCells_Cluster_label"], marm_anno["Cluster_id"])

mouse_anno = pd.read_csv(data_dir / "siletti_on_mouse_abc_anno_table.csv", index_col=0)
mouse_anno = mouse_anno[~mouse_anno["mouse_abc_class"].isin(MOUSE_EXCLUDED_CLASSES)]
mouse_homolog = homologous_siletti_cluster(mouse_anno["Siletti_MapMyCells_Cluster name_name"],
                                           mouse_anno["mouse_abc_supertype"])

# the same Siletti clusters must be used for both species
shared = set(marm_homolog) & set(mouse_homolog)
marm_homolog = marm_homolog[marm_homolog.isin(shared)]
mouse_homolog = mouse_homolog[mouse_homolog.isin(shared)]
# the mouse table has one row per cluster; a supertype takes the Siletti cluster of its last row
mouse_homolog = mouse_homolog[~mouse_homolog.index.duplicated(keep="last")]

# --------------------------------------------------
# Metacells of the three species on the homologous genes
# --------------------------------------------------
siletti_genes = pd.read_csv(data_dir / "siletti_genes.csv", index_col=0)["Gene"]
mouse_genes = pd.read_csv(data_dir / "human_marm_mouse_genes.csv")
mouse_genes["ensembl_id"] = mouse_genes["ensembl_id-0"].str.split(".").str[0]  # drop the version

marm = read_raw_expression(data_dir / "marm_Cluster_v4_metacell.h5ad")
human = read_raw_expression(data_dir / "human_cluster_id_metacell.h5ad", gene_names=siletti_genes)
mouse = read_raw_expression(data_dir / "mouse_abc_supertype_metacell.h5ad",
                            gene_names=mouse_genes.set_index("ensembl_id")["gene"], drop_unmapped=True)

hvg_genes = set((data_dir / "hvg_genes.txt").read_text().split())
genes = sorted(set(marm.columns) & set(human.columns) & set(mouse.columns) & hvg_genes)
for name, expr in (("marmoset", marm), ("human", human), ("mouse", mouse)):
    if expr[genes].columns.duplicated().any():
        raise ValueError(f"Duplicated gene symbols in the {name} metacells")

# --------------------------------------------------
# Pseudo Siletti cluster metacells
# --------------------------------------------------
tables = {
    "marmoset": mean_by_homolog(marm[genes], marm_homolog),
    "mouse": mean_by_homolog(mouse[genes], mouse_homolog),
    "human": human[genes].sort_index(),
}
tables["human"].index = "siletti_" + tables["human"].index

tables["marmoset"] = tables["marmoset"][tables["marmoset"].index.isin(tables["human"].index)]
tables["mouse"] = tables["mouse"][tables["mouse"].index.isin(tables["human"].index)]
tables["human"] = tables["human"][tables["human"].index.isin(tables["marmoset"].index)]

# clusters present in all three species, in the same order
common_clusters = [c for c in tables["mouse"].index if c in tables["marmoset"].index and c in tables["human"].index]
for species, table in tables.items():
    table.loc[common_clusters].to_csv(out_dir / f"homology_clusters_{species}_all.csv")
    print(f"{species}: {len(common_clusters)} clusters x {len(genes)} genes")

"""Shared paths, constants and plotting style for the HPF anterior-posterior gradient analysis.

Edit DATA_DIR / RESULTS_DIR here if your files live somewhere else.
"""
from pathlib import Path

import matplotlib.pyplot as plt

REPO_DIR = Path(__file__).resolve().parent
DATA_DIR = REPO_DIR / "data"        # processed inputs for notebooks 1-2 (download, see README)
RAW_DIR = DATA_DIR / "raw"          # full datasets, only needed to rerun notebook 0
RESULTS_DIR = REPO_DIR / "results"  # tables and figures written by the notebooks
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ---- processed inputs (used by notebooks 1-2; created by notebook 0) -----------------------
SPATIAL_HPF_H5AD = DATA_DIR / "spatial_hpf_data_subset_spatially_variable_genes.h5ad"
RNA_HPF_H5AD = DATA_DIR / "rna_hpf_data_subset.h5ad"
ATAC_PEAKS_H5AD = DATA_DIR / "hpf_atac_cells_by_eregulon_peaks.h5ad"
EREGULON_TSV = DATA_DIR / "eRegulon_direct.tsv"
UMAP1D_DIR = DATA_DIR / "umap1d"    # precomputed 1D UMAP embeddings: "<group>_umap1d.csv"

# ---- raw inputs (used by notebook 0 only) --------------------------------------------------
SPATIAL_RAW_H5AD = RAW_DIR / "Xenium-CJ23.56.004-marmoset-subcortex-cluster-v4_PREPRINT_FREEZE_20251022_subcortex_neurons.h5ad"
RNA_RAW_H5AD = RAW_DIR / "rna_revision_20260226.h5ad"
ANNO_TABLE_CSV = RAW_DIR / "cluster_anno_table.csv"
ALIGNMENT_CSV = RAW_DIR / "1003-marmoset-full.csv"
ATAC_PEAKS_FULL_H5AD = RAW_DIR / "hpf_atac_cells_by_462k_peaks.h5ad"

# ---- shared constants ----------------------------------------------------------------------
GROUP_PALETTE = {
    "CA3 Glut": "#c63631",
    "CA1-ProS Glut": "#284C83",
    "SUB-ProS Glut": "#705aa1",
    "CA2-CA3-IG Glut": "#e17990",
    "DG Glut": "#ADD7A1",
}

PLOT_STYLE = {
    "figure.dpi": 100,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "lines.linewidth": 0.3,
    "font.size": 8,
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
}


def set_plot_style(**overrides):
    """Apply the shared matplotlib style (optionally overriding individual rcParams)."""
    plt.rcParams.update({**PLOT_STYLE, **overrides})

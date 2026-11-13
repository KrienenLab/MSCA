#!/usr/bin/env python
#
# Generate spatial plots to show THM clusters with MB regional peaks,
# which highlight the localization of those clusters.
#
# There is a matrix of z-scores for that peakset in all clusters, but
# in practice I prefer the peakset eigenvalues.
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Patch
from pathlib import Path
import re

pdir         = Path("/jukebox/krienen/lab_members/mike/251007_HMBA_Spatial_ATAC")
path_h5ad    = Path(pdir, "common_data/Cluster_v4_by_462k_peaks_modcopy.h5ad")   ## used for annotations
path_spatial = Path("/jukebox/krienen/lab_members/daisy/spatial_for_Fenna/Xenium_update_annot.h5ad")

# these are the MB regional peaks that were seen in THM GABA clusters
path_mb_peaks_z = Path(pdir, "RegionalPeaksBySpecificity/data/MB_peaks_in_THM_clusters_zscore_mat.csv")  ## [not used]
path_mb_eigen   = Path(pdir, "RegionalPeaksBySpecificity/data/MB_peaks_in_THM_clusters_cl_eigen.csv")
path_mb_peaks   = Path(pdir, "RegionalPeaksBySpecificity/data/MB_peaks_in_THM_clusters.txt")
path_thm_clst   = Path(pdir, "RegionalPeaksBySpecificity/data/THM_clusters_with_MB_peaks.txt")

pdir_plots      = Path(pdir, "RegionalPeaksBySpecificity/plots/plots_spatial")
pdir_plots.mkdir(parents=False, exist_ok=True)


# --------------------------------------------------
# Imports
# --------------------------------------------------

adata       = ad.read_h5ad(path_spatial, backed='r')
anno        = ad.read_h5ad(path_h5ad, backed='r').obs

mb_peaks_z  = pd.read_csv(path_mb_peaks_z, index_col=0)  # peak-by-cluster [not used]
mb_peaks_ei = pd.read_csv(path_mb_eigen, index_col=0)    # vals for each cluster
mb_peaks_ei.rename(columns={"PC1": "MB_peaks_in_THM_eigen"}, inplace=True) 


# --------------------------------------------------
# Note on a cluster missing in ATAC pseudobulk
# --------------------------------------------------

## Not all clusters in spatial in ATAC pseudobulk:
# >>> all(adata.obs['CDM_Cluster_v4_label'].isin(mb_peaks_ei.index))
# False

## Which cluster?
# >>> np.setdiff1d(adata.obs['CDM_Cluster_v4_label'], mb_peaks_ei.index)
# array(['cluster_628'], dtype=object)

# >>> 'cluster_628' in anno.index
# True

## It's class THM Glut, group AD Glut.

## All from a single donor_slab_tile 
##  [this from a different 'anno' though, I dropped that column]
# >>> anno.loc['cluster_628', 'donor_slab_tile']
# 'cjJellybean_Slab6_Tile4(1.0)

## 1540 spatial cells map to it
## 
## Assume this was an ATAC library QC fail
## 

##
## [!] Note, because I'm using the ATAC pseudobulk as the source
##     of the annotations, those cells will just not appear at all
##     in my plots anywhere. I think that's fair because I'm 
##     specifically plotting the ATAC data.
##


# --------------------------------------------------
# Make simplified dataframe for plotting
# --------------------------------------------------

# Won't do a bunch of joins or anything with the anndatas, will just
# extract and collect what I need for plotting

# Want dataframe indexed by cluster
spatial_cols_keep = [
    'z_order',
    'brain_section_label',
    'x',
    'y',
    'CDM_Cluster_v4_label'
]

### Note the 'x','y' in obs are original coordinates, but I think we want the 'rotated_spatial'
# >>> adata.obsm
# AxisArrays with keys: rotated_spatial, spatial_cirro, spatial_cirro_grid, spatial_cirro_grid_4col, spatial_cirro_grid_6col
#
# >>> all(adata.obs['x'] == adata.obsm['spatial_cirro'][:,0])
# True
# >>> all(adata.obs['y'] == adata.obsm['spatial_cirro'][:,1])
# True
#
df_plt = adata.obs.loc[:, spatial_cols_keep]
df_plt.rename(columns={'x':'x_orig', 'y':'y_orig'}, inplace=True)
df_plt['x'] = adata.obsm['rotated_spatial'][:,0]
df_plt['y'] = adata.obsm['rotated_spatial'][:,1]

## flip coordinates! [upside-down & backwards]
df_plt['x'] *= -1
df_plt['y'] *= -1


## drop the 1540 cells mapping to cluster_628 [like 3-in-10k cells]
df_plt = df_plt.loc[df_plt.CDM_Cluster_v4_label.isin(anno.index), :]

## Join in annotations
anno_cols_keep = ['Subcortex_Class_v4', 'Subcortex_Group_v4', 'NT', 'Dominant_Region']
df_plt = df_plt.join(anno.loc[:, anno_cols_keep], on = 'CDM_Cluster_v4_label')

## For plotting
df_plt['Region_NT'] = df_plt['Dominant_Region'].astype(str) + '_' + df_plt['NT'].astype(str)


# --------------------------------------------------
# Add eigenvalues for the peak-set (for each cluster)
# --------------------------------------------------

df_plt = df_plt.join(mb_peaks_ei, on="CDM_Cluster_v4_label")


# --------------------------------------------------
# [skip] Add mean z-scores for peak set to matrix
# --------------------------------------------------

# # mean across peaks for each cluster -> Series indexed by cluster id
# mean_z_by_cluster = mb_peaks_z.mean(axis=0)
# mean_z_by_cluster.name = "MB_peaks_in_THM_mean_z"
# mean_z_df = mean_z_by_cluster.to_frame()

# # Join onto the plotting obs frame
# df_plt = df_plt.join(mean_z_df, on="CDM_Cluster_v4_label")


# --------------------------------------------------
# Drop cells mapping to Junk
# --------------------------------------------------

# >>> sum(df_plt['Subcortex_Class_v4'] == 'Junk')
# 235018
df_plt = df_plt.loc[df_plt['Subcortex_Class_v4'] != 'Junk', ]




# --------------------------------------------------
# z position vs. slice name
# --------------------------------------------------

df_z_slice = adata.obs.loc[:, ['z_order', 'brain_section_label']].drop_duplicates().reset_index().drop(columns="cell_label",axis=1)

# >>> df_z_slice
#    z_order      brain_section_label
# 0        1  CJ23.56.004.CX.42.01.06
# 1        2  CJ23.56.004.CX.42.02.06
# 2        3  CJ23.56.004.CX.42.03.07
# 3        4  CJ23.56.004.CX.42.04.06
# 4        5  CJ23.56.004.CX.42.05.06
# 5        6  CJ23.56.004.CX.42.06.06
# 6        7  CJ23.56.004.CX.42.07.06
# 7        8  CJ23.56.004.CX.42.08.06
# 8        9  CJ23.56.004.CX.42.09.07
# 9       10  CJ23.56.004.CX.42.10.07
# 10      11  CJ23.56.004.CX.42.13.07
# 11      12  CJ23.56.004.CX.42.14.06
# 12      13  CJ23.56.004.CX.42.15.06
# 13      14  CJ23.56.004.CX.42.16.06
# 14      15  CJ23.56.004.CX.42.17.06
# 15      16  CJ23.56.004.CX.42.18.07
# 16      17  CJ23.56.004.CX.42.19.06
# 17      18  CJ23.56.004.CX.42.20.07
# 18      19  CJ23.56.004.CX.42.22.06
# 19      20  CJ23.56.004.CX.43.01.07
# 20      21  CJ23.56.004.CX.43.04.06
# 21      22  CJ23.56.004.CX.43.05.07
# 22      23  CJ23.56.004.CX.43.06.06
# 23      24  CJ23.56.004.CX.43.08.07
# 24      25  CJ23.56.004.CX.43.10.06
# 25      26  CJ23.56.004.CX.43.11.06
# 26      27  CJ23.56.004.CX.43.12.06
# 27      28  CJ23.56.004.CX.43.13.07
# 28      29  CJ23.56.004.CX.43.15.06
# 29      30  CJ23.56.004.CX.43.17.06
# 30      31  CJ23.56.004.CX.43.19.07
# 31      32  CJ23.56.004.CX.43.20.06
# 32      33  CJ23.56.004.CX.43.21.07
# 33      34  CJ23.56.004.CX.44.01.06
# 34      35  CJ23.56.004.CX.44.03.07
# 35      36  CJ23.56.004.CX.44.04.07
# 36      37  CJ23.56.004.CX.44.05.06
# 37      38  CJ23.56.004.CX.44.06.06
# 38      39  CJ23.56.004.CX.44.07.07
# 39      40  CJ23.56.004.CX.44.08.07
# 40      41  CJ23.56.004.CX.44.09.06
# 41      42  CJ23.56.004.CX.44.10.06
# 42      43  CJ23.56.004.CX.44.11.07
# 43      44  CJ23.56.004.CX.44.12.06
# 44      45  CJ23.56.004.CX.44.13.06
# 45      46  CJ23.56.004.CX.44.14.07
# 46      47  CJ23.56.004.CX.44.15.07

z_to_sl = {
    df_z_slice.loc[i, 'z_order'] : df_z_slice.loc[i, 'brain_section_label']
    for i in range(df_z_slice.shape[0])
}


# --------------------------------------------------
# Get slices of interest
# --------------------------------------------------

# adata.obs.brain_section_label
# >>> sorted(list(set(adata.obs.brain_section_label)))
# ['CJ23.56.004.CX.42.01.06', 'CJ23.56.004.CX.42.02.06', 'CJ23.56.004.CX.42.03.07', 'CJ23.56.004.CX.42.04.06', 'CJ23.56.004.CX.42.05.06', 'CJ23.56.004.CX.42.06.06', 'CJ23.56.004.CX.42.07.06', 'CJ23.56.004.CX.42.08.06', 'CJ23.56.004.CX.42.09.07', 'CJ23.56.004.CX.42.10.07', 'CJ23.56.004.CX.42.13.07', 'CJ23.56.004.CX.42.14.06', 'CJ23.56.004.CX.42.15.06', 'CJ23.56.004.CX.42.16.06', 'CJ23.56.004.CX.42.17.06', 'CJ23.56.004.CX.42.18.07', 'CJ23.56.004.CX.42.19.06', 'CJ23.56.004.CX.42.20.07', 'CJ23.56.004.CX.42.22.06', 'CJ23.56.004.CX.43.01.07', 'CJ23.56.004.CX.43.04.06', 'CJ23.56.004.CX.43.05.07', 'CJ23.56.004.CX.43.06.06', 'CJ23.56.004.CX.43.08.07', 'CJ23.56.004.CX.43.10.06', 'CJ23.56.004.CX.43.11.06', 'CJ23.56.004.CX.43.12.06', 'CJ23.56.004.CX.43.13.07', 'CJ23.56.004.CX.43.15.06', 'CJ23.56.004.CX.43.17.06', 'CJ23.56.004.CX.43.19.07', 'CJ23.56.004.CX.43.20.06', 'CJ23.56.004.CX.43.21.07', 'CJ23.56.004.CX.44.01.06', 'CJ23.56.004.CX.44.03.07', 'CJ23.56.004.CX.44.04.07', 'CJ23.56.004.CX.44.05.06', 'CJ23.56.004.CX.44.06.06', 'CJ23.56.004.CX.44.07.07', 'CJ23.56.004.CX.44.08.07', 'CJ23.56.004.CX.44.09.06', 'CJ23.56.004.CX.44.10.06', 'CJ23.56.004.CX.44.11.07', 'CJ23.56.004.CX.44.12.06', 'CJ23.56.004.CX.44.13.06', 'CJ23.56.004.CX.44.14.07', 'CJ23.56.004.CX.44.15.07']

# taking a big range of 19 slices over which to plot
# focus_slices = [
#     'CJ23.56.004.CX.43.10.06', 'CJ23.56.004.CX.43.11.06', 'CJ23.56.004.CX.43.12.06', 'CJ23.56.004.CX.43.13.07', 'CJ23.56.004.CX.43.15.06', 
#     'CJ23.56.004.CX.43.17.06', 'CJ23.56.004.CX.43.19.07', 'CJ23.56.004.CX.43.20.06', 'CJ23.56.004.CX.43.21.07', 'CJ23.56.004.CX.44.01.06', 
#     'CJ23.56.004.CX.44.03.07', 'CJ23.56.004.CX.44.04.07', 'CJ23.56.004.CX.44.05.06', 'CJ23.56.004.CX.44.06.06', 'CJ23.56.004.CX.44.07.07', 
#     'CJ23.56.004.CX.44.08.07', 'CJ23.56.004.CX.44.09.06', 'CJ23.56.004.CX.44.10.06', 'CJ23.56.004.CX.44.11.07'
# ]

## just do all
# focus_slices = sorted(list(adata.obs.brain_section_label.cat.categories))

## and use z number for this for now on
focus_slices = z_to_sl.keys()


# --------------------------------------------------
# Plot eigenvalues for the peakset
# --------------------------------------------------

# plot folders:
path_ei_by_slice = Path(pdir_plots, "mb_peaks_eigenvalue_by_slice")
path_ei_by_slice.mkdir(exist_ok=True)

# for each slice ('brain_section_label') in focus_slices, plot the eigenvalue for the peak set (so x,y color=eigen);
# should be in that folder, file name should be slice name
# robust_low, robust_high = (
#     mb_peaks_ei["MB_peaks_in_THM_eigen"].quantile(0.01),  # bottom 18 clusters
#     mb_peaks_ei["MB_peaks_in_THM_eigen"].quantile(0.995), # like top 9 or 10 clusters
# )

## vs just min and max, pretty similar
low, high = (
    mb_peaks_ei["MB_peaks_in_THM_eigen"].min(),
    mb_peaks_ei["MB_peaks_in_THM_eigen"].max()
)

plt.rcParams.update({
    "figure.dpi": 400,
    "savefig.dpi": 400,
    "axes.facecolor": "white"
})


for z in focus_slices:
    sl = z_to_sl[z]
    d = df_plt[df_plt['z_order'] == z]
    fig, ax = plt.subplots(figsize=(6, 6))
    # pull to numpy once (faster than repeated .loc)
    x = d["x"].to_numpy()
    y = d["y"].to_numpy()
    m = d["MB_peaks_in_THM_eigen"].to_numpy()
    isnan = np.isnan(m)
    # 1) Draw NaNs first -> stays underneath
    if isnan.any():
        ax.scatter(
            x[isnan], y[isnan],
            s=0.6, lw=0,
            c="lightgray", alpha=0.1,
            rasterized=True,
        )
    # 2) Sort valid points by value so higher scores end up on top
    valid = ~isnan
    vals = m[valid]
    order = np.argsort(vals)  # ascending: low first, high last (on top)
    sc = ax.scatter(
        x[valid][order],
        y[valid][order],
        c=vals[order],
        s=1.2,              # ~20% larger than before
        lw=0,
        cmap="viridis",
        vmin=low, vmax=high,
        rasterized=True,
    )
    # cosmetics / layout
    ax.set_title(f"z = {z} / 47\n{sl}\nMB peaks in THM GABA\n(eignvalue of peak set for cluster)", fontsize=10)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.margins(x=0.05, y=0.05)
    # flip flip
    ax.invert_yaxis()
    ax.invert_xaxis()
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("eigenvalue")
    out = path_ei_by_slice / f"z{z}_{sl}.png"
    fig.savefig(out, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)



# --------------------------------------------------
# Plot categorical:
# --------------------------------------------------

# 
# Could use both a little refactoring (to make less "class" specific
# (particularly the synonym remapping needs to go) and might better
# be spun out into a utils script for palettes
#


# Existing palette for class
palette_class_supplied_hex = {
    "BG":               "#CBC3E3",
    "CNU-HTH GABA":     "#E69F00",
    "CNU-HTH Glut":     "#56B4E9",
    "CTX GABA":         "#009E73",
    "HTH GABA":         "#F0E442",
    "HTH Glut":         "#0072B2",
    "HTH Hist":         "#D55E00",
    "MB-HB GABA":       "#CC79A7",
    "MB-HB Glut":       "#44AA99",
    "MB-HTH-HB Sero":   "#AA4499",
    "Nonneuron":        "#999999",
    "Pineal Gland":     "#999999",
    "SEP GABA":         "#882255",
    "TEL Glut":         "#88CCEE",
    "THM GABA":         "#DDCC77",
    "THM Glut":         "#332288",
}

# Map some labels in your data to palette keys
class_synonyms = {
    "CTX-MGE GABA": "CTX GABA",
    "CTX-CGE GABA": "CTX GABA",
    # "Junk":         "Nonneuron",
}

def _hex_to_rgba(hex_str):
    hex_str = hex_str.lstrip("#")
    r = int(hex_str[0:2], 16) / 255.0
    g = int(hex_str[2:4], 16) / 255.0
    b = int(hex_str[4:6], 16) / 255.0
    return (r, g, b, 1.0)

def _rgba_to_hex(rgba):
    r,g,b, *_ = rgba
    return "#{:02X}{:02X}{:02X}".format(int(round(r*255)), int(round(g*255)), int(round(b*255)))

def _discrete_color_pool(names=("tab20", "tab20b", "tab20c", "Set3")):
    pool = []
    for name in names:
        cmap = cm.get_cmap(name)
        pool.extend([cmap(i) for i in range(cmap.N)])
    return pool  # list of RGBA tuples

# not used but i like
def make_palette_hex(labels):
    # preserve input order, drop dups
    labels = list(dict.fromkeys(labels))
    pool = _discrete_color_pool()
    # repeat colors if you have more labels than the pool size
    if len(labels) > len(pool):
        from itertools import islice, cycle
        colors = list(islice(cycle(pool), len(labels)))
    else:
        colors = pool[:len(labels)]
    return {lab: _rgba_to_hex(col) for lab, col in zip(labels, colors)}


def _build_class_palette(df, col="Subcortex_Class_v4"):
    cats = pd.Series(df[col].dropna().unique()).tolist()
    # Apply synonyms
    mapped_keys = [class_synonyms.get(c, c) for c in cats]
    # Start with supplied hex colors (converted to RGBA)
    pal = {k: _hex_to_rgba(v) for k, v in palette_class_supplied_hex.items()}
    # Determine which categories are missing from the supplied palette
    missing = []
    for c, k in zip(cats, mapped_keys):
        if k not in pal:
            missing.append(c)   # note: use the original label for display
    # Assign fallback colors only for the missing labels (stable order)
    if missing:
        pool = _discrete_color_pool()
        # Avoid using colors equal to supplied ones (best-effort)
        used = set(pal.values())
        fallback = []
        for color in pool:
            if color not in used:
                fallback.append(color)
            if len(fallback) >= len(missing):
                break
        pal.update({lab: col for lab, col in zip(missing, fallback)})
    # Finally, make sure every encountered label has a color:
    # - use synonym color if defined
    # - else use its own assigned fallback (above)
    final = {}
    for c in cats:
        key = class_synonyms.get(c, c)
        final[c] = pal[key] if key in pal else pal[c]
    return final

# Global, stable palettes
palette_class_rgba = _build_class_palette(df_plt, "Subcortex_Class_v4")

def _plot_categorical_by_slice(col, outdir, title_label, palette_rgba, legend_cap=20):
    for z in focus_slices:
        sl = z_to_sl[z]
        d = df_plt[df_plt['z_order'] == z]
        if d.empty:
            continue
        x = d["x"].to_numpy()
        y = d["y"].to_numpy()
        cats = d[col].to_numpy()
        # Colors per point
        isnan = pd.isna(cats)
        colors = np.empty((len(cats), 4), dtype=float)
        colors[isnan] = (0.8, 0.8, 0.8, 0.10)
        valid = ~isnan
        valid_cats = cats[valid]
        colors_valid = [palette_rgba.get(c, (0.5, 0.5, 0.5, 1.0)) for c in valid_cats]
        colors[valid] = np.array(colors_valid)
        # Draw rare categories last
        # size_counts = pd.Series(valid_cats).value_counts()
        # priority = np.empty(len(cats), dtype=np.int64)
        # priority[isnan] = -1
        # priority[valid] = [size_counts.get(c, 0) for c in valid_cats]
        # order = np.argsort(priority)  # NaN first, then rare → common
        fig, ax = plt.subplots(figsize=(6, 6))
        # ax.scatter(x[order], y[order], c=colors[order], s=1.2, lw=0, rasterized=True)
        ax.scatter(x, y, c=colors, s=1.2, lw=0, rasterized=True)
        ax.set_title(f"z = {z} / 47\n{sl}\n{title_label}", fontsize=10)
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.margins(x=0.05, y=0.05)
        ax.invert_yaxis()
        ax.invert_xaxis()
        # Legend only for categories present in THIS slice (and cap length)
        present = pd.Series(valid_cats).dropna().unique().tolist()
        present = sorted(present)
        if len(present) <= legend_cap:
            handles = [Patch(facecolor=palette_rgba.get(c, (0.5,0.5,0.5,1.0)),
                             edgecolor="none", label=str(c)) for c in present]
            ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5),
                      frameon=False, fontsize=7, title=col)
        out = outdir / f"z{z}_{sl}.png"
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)


# --------------------------------------------------
# Plot Class
# --------------------------------------------------

path_class_by_slice = Path(pdir_plots, "Subcortex_Class_v4_by_slice")
path_class_by_slice.mkdir(exist_ok=True)

_plot_categorical_by_slice("Subcortex_Class_v4", path_class_by_slice,
                           "Subcortex_Class_v4", palette_class_rgba)


# --------------------------------------------------
# Plot Region
# --------------------------------------------------

path_region_by_slice = Path(pdir_plots, "Dominant_Region_by_slice")
path_region_by_slice.mkdir(exist_ok=True)

palette_region_auto = _build_class_palette(df_plt, "Dominant_Region")
_plot_categorical_by_slice("Dominant_Region", path_region_by_slice,
                           "Dominant_Region", palette_region_auto)



# --------------------------------------------------
# Plot NT
# --------------------------------------------------

path_nt_by_slice = Path(pdir_plots, "NT_by_slice")
path_nt_by_slice.mkdir(exist_ok=True)

palette_nt_auto = _build_class_palette(df_plt, "NT")
_plot_categorical_by_slice("NT", path_nt_by_slice,
                           "NT", palette_nt_auto)


# --------------------------------------------------
# Eigen by NT (per-slice): show only one NT colored, others gray
# Files like: {slice}_{NT}.png
# --------------------------------------------------

path_ei_by_slice_nt = Path(pdir_plots, "mb_peaks_eigenvalue_by_slice_and_NT")
path_ei_by_slice_nt.mkdir(exist_ok=True)

ei_low, ei_high = (
    mb_peaks_ei["MB_peaks_in_THM_eigen"].min(),
    mb_peaks_ei["MB_peaks_in_THM_eigen"].max(),
)

def _slugify(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", str(s))

# Unique NT values present in the data (skip NaN)
all_nts = sorted(pd.Series(df_plt["NT"]).dropna().unique().tolist())

for z in focus_slices:
    sl = z_to_sl[z]
    d = df_plt[df_plt['z_order'] == z]
    if d.empty:
        continue
    # NTs actually present in this slice
    nts_in_slice = sorted(pd.Series(d["NT"]).dropna().unique().tolist())
    if not nts_in_slice:
        continue
    # Pull arrays once per slice
    x_all = d["x"].to_numpy()
    y_all = d["y"].to_numpy()
    eig_all = d["MB_peaks_in_THM_eigen"].to_numpy()
    nt_all = d["NT"].to_numpy()
    for nt in nts_in_slice:
        target = (nt_all == nt)
        fig, ax = plt.subplots(figsize=(6, 6))
        # 1) Background: all non-target cells, plotted first in faint gray
        bg = ~target
        if np.any(bg):
            ax.scatter(
                x_all[bg], y_all[bg],
                s=0.6, lw=0,
                c="lightgray", alpha=0.10,
                rasterized=True,
            )
        # 2) Target cells: color by eigen (NaNs as faint gray), highs on top
        m_t = eig_all[target]
        x_t = x_all[target]
        y_t = y_all[target]
        if m_t.size:
            isnan_t = np.isnan(m_t)
            # NaN target points as faint gray underlay
            if np.any(isnan_t):
                ax.scatter(
                    x_t[isnan_t], y_t[isnan_t],
                    s=0.6, lw=0,
                    c="lightgray", alpha=0.10,
                    rasterized=True,
                )
            # Valid target points, sorted so high eigen draws last (on top)
            valid_t = ~isnan_t
            if np.any(valid_t):
                vals = m_t[valid_t]
                order = np.argsort(vals)  # low first, high last
                sc = ax.scatter(
                    x_t[valid_t][order],
                    y_t[valid_t][order],
                    c=vals[order],
                    s=1.2,
                    lw=0,
                    cmap="viridis",
                    vmin=ei_low, vmax=ei_high,
                    rasterized=True,
                )
                # Colorbar only if we had valid (non-NaN) eigen values
                cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
                cbar.set_label("eigenvalue")
        # Cosmetics / layout
        ax.set_title(f"z = {z} / 47\n{sl}\nPeak-set eigenvalue — NT: {nt}", fontsize=10)
        ax.set_aspect("equal")
        ax.set_axis_off()
        ax.margins(x=0.05, y=0.05)
        ax.invert_yaxis()
        ax.invert_xaxis()
        out = path_ei_by_slice_nt / f"z{z}_{_slugify(sl)}_{_slugify(nt)}.png"
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)


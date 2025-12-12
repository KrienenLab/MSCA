#!/usr/bin/env python3
"""
Make bigWig (actually bedGraph here...) files for any grouping variable(s).

--> Something odd is happening where the bigWig writer quits unexpectedly (gets 'disconnected')
    after a while (sometimes a day). As a workaround, this script will export bedgraph files
    instead, and then I'll use another script to convert the bedgraphs to bigwigs
    
    [*] Also the combination of writing bedGraphs and then using bedGraphToBigWig is faster
"""
import anndata as ad
from pathlib import Path
import os
import shutil
import snapatac2 as snap
import pandas as pd
import polars as pl

from functools import partial
print = partial(print, flush=True)

# ================== CONFIG ====================
# if multiple grouping variables, will make a new column that concatenates them by 'collapse';
# grouping_vars  = ["Subcortex_Class_v4"]
grouping_vars = ["Subcortex_Group_v4"]
# grouping_vars  = ['Cluster_Dominant_Region', 'Cluster_NT']
collapse       = "__"
assert isinstance(grouping_vars, list)
outname        = collapse.join(grouping_vars)
print(f"Running make bigWigs for column: {outname}")

pdir           = '/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster'
path_atac_h5ad = Path(pdir, 'data/atac/atac_h5ads/_dataset.h5ads')
# path_anno      = Path(pdir, 'metadata/subcortex_anno_table.csv')
# path_anno      = Path(pdir, 'metadata/subcortex_anno_table_archive/subcortex_anno_table_20251008.csv')
path_anno      = Path(pdir, 'metadata/preprint_freeze/subcortex_anno_table_preprint.csv')
# path_reg_anno  = Path(pdir, 'metadata/subcortex_anno_table_archive/subcortex_anno_table_cluster_region_20251008.csv')  # <--- TEMPORARY
out_dir        = Path(pdir, f"data/atac/bw/{outname}/cpm_bedgraph") # <---- WHILE USING BEDGRAPH
out_dir_raw    = Path(pdir, f"data/atac/bw/{outname}/raw_bedgraph") # <---- WHILE USING BEDGRAPH
out_dir_cov    = Path(pdir, f"data/atac/bw/{outname}/coverage_rpkm_bedgraph") # <---- WHILE USING BEDGRAPH
out_dir_cov_raw= Path(pdir, f"data/atac/bw/{outname}/coverage_raw_bedgraph")  # <---- WHILE USING BEDGRAPH
tempdir        = Path(pdir, ".make_bw_tmp")
# ==============================================

# Rust diagnostics
os.environ.setdefault("RUST_BACKTRACE", "1")
os.environ.setdefault("RUST_LOG", "bigtools=debug,snapatac2=info")

print("=== Paths ===")
print("pdir            :", pdir)
print("h5ad            :", path_atac_h5ad)
print("anno            :", path_anno)
# print("reg         :", path_reg_anno)
print("out_dir.        :", out_dir)
print("out_dir_raw     :", out_dir_raw)
print("out_dir_cov     :", out_dir_cov)
print("out_dir_cov_raw :", out_dir_cov_raw)
print("tempdir         :", tempdir)

os.makedirs(out_dir, exist_ok=True)
os.makedirs(out_dir_raw, exist_ok=True)
os.makedirs(out_dir_cov, exist_ok=True)
os.makedirs(out_dir_cov_raw, exist_ok=True)
os.makedirs(tempdir, exist_ok=True)

print("\n=== Load dataset (r+) ===")
adata = snap.read_dataset(path_atac_h5ad, mode='r+')

print("\n=== Load current annotations ===")
anno = pd.read_csv(path_anno).rename(columns={'marm_cluster_id': 'Cluster_v4'})
anno.set_index('Cluster_v4', inplace=True)

# this script will modify the `_dataset.h5ads` file, but will at least not permanently save
# the added metadata columns into it.
# (!) NOTE: Existing columns that get updated will maintain those updates!!
cols_added = []

def fix_pandas_types_NAs(col, nm):
    orig_dtype = col.dtype
    dtype_str = str(orig_dtype)
    # Attempt to coerce object columns into atomic types
    if pd.api.types.is_object_dtype(orig_dtype):
        try:
            col_converted = pd.to_numeric(col, errors='raise')
            col = col_converted.astype(float)  # preserves NaNs
        except Exception:
            try:
                col_lower = col.astype(str).str.lower()
                if col_lower.isin(['true', 'false', 'nan', 'none']).all():
                    col = col_lower.map({'true': True, 'false': False}).astype('boolean')
                else:
                    raise ValueError
            except Exception:
                col = col.astype(str).replace({'nan': pd.NA, 'None': pd.NA})
    # Ensure type is Polars/Arrow-safe
    if pd.api.types.is_numeric_dtype(col):
        col = col.astype(float)  # ensure float for missing
    elif pd.api.types.is_bool_dtype(col):
        col = col.astype('boolean')  # Arrow-compatible
    elif pd.api.types.is_string_dtype(col):
        col = col.replace({'nan': pd.NA, 'None': pd.NA})
    elif pd.api.types.is_categorical_dtype(col):
        col = col.astype(str).replace({'nan': pd.NA})
    else:
        print(f"[WARNING] grouping_var '{nm}' has unsupported dtype ({dtype_str}); converting to string")
        col = col.astype(str).replace({'nan': pd.NA})
    return col

print("\n=== Sync grouping vars into adata.obs ===")
for gv in grouping_vars:
    if gv in anno.columns:
        # assume the grouping_var is applied/consistent at the cluster-level
        new_entries = anno.loc[adata.obs['Cluster_v4'], gv]
        new_entries = fix_pandas_types_NAs(new_entries, gv)
        new_entries = pl.Series(new_entries)
        if not gv in adata.obs:
            print(
                f"  ~ Found `{gv}` in annotation file, but not in adata.obs. Will update the file now but remove it later in script.",
                flush=True
            )
            adata.obs[gv] = new_entries
            cols_added += [gv]
        elif not all(new_entries == adata.obs[gv]):
            print(f"  ~ Found `{gv}` in adata.obs, but different values are present in the annotation table. " \
            "UPDATING THE `_dataset.h5ads` FILE!!")
            adata.obs[gv] = new_entries
        else:
            print(f"  ~ Found `{gv}` in adata.obs and annotation table, and both identical")
    elif gv in adata.obs:
        print(f"  ~ Found `{gv}` in adata.obs (but not in annotation table)")
    else:
        assert False, f"[!] Didn't find `{gv}` in adata or in annotation table."

# if multiple, make a new concatenated column
if len(grouping_vars) > 1:
    print("\n=== Build concatenated column ===")
    out = adata.obs[grouping_vars[0]]
    for name in grouping_vars[1:]:
        out = out + collapse + adata.obs[name]
    adata.obs[outname] = out
    cols_added += [outname]


# also be sure no NA values that will fail when making bigWigs
# 
# ----> [!] maybe this should always be done, but would result in always updating the input file
# 
print("\n=== Sanitize group column for filenames ===")
if outname in cols_added:
    adata.obs[outname] = (
        adata.obs[outname]
        .cast(pl.Utf8)
        .fill_null("NA")
        .str.replace_all(r"[^\w\-.]+", "_")  # make filenames safe
    )


print("\n=== Make fragment files ===")
print("Making CPM bedGraphs (in lieu of bigWigs) -> need to be converted...")
snap.ex.export_coverage(
    adata, 
    groupby           = outname, 
    normalization     = "CPM", 
    counting_strategy = "insertion", 
    bin_size          = 1,
    suffix            = '.bedgraph',   # filetype inferred from suffix
    out_dir           = out_dir, 
    n_jobs            = 1,
    tempdir           = tempdir
)

print("Making raw bedGraphs (in lieu of bigWigs) -> need to be converted...")
snap.ex.export_coverage(
    adata, 
    groupby           = outname, 
    normalization     = None, 
    counting_strategy = "insertion", 
    bin_size          = 1,
    suffix            = '.bedgraph',   # filetype inferred from suffix
    out_dir           = out_dir_raw, 
    n_jobs            = 1,
    tempdir           = tempdir
)

# print("Making RAW WHOLE READ COVERAGE bedGraphs (in lieu of bigWigs) -> need to be converted...")
# snap.ex.export_coverage(
#     adata, 
#     groupby           = outname, 
#     normalization     = None, 
#     counting_strategy = "fragment", 
#     bin_size          = 1,
#     suffix            = '.bedgraph',   # filetype inferred from suffix
#     out_dir           = out_dir_cov_raw, 
#     n_jobs            = 1,
#     tempdir           = tempdir
# )

# print("Making RPKM WHOLE READ COVERAGE bedGraphs (in lieu of bigWigs) -> need to be converted...")
# snap.ex.export_coverage(
#     adata, 
#     groupby           = outname, 
#     normalization     = "RPKM", 
#     counting_strategy = "fragment", 
#     bin_size          = 1,
#     suffix            = '.bedgraph',   # filetype inferred from suffix
#     out_dir           = out_dir_cov, 
#     n_jobs            = 1,
#     tempdir           = tempdir
# )


if cols_added:
    print("\n=== Cleanup ===")
    print(f'Removing added columns from `_dataset.h5ads`: {cols_added}')
    for col in cols_added:
        try:
            del adata.obs[col]
        except Exception as e:
            print(f"Warning: couldn't remove {col}: {e}")

adata.close()
shutil.rmtree(tempdir)
print("\nScript finished.")

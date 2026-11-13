#!/usr/bin/env python3
"""
Written and documented by an LLM but I read through it reasonably carefully to verify.

Pseudobulk Xenium spatial counts by Cluster_v4_name (MMC Cluster_v4), summing panel-gene reads
across qc_pass cells.

Outputs:
  1) genes x clusters CSV: rows = genes (var_names), columns = Cluster_v4_name, values = summed counts
  2) cluster stats CSV: one row per cluster with qc_pass cell counts (and some handy totals)

Notes:
- Uses a backed .h5ad (does not load the full matrix).
- Aggregates ONLY cells with qc_pass==True AND a non-null cluster label.
- By default uses spatial.layers["raw"].

Compression:
- pandas `to_csv("file.csv.gz")` wILL gzip-compress automatically based on the filename suffix.
"""

import argparse
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def _as_dense_int64(x):
    """
    Convert a small (chunk) matrix x to a dense int64 numpy array.
    Works for scipy sparse, numpy arrays, and many array-likes.
    """
    if sp.issparse(x):
        arr = x.toarray()
    else:
        arr = np.asarray(x)

    # raw counts may be stored as float32; round defensively then cast
    if not np.issubdtype(arr.dtype, np.integer):
        arr = np.rint(arr)

    return arr.astype(np.int64, copy=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True, help="Input Xenium spatial .h5ad")
    ap.add_argument("--out-counts", required=True,
                    help="Output genes×clusters counts CSV (optionally .gz)")
    ap.add_argument("--out-cluster-stats", required=True,
                    help="Output per-cluster stats CSV (optionally .gz)")
    ap.add_argument("--layer", default="raw", help="Layer to use (default: raw)")
    ap.add_argument("--cluster-col", default="Cluster_v4_name", help="obs column for clusters")
    ap.add_argument("--qc-col", default="qc_pass", help="obs QC pass boolean column")
    ap.add_argument("--chunk-size", type=int, default=50000, help="Cells per chunk (default: 50k)")
    ap.add_argument("--drop-empty-clusters", action="store_true",
                    help="Drop clusters with 0 qc_pass assigned cells (default: keep all)")
    args = ap.parse_args()

    # Load backed
    spatial = ad.read_h5ad(args.h5ad, backed="r")

    n_obs, n_vars = spatial.n_obs, spatial.n_vars
    genes = spatial.var_names.astype(str).to_list()

    obs = spatial.obs
    if args.qc_col not in obs.columns:
        raise KeyError(f"qc col '{args.qc_col}' not found in spatial.obs")
    if args.cluster_col not in obs.columns:
        raise KeyError(f"cluster col '{args.cluster_col}' not found in spatial.obs")

    qc = obs[args.qc_col].to_numpy(dtype=bool)

    # Use categorical for stable cluster ordering + integer codes
    cl_cat = pd.Categorical(obs[args.cluster_col])
    cluster_names = cl_cat.categories.astype(str).to_list()
    cluster_codes = cl_cat.codes.astype(np.int32)  # -1 means NaN/unassigned
    n_clusters = len(cluster_names)

    # Valid rows for aggregation
    valid = qc & (cluster_codes >= 0)

    # Access layer lazily per chunk
    if args.layer not in spatial.layers.keys():
        raise KeyError(
            f"layer '{args.layer}' not found in spatial.layers "
            f"(available: {list(spatial.layers.keys())})"
        )
    layer = spatial.layers[args.layer]

    chunk = args.chunk_size
    if chunk <= 0:
        raise ValueError("--chunk-size must be > 0")

    # Progress summary
    n_qc = int(qc.sum())
    n_valid = int(valid.sum())
    n_unassigned_qc = int((qc & (cluster_codes < 0)).sum())
    print(f"[info] n_obs={n_obs:,}  n_vars={n_vars:,}  n_clusters={n_clusters:,}", file=sys.stderr)
    print(f"[info] qc_pass cells: {n_qc:,}", file=sys.stderr)
    print(f"[info] qc_pass & assigned cells: {n_valid:,}", file=sys.stderr)
    print(f"[info] qc_pass but unassigned cluster: {n_unassigned_qc:,}", file=sys.stderr)
    print(f"[info] chunk_size={chunk:,}  layer={args.layer}  cluster_col={args.cluster_col}", file=sys.stderr)

    # Accumulators:
    # - sums: clusters × genes
    # - cell_counts: qc_pass assigned cells per cluster
    # - total_counts_per_cluster: summed counts across all genes per cluster (qc_pass assigned only)
    sums = np.zeros((n_clusters, n_vars), dtype=np.int64)
    cell_counts = np.zeros(n_clusters, dtype=np.int64)
    total_counts_per_cluster = np.zeros(n_clusters, dtype=np.int64)

    # Precompute counts per cluster (qc_pass & assigned), cheap and exact
    if valid.any():
        cell_counts += np.bincount(cluster_codes[valid], minlength=n_clusters).astype(np.int64)

    # Iterate chunks and aggregate counts
    for start in range(0, n_obs, chunk):
        stop = min(n_obs, start + chunk)
        vmask = valid[start:stop]
        if not vmask.any():
            continue

        codes = cluster_codes[start:stop][vmask]
        n_pass = int(vmask.sum())

        # Slice the matrix block and then subset to qc_pass rows
        block = layer[start:stop, :]
        try:
            block = block[vmask, :]
        except Exception:
            block = np.asarray(block)[vmask, :]

        X = _as_dense_int64(block)  # (n_pass, n_genes)

        # Build a "one nonzero per row" CSR matrix mapping rows -> cluster code
        data = np.ones(n_pass, dtype=np.int8)
        indptr = np.arange(n_pass + 1, dtype=np.int64)
        G = sp.csr_matrix((data, codes, indptr), shape=(n_pass, n_clusters), dtype=np.int8)

        # Aggregate: (n_clusters, n_genes) += G.T @ X
        chunk_sums = (G.T @ X)  # dense (n_clusters, n_genes) because X is dense
        sums += chunk_sums

        # Also accumulate per-cluster total counts (across genes)
        total_counts_per_cluster += chunk_sums.sum(axis=1).astype(np.int64)

        if (start // chunk) % 20 == 0:
            print(f"[progress] processed {stop:,}/{n_obs:,} cells ...", file=sys.stderr)

    # Optionally drop empty clusters
    if args.drop_empty_clusters:
        keep = cell_counts > 0
        if not keep.all():
            sums = sums[keep, :]
            total_counts_per_cluster = total_counts_per_cluster[keep]
            cell_counts = cell_counts[keep]
            cluster_names = [nm for nm, k in zip(cluster_names, keep) if k]
            n_clusters = len(cluster_names)
            print(f"[info] dropped empty clusters; keeping {n_clusters:,}", file=sys.stderr)

    # Build outputs
    out_counts = pd.DataFrame(sums.T, index=genes, columns=cluster_names)

    out_stats = pd.DataFrame({
        "Cluster_v4_name": cluster_names,
        "qc_pass_assigned_cells": cell_counts,
        "qc_pass_assigned_total_counts": total_counts_per_cluster,
    })
    # Optional convenience metric (avoid divide-by-zero)
    out_stats["qc_pass_assigned_mean_counts_per_cell"] = np.where(
        out_stats["qc_pass_assigned_cells"].to_numpy() > 0,
        out_stats["qc_pass_assigned_total_counts"].to_numpy() / out_stats["qc_pass_assigned_cells"].to_numpy(),
        np.nan,
    )

    # Write
    for pth in [args.out_counts, args.out_cluster_stats]:
        out_dir = os.path.dirname(os.path.abspath(pth))
        if out_dir and not os.path.isdir(out_dir):
            os.makedirs(out_dir, exist_ok=True)

    out_counts.to_csv(args.out_counts, index=True)
    out_stats.to_csv(args.out_cluster_stats, index=False)

    # Sanity prints
    print(f"[done] wrote counts:        {args.out_counts}", file=sys.stderr)
    print(f"[done] wrote cluster stats: {args.out_cluster_stats}", file=sys.stderr)
    print(f"[sanity] aggregated qc_pass assigned cells: {int(cell_counts.sum()):,}", file=sys.stderr)
    print(f"[sanity] summed counts (all genes, all clusters): {int(out_counts.to_numpy().sum()):,}", file=sys.stderr)


if __name__ == "__main__":
    main()

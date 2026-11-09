#!/usr/bin/env python
"""Summarize cell-level metadata per cluster (or group): the annotation tables behind S4 / S5.

For every categorical column the labels of a cluster are written as 'label(fraction)| label(fraction)', most
frequent first, fraction = share of the cluster's cells, labels below MIN_PERCENT omitted. For MapMyCells columns
this gives e.g. 'Splat_382(0.85)| Splat_383(0.1)', the format read by 3a / 3b. Numeric columns are summarized by
their median, and the number of cells of each cluster is added.

Usage
    python 2_make_annotation_table.py --metadata cell_metadata.csv --group-by Cluster_v4 \
        --output cluster_anno_table.csv \
        --categorical donor_name Siletti_MapMyCells_Cluster_label Siletti_MapMyCells_Supercluster_label \
        --numeric nCount_RNA nFeature_RNA percent.mt scDblFinder.score

    --metadata     cell x metadata table (csv, first column = cell id); one row per cell
    --group-by     column defining the clusters (or groups) to summarize
    --categorical  categorical columns to summarize (values are converted to strings)
    --numeric      numeric columns to summarize (median)
    --min-percent  smallest fraction of a cluster that is listed for a label (default 0.05)
    --doublet-flag optional: flag clusters whose median of this numeric column exceeds --doublet-cutoff
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def build_annotation_table(metadata, group_by, categorical_annotations, numeric_annotations, min_percent=0.05):
    """One row per value of `group_by`; see the module docstring for the columns."""
    summary = {}
    clusters = pd.unique(metadata[group_by].dropna())

    for annotation in categorical_annotations:
        print("Processing annotation:", annotation)
        counts = metadata.groupby([group_by, annotation])[annotation].size()
        summary[annotation] = {}
        for cluster in clusters:
            cluster_counts = counts.loc[counts.index.get_level_values(group_by) == cluster]
            cluster_counts = cluster_counts.sort_values(ascending=False)
            fractions = np.round(cluster_counts / cluster_counts.sum(), 2)
            labels = fractions.index.get_level_values(annotation)
            terms = [f"{label}({fraction})" for label, fraction in zip(labels, fractions) if fraction > min_percent]
            summary[annotation][cluster] = "| ".join(terms)

    summary[f"{group_by}_size"] = metadata.groupby(group_by).size()
    for column in numeric_annotations:
        print("Processing annotation:", column)
        summary[f"{column}_median"] = metadata.groupby(group_by)[column].median().round(4)

    return pd.DataFrame.from_dict(summary)


def sort_by_trailing_number(table):
    """Order ids such as 'cluster_2', 'cluster_10' numerically; leave other ids unchanged."""
    numbers = table.index.to_series().astype(str).str.extract(r"(\d+)$")[0]
    if numbers.isna().any():
        return table
    return table.iloc[np.argsort(numbers.astype(int).to_numpy(), kind="stable")]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--group-by", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--categorical", nargs="*", default=[])
    parser.add_argument("--numeric", nargs="*", default=[])
    parser.add_argument("--min-percent", type=float, default=0.05)
    parser.add_argument("--doublet-flag", default=None, help="numeric column used to flag clusters")
    parser.add_argument("--doublet-cutoff", type=float, default=0.1)
    args = parser.parse_args()

    metadata = pd.read_csv(args.metadata, index_col=0, low_memory=False)
    missing = [c for c in [args.group_by, *args.categorical, *args.numeric] if c not in metadata.columns]
    if missing:
        raise ValueError(f"Columns not found in {args.metadata.name}: {missing}")
    metadata[args.categorical] = metadata[args.categorical].astype(str)

    table = build_annotation_table(metadata, args.group_by, args.categorical, args.numeric, args.min_percent)
    if args.doublet_flag:
        table["caution_high_doublet_score"] = table[f"{args.doublet_flag}_median"] > args.doublet_cutoff
    table = sort_by_trailing_number(table)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output)
    print(f"{len(table)} rows written to {args.output}")


if __name__ == "__main__":
    main()

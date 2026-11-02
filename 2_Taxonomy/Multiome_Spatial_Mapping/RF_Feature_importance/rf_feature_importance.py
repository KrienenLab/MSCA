"""Rank genes by how much they help random forests predict cell type clusters.

Trains N_FORESTS sklearn random forests (TREES_PER_FOREST trees each, seeds 0..N_FORESTS-1) on all
cells of the input AnnData and averages their impurity-based feature importances (mean decrease in
Gini impurity). Writes one table, sorted from most to least important gene, to the output CSV.

Usage:
    python rf_feature_importance.py --input <cells_x_genes_counts.h5ad> --output <feature_importances.csv>
"""
import argparse
import time
from contextlib import contextmanager
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import LabelEncoder

LABEL_KEY = 'Cluster_v4'  # obs column the forests learn to predict
N_FORESTS = 40            # forests trained with seeds 0..N_FORESTS-1
TREES_PER_FOREST = 25     # 40 x 25 = 1000 trees in total


@contextmanager
def time_code(label='Code Block'):
    """Print when a block starts and how long it ran."""
    print(f'[{label}] start @: {datetime.now().strftime("%m_%d_%Y_%H:%M:%S")}', flush=True)
    start = time.time()
    yield
    elapsed = time.time() - start
    print(f'[{label}] runtime: {int(elapsed / 60)} min {int(elapsed % 60)} sec', flush=True)


def mean_feature_importance(X, y, n_forests, trees_per_forest):
    """Average the feature importances of n_forests random forests, each seeded by its index."""
    importances = []
    for seed in range(n_forests):
        with time_code(f'  - Train RF: {seed}'):
            rf = RandomForestClassifier(n_estimators=trees_per_forest, n_jobs=-1, random_state=seed)
            rf.fit(X, y)
            importances.append(rf.feature_importances_)
    # Forests have equal size, so the mean over forests is the mean over all trees
    return np.mean(importances, axis=0)


def importance_table(genes, importances):
    """Table of genes and their importances, most important first."""
    table = pd.DataFrame({'Gene': genes, 'Feature Importance': importances}, index=genes)
    return table.sort_values(by='Feature Importance', ascending=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--input', required=True, help=f'h5ad of raw counts (cells x genes) with obs["{LABEL_KEY}"]')
    parser.add_argument('--output', required=True, help='CSV to write the gene ranking to')
    args = parser.parse_args()

    # Load counts; uint16 holds the counts exactly and keeps memory down on ~680k cells
    adata = sc.read_h5ad(args.input)
    adata.X = adata.X.astype(np.uint16)
    print(adata)
    y = LabelEncoder().fit_transform(adata.obs[LABEL_KEY])

    # Train all forests on all cells: the goal is a gene ranking, not held-out accuracy
    with time_code('Train RFs'):
        importances = mean_feature_importance(adata.X, y, N_FORESTS, TREES_PER_FOREST)

    # Save the ranking
    table = importance_table(adata.var_names, importances)
    print(table.head(10))
    table.to_csv(args.output)


if __name__ == '__main__':
    main()

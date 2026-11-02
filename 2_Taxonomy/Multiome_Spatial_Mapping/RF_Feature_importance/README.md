# Random forest feature importance

Ranks genes by how much they help random forests predict `Cluster_v4` labels in the scRNA-seq data.
`rf_feature_importance.py` trains 40 sklearn random forests of 25 trees each (1000 trees, seeds 0-39) on all cells
and averages their impurity-based feature importances.

## Run

```bash
python rf_feature_importance.py --input <counts.h5ad> --output <feature_importances.csv>
```

`--input` is an h5ad of raw counts (cells x genes) with cluster labels in `obs["Cluster_v4"]`.

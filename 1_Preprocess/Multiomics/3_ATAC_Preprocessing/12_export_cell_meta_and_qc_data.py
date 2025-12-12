#!/usr/bin/env python
import anndata as ad
import pandas as pd
from pathlib import Path

pdir           = '/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster'
dir_count_mat  = Path(pdir, "data/atac/atac_count_matrix")
path_in        = Path(dir_count_mat, 'cell_by_consensus_peaks.h5ad')
path_out       = Path(pdir, "data/atac/atac_cell_meta.csv")

adata = ad.read_h5ad(path_in, backed='r')
adata.obs.to_csv(path_out)

# Also calculate total reads in peaks (for later FRIP calculation)
# -> (just do at cluster-level here; aggregate cells later for the total reads/denominator)
path_in = Path(dir_count_mat, 'Cluster_v4_by_consensus_peaks.h5ad')
adata = ad.read_h5ad(path_in)

Path(pdir, "qc/data").mkdir(parents=False, exist_ok=True)

cluster_sums = adata.X.sum(axis=1)
pd.DataFrame({
    "marm_cluster_id": adata.obs_names,
    "reads_in_peaks" : cluster_sums.flatten().tolist()[0]
}).to_csv(Path(pdir, "qc/data/atac_reads_in_consensus_peaks_by_cluster.csv"), index=False)


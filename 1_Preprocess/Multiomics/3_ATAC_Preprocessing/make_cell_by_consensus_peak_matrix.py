#!/usr/bin/env python3
from pathlib import Path
import os
import snapatac2 as snap

pdir           = '/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster'
path_atac_h5ad = Path(pdir, 'data/atac/atac_h5ads/_dataset.h5ads')
path_peaks     = Path(pdir, 'data/atac/atac_consensus_peaks.bed')
dir_count_mat  = Path(pdir, "data/atac/atac_count_matrix")
path_out       = Path(dir_count_mat, 'cell_by_consensus_peaks.h5ad')

if not dir_count_mat.is_dir():
    os.mkdir(dir_count_mat)

print("Importing concatenated h5ad...")
adata = snap.read_dataset(path_atac_h5ad)

print("Generating peak matrix...")
atac_reads_by_peaks = snap.pp.make_peak_matrix(
    adata, 
    file              = str(path_out), 
    peak_file         = str(path_peaks),
    chunk_size        = int(5e4),
    counting_strategy = "insertion"
)
adata.close()
print("Script finished.")

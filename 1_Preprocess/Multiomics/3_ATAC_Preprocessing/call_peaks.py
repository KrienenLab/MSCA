#!/usr/bin/env python3
#
# Note: with slurm job & 1 cpu core for calling peaks, took 2 days
import anndata as ad
import pandas as pd
import numpy as np
from pathlib import Path
import os
import shutil
import snapatac2 as snap
import re
import polars as pl

# NUM_WORKERS=4 # this doesn't work, it's possible with mode='r' and getting rid of any output writing it could work but no time...

pdir = '/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster'
path_atac_h5ad = Path(pdir, 'data/atac/atac_h5ads/_dataset.h5ads')
path_cellr_ref = '/jukebox/krienen/genomic_annotations/cellranger_ref/mCalJa1.2.pat.X_rs_mitos2/mCalJac1-pat-X_rs_mitos2'
# see bed file output names at end of script

dir_macs3_tmp = Path(pdir, ".macs3_tmp")

if not dir_macs3_tmp.is_dir():
    os.mkdir(dir_macs3_tmp)

# --------------------------------------------------
# Genome info
# --------------------------------------------------

print("Importing genome info...")

chrom_sizes = pd.read_csv(Path(path_cellr_ref, 'star/chrNameLength.txt'), sep='\t', header=None)
chrom_sizes = dict(zip(chrom_sizes[0], chrom_sizes[1]))
path_fasta = Path(path_cellr_ref, 'fasta/genome.fa')
if not path_fasta.is_file():
    path_fasta = Path(path_cellr_ref, 'fasta/genome.fa.gz')

path_gtf = Path(path_cellr_ref, 'genes/genes.gtf')
if not path_gtf.is_file():
    path_gtf = Path(path_cellr_ref, 'genes/genes.gtf.gz')

genome = snap.genome.Genome(fasta = path_fasta, annotation = path_gtf, chrom_sizes = chrom_sizes)


# --------------------------------------------------
# Data import
# --------------------------------------------------

print("Importing concatenated h5ad...")
adata = snap.read_dataset(path_atac_h5ad)


# --------------------------------------------------
# Sanitize grouping labels
# --------------------------------------------------

#
# snapatac2 will make temporary files with the literal name from the
# labels used to group the cells, so "/" are essentially forbidden
# and create errors
#
# so we will have to remove slashes, keep a safe label, and then
# revert to the original labels after
#
# actually nevermind, i'm going to overwrite the group labels, and these
# should be permanently changed eveywhere else
#

adata.obs['Subcortex_Group_v4'] = (
    adata.obs['Subcortex_Group_v4']
    .cast(pl.Utf8)
    .str.replace('/', '-')
)


# --------------------------------------------------
# Call peaks
# --------------------------------------------------

print("Calling peaks...")
snap.tl.macs3(
    adata, 
    groupby = 'Subcortex_Group_v4', 
    key_added = 'macs3_Subcortex_Group_v4', 
    n_jobs = 1,
    tempdir = dir_macs3_tmp
)

# I meant to keep the unmerged peak set and set a new _merged one, but instead it's overwritten
adata.uns['macs3_Subcortex_Group_v4'] = snap.tl.merge_peaks(adata.uns['macs3_Subcortex_Group_v4'], chrom_sizes=genome)

## Setup bed file
consensus_peaks = adata.uns['macs3_Subcortex_Group_v4'].to_pandas()
consensus_peaks["chr"]   = [peak.split(":")[0] for peak in consensus_peaks["Peaks"]]
consensus_peaks["start"] = [peak.split(":")[1].split("-")[0] for peak in consensus_peaks["Peaks"]]
consensus_peaks["end"]   = [peak.split(":")[1].split("-")[1] for peak in consensus_peaks["Peaks"]]

# >>> consensus_peaks.shape
# (1777056, 201)

# Write to a BED file
with open(Path(pdir, "data/atac/atac_consensus_peaks.bed"), "w") as f:
    consensus_peaks[["chr", "start", "end"]].to_csv(f, sep="\t", header=False, index=False)

with open(Path(pdir, "data/atac/atac_consensus_peaks_celltypes.bed"), "w") as f:
    consensus_peaks.to_csv(f, sep="\t", header=True, index=False)

adata.close()
shutil.rmtree(dir_macs3_tmp)
print("Script finished.")

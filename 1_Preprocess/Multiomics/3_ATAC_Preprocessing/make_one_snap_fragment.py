#!/usr/bin/env python3
#
# Run with positional argument of `sample_name`. Script will find the fragment path and the barcoded_cell_sample_name
# based on the library qc table, and will curated cell barcodes from the RNA.
import sys
import anndata as ad
import pandas as pd
import re
import polars
from pathlib import Path
import snapatac2 as snap

sample = sys.argv[1]
print(f"Script started for sample: {sample}")

# Project setup
pdir = '/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster'
dir_anndatas = Path(pdir, 'data/atac/atac_h5ads/anndatas')
dir_anndatas.mkdir(parents=True, exist_ok=True)

# Paths
# rna_path = Path(pdir, 'data/rna/rna_raw.h5ad')
rna_path = Path(pdir, 'data/rna/rna_raw_v4_freeze.h5ad')
qc_csv = Path(pdir, 'qc/fragment_jobs_to_run.csv')
cellr_ref = '/jukebox/krienen/genomic_annotations/cellranger_ref/mCalJa1.2.pat.X_rs_mitos2/mCalJac1-pat-X_rs_mitos2'

# Metadata
library_qc = pd.read_csv(qc_csv).set_index('sample_name')
fragment_path = library_qc.loc[sample, 'fragment_path']
barcoded_name = library_qc.loc[sample, 'barcoded_cell_sample_name']

# keeping only the pipeline standardized columns, plus the "Cluster_v4",
# which were manually added
cols_keep = [
    'barcoded_cell_sample_name', 'tissue_name', 'dissociated_cell_sample_name', 
    'donor_name', 'mit_name', 'sample_name', 'donor_slab_tile', 'hemisphere', 
    'alignment_job_id', 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.Rpl', 
    'percent.Hbb', 'scDblFinder.class', 'scDblFinder.score', 'scDblFinder.weighted', 
    'scDblFinder.cxds_score', 'S.Score', 'G2M.Score', 'Phase',
    'Cluster_v4'
]

# Genome
chrom_sizes = pd.read_csv(Path(cellr_ref, 'star/chrNameLength.txt'), sep='\t', header=None)
chrom_sizes = dict(zip(chrom_sizes[0], chrom_sizes[1]))
fasta = Path(cellr_ref, 'fasta/genome.fa')
gtf = Path(cellr_ref, 'genes/genes.gtf')
if not fasta.is_file():
    fasta = fasta.with_suffix('.fa.gz')
if not gtf.is_file():
    gtf = gtf.with_suffix('.gtf.gz')
genome = snap.genome.Genome(fasta=fasta, annotation=gtf, chrom_sizes=chrom_sizes)

# RNA metadata
print("Importing RNA cell metadata...")
rna_meta = ad.read_h5ad(rna_path, backed='r').obs
sample_meta = rna_meta.loc[rna_meta.sample_name == sample, cols_keep]

print("Processing and checking cell metadata...")
for cname in sample_meta.columns:
    col = sample_meta[cname]
    orig_dtype = col.dtype
    dtype_str = str(orig_dtype)

    # Attempt to coerce object columns into atomic types
    if pd.api.types.is_object_dtype(orig_dtype):
        print(f"→ Attempting to infer dtype for object column '{cname}' …")
        try:
            col_converted = pd.to_numeric(col, errors='raise')
            sample_meta[cname] = col_converted.astype(float)  # preserves NaNs
            print(f"   ↳ Inferred numeric (float) for '{cname}'")
        except Exception:
            try:
                col_lower = col.astype(str).str.lower()
                if col_lower.isin(['true', 'false', 'nan', 'none']).all():
                    sample_meta[cname] = col_lower.map({'true': True, 'false': False}).astype('boolean')
                    print(f"   ↳ Inferred boolean for '{cname}'")
                else:
                    raise ValueError
            except Exception:
                sample_meta[cname] = col.astype(str).replace({'nan': pd.NA, 'None': pd.NA})
                print(f"   ↳ Kept as string for '{cname}' (converted object to string)")

    # Ensure type is Polars/Arrow-safe
    col = sample_meta[cname]
    if pd.api.types.is_numeric_dtype(col):
        print(f"→ Ensuring numeric column '{cname}' (dtype: {dtype_str}) allows NaN")
        sample_meta[cname] = col.astype(float)  # ensure float for missing

    elif pd.api.types.is_bool_dtype(col):
        print(f"→ Converting boolean column '{cname}' to nullable Boolean")
        sample_meta[cname] = col.astype('boolean')  # Arrow-compatible

    elif pd.api.types.is_string_dtype(col):
        sample_meta[cname] = col.replace({'nan': pd.NA, 'None': pd.NA})
        print(f"→ Ensuring string column '{cname}' is clean")

    elif pd.api.types.is_categorical_dtype(col):
        sample_meta[cname] = col.astype(str).replace({'nan': pd.NA})
        print(f"→ Converted categorical '{cname}' to string")

    else:
        print(f"[WARNING] Column '{cname}' has unsupported dtype ({dtype_str}) → converting to string")
        sample_meta[cname] = col.astype(str).replace({'nan': pd.NA})

for colname, col in sample_meta.items():
    if sample_meta[colname].isna().all():
        print(f"[WARNING] Entire column '{colname}' is NA — this may cause failure.")

# Fragment import
print("Importing fragments...")
adata = snap.pp.import_fragments(
    fragment_path,
    chrom_sizes=genome,
    file=str(dir_anndatas / f"{barcoded_name}.h5ad"),
    sorted_by_barcode=False,
    whitelist=[re.sub("-.*", "-1", x) for x in sample_meta.index]
)
adata.obs_names = [re.sub("1$", barcoded_name, x) for x in adata.obs_names]
# Don't do this:
#   new_meta = polars.DataFrame(sample_meta.loc[adata.obs_names, :])
# Instead:
#   Select subset and convert to dict-of-lists to avoid dtype inference ambiguity,
#   and explicitly construct polars.DataFrame from well-typed dict
meta_subset = sample_meta.loc[adata.obs_names].to_dict(orient="list")
new_meta = polars.DataFrame({k: polars.Series(k, v) for k, v in meta_subset.items()})

adata.obs = new_meta.hstack(adata.obs[:])
adata.uns['genome'] = 'mCalJac1-pat-X_rs_mitos2'
adata.uns['species'] = 'Callithrix jacchus'

print("Calculating metrics...")
snap.metrics.tsse(adata, genome)
snap.metrics.frag_size_distr(adata)
print("Closing file...")
adata.close()
print("Script finished.")

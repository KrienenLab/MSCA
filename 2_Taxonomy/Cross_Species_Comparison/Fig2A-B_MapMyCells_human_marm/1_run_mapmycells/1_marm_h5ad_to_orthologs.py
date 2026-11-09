# Subset a marmoset h5ad to the genes with orthologs in another species and re-index the genes
# by the ortholog's gene ID or gene name, so the data can be mapped with MapMyCells to a
# reference from that species.
#
# Args:
# 1) ortholog tables dir
# 2) ortholog species name
# 3) input h5ad path
# 4) output directory
# 5) gene type used for the new index ('gene_id' or 'gene_name')

REMOVE_METADATA = False

print("##### Begin script marm_h5ad_to_orthologs.py #####")

from pathlib import Path
import sys
import anndata as ad
import pandas as pd

ORTHO_DIR = Path(sys.argv[1])
ORTHO_SPECIES = sys.argv[2]
in_path = Path(sys.argv[3])
out_path = Path(sys.argv[4])
gene_type = sys.argv[5]

if out_path.exists() is False:
    out_path.mkdir(parents=True, exist_ok=True)

out_data_path = Path(out_path, f"rna_{ORTHO_SPECIES}_ortho_{gene_type}.h5ad")
if out_data_path.exists():
    print(f"Output file {out_data_path} already exists. Exiting...")
    sys.exit()

# Ortholog table: marmoset gene name -> ortholog gene stable ID and gene name
ORTHO_PATH = Path(ORTHO_DIR, f"marm_{ORTHO_SPECIES}_ortholog_tables/marm_{ORTHO_SPECIES}_orthologs.csv")

# the table's column names use the capitalized species name
ORTHO_SPECIES_CAP = ORTHO_SPECIES.capitalize()

df_ortho = pd.read_csv(ORTHO_PATH)
df_ortho.set_index('Gene name', inplace = True)
df_ortho['Gene name_id'] = df_ortho.index
df_ortho = df_ortho[[f'{ORTHO_SPECIES_CAP} gene stable ID',f'{ORTHO_SPECIES_CAP} gene name','Gene name_id']]

# Import data
print("Importing h5ad...")
adata = ad.read_h5ad(in_path)

# optionally drop all obs columns to save space
if REMOVE_METADATA:
    adata.obs = pd.DataFrame(index = adata.obs.index)

# keep genes with orthologs, add the ortholog IDs/names and re-index
print("Subsetting and re-indexing...")
gene_idx_keep = adata.var.index.isin(df_ortho.index)
adata = adata[:, gene_idx_keep] 

if REMOVE_METADATA:
    adata.var = pd.DataFrame(index = adata.var.index)
    
adata.var = adata.var.join(df_ortho, how="left")

if gene_type == "gene_id":
    adata.var.set_index(f'{ORTHO_SPECIES_CAP} gene stable ID', inplace = True)
elif gene_type == "gene_name":
    adata.var.set_index('Gene name_id', inplace = True)

# several marmoset genes can share an ortholog; make the index unique
adata.var_names_make_unique()

print(f"Exporting finished h5ad file, with {adata.shape[0]} cells and {adata.shape[1]} genes.")
adata.write(out_data_path)

print("Script finished.")

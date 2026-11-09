#!/usr/bin/env python
"""scVI integration of single-nucleus datasets across species (template used for Fig. S6).

Model parameters (n_hidden, n_latent, n_layers, number of highly variable genes, training epochs) of each
integrated dataset are listed in the heading of Supplementary Table S4; set them with the arguments below.
The defaults are those of the basal ganglia integration of marmoset, macaque and human.

Steps
    1. Gene names of every dataset are replaced by the symbol of the orthologous gene (all_gene_ids.csv) and the
       datasets are restricted to the genes they share. Raw counts are taken from .raw if present, else from .X.
    2. The datasets are concatenated (obs column 'Study' = dataset name).
    3. Highly variable genes are selected (seurat_v3, per Study) and scVI is trained with 'Study' as batch (gene
       dispersion is modeled per Study) and the donor column as categorical covariate.
    4. The latent space (X_scVI) and a UMAP computed from it (X_umap_integrated) are saved.

Usage
    python scvi_integration.py --out-dir scvi_output \
        --dataset marmoset marmoset.h5ad Callithrix_jacchus \
        --dataset macaque macaque.h5ad Macaca_nemestrina \
        --dataset human human.h5ad Homo_sapiens \
        --donor-column Donor_label

    --dataset NAME H5AD SPECIES   repeat for every dataset; SPECIES is the name used in the 'species' column of
                                  all_gene_ids.csv (e.g. Homo_sapiens, Callithrix_jacchus, Mus_musculus)
    --gene-table                  ortholog table (default: all_gene_ids.csv next to this script)
    --donor-column                obs column with the donor of every cell, present in every dataset
Outputs (--out-dir): model/, X_scVI.csv, X_umap_integrated.csv, adata_int.h5ad
"""

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
import scvi
from scipy.sparse import csr_matrix

script_dir = Path(__file__).parent.resolve()


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dataset", nargs=3, action="append", required=True, metavar=("NAME", "H5AD", "SPECIES"))
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--gene-table", type=Path, default=script_dir / "all_gene_ids.csv")
    parser.add_argument("--donor-column", default="Donor_label")
    parser.add_argument("--n-top-genes", type=int, default=4000)
    parser.add_argument("--n-hidden", type=int, default=256)
    parser.add_argument("--n-latent", type=int, default=64)
    parser.add_argument("--n-layers", type=int, default=3)
    parser.add_argument("--max-epochs", type=int, default=300)
    parser.add_argument("--umap-min-dist", type=float, default=0.3)
    return parser.parse_args()


def to_orthologs(adata, ortholog_table, species):
    """Rename genes to the symbol of their ortholog; genes without one keep their name."""
    rows = ortholog_table[ortholog_table["species"] == species]
    if rows.empty:
        raise ValueError(f"Species '{species}' not found in the gene table")
    symbols = dict(zip(rows["gene"], rows["ortholog_symbol"]))
    adata.var["original_symbol"] = adata.var_names
    adata.var_names = pd.Index(adata.var_names.map(lambda gene: symbols.get(gene, gene)))
    adata.var_names_make_unique()
    return adata


def main():
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    ortholog_table = pd.read_csv(args.gene_table, usecols=["gene", "species", "ortholog_symbol"])

    # ----------------------------------------
    # Load data and map genes to orthologs
    # ----------------------------------------
    datasets = {}
    for name, path, species in args.dataset:
        adata = sc.read_h5ad(path)
        if args.donor_column not in adata.obs.columns:
            raise ValueError(f"Column '{args.donor_column}' not found in the obs of {name}")
        print(f"{name}: {adata.n_obs} cells x {adata.n_vars} genes ({species})")
        if adata.raw is not None:  # raw counts
            adata = ad.AnnData(X=adata.raw.X, obs=adata.obs, var=adata.raw.var)
        datasets[name] = to_orthologs(adata, ortholog_table, species)

    common_genes = sorted(set.intersection(*(set(a.var_names) for a in datasets.values())))
    print(f"{len(common_genes)} genes shared by all datasets")
    for name, adata in datasets.items():
        adata = adata[:, common_genes].copy()
        adata.X = csr_matrix(adata.X)
        adata.obs["Donor_label"] = adata.obs[args.donor_column].astype(str)
        datasets[name] = adata

    # ----------------------------------------
    # Merge; keep the counts in a layer
    # ----------------------------------------
    adata = ad.concat(datasets, join="inner", label="Study", index_unique="-")
    adata.layers["UMIs"] = adata.X.copy()

    # ----------------------------------------
    # scVI
    # ----------------------------------------
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=args.n_top_genes, layer="UMIs",
                                batch_key="Study", subset=False)
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    scvi.model.SCVI.setup_anndata(adata_hvg, layer="UMIs", batch_key="Study",
                                  categorical_covariate_keys=["Donor_label"])
    model = scvi.model.SCVI(adata_hvg, dispersion="gene-batch", n_hidden=args.n_hidden,
                            n_latent=args.n_latent, n_layers=args.n_layers)
    model.train(max_epochs=args.max_epochs)
    model.save(args.out_dir / "model", overwrite=True)

    # ----------------------------------------
    # Integrated space and UMAP
    # ----------------------------------------
    adata.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_scVI")
    adata.obsm["X_umap_integrated"] = sc.tl.umap(adata, min_dist=args.umap_min_dist, copy=True).obsm["X_umap"]

    pd.DataFrame(adata.obsm["X_scVI"], index=adata.obs_names).to_csv(args.out_dir / "X_scVI.csv")
    pd.DataFrame(adata.obsm["X_umap_integrated"], index=adata.obs_names,
                 columns=["UMAP1", "UMAP2"]).to_csv(args.out_dir / "X_umap_integrated.csv")

    # obs columns with mixed types (e.g. NA and strings) cannot be written to h5ad
    for column in adata.obs.columns[adata.obs.dtypes == object]:
        adata.obs[column] = adata.obs[column].astype(str)
    adata.write_h5ad(args.out_dir / "adata_int.h5ad")
    print(f"Saved to {args.out_dir}")


if __name__ == "__main__":
    main()

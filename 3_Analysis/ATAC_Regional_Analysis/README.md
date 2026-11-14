# Regional chromatin accessibility analysis

Analysis and plotting code for two analyses of ATAC signal by brain region:

- **RNA and ATAC co-enrichment by region and neurotransmitter** (Figure S7E): whether GABAergic and
  glutamatergic neurons of the same region share enriched genes and peaks.
- **Regional ATAC marker peaks** (Figures 7C and S16G), and the per-cluster score of the midbrain peak set
  handed to the spatial plots in Figures 7D-E.

The eRegulon analyses are in [`../ATAC_Regulon_Analysis`](../ATAC_Regulon_Analysis).

## Reading this code

- The code is published as it was run. It is not a pipeline that runs end to end: the scripts were run once
  each and the notebooks were run interactively. The notebooks include exploratory analyses that are not in
  the paper.
- Paths to files in this repository are built from the repository root: `here::here()` in R (the `here`
  package), `Path(__file__)` in Python and `git rev-parse --show-toplevel` in shell, so run the code from
  inside a clone. Large inputs that are not in the repository (see "Inputs") keep the path they had on the
  lab volume, `/jukebox/krienen` on the cluster or `/Volumes/krienen` on a laptop; point those at your copies.
- Final figure panels were edited in Illustrator.

## Layout

```
common_pipeline/               step 1: annotate and export the Cluster-level count matrices
DESeq2_AllRegionNT/            steps 2-3 for ATAC: differential accessibility, then co-enrichment
DESeq2_AllRegionNT_RNA/        steps 2-3 for RNA, and the notebook that makes Figure S7E
RegionalPeaksBySpecificity/    regional marker peaks
```

## Run order

| Step | Script | What it does |
| --- | --- | --- |
| 1 | `common_pipeline/0_update_cluster_region_annotation.py` | Bring the Class, Group and neurotransmitter columns of the cluster-level region file in line with the annotation table. |
| 1 | `common_pipeline/1_copy_modify_atac_pseudobulk_h5ad.py`, `2_copy_modify_rna_pseudobulk_h5ad.py` | Copy the Cluster-level count matrices, join the annotation table and one dominant region per cluster, and export the matrix and cluster table as csv to `common_data/`. |
| 2 | `DESeq2_AllRegionNT[_RNA]/pipeline/1_*.py`, `2_run_deseq2.sh`, `batchscript_region_nt_deseq2.sh` | Code each region+neurotransmitter grouping vs. neurons from other regions, then run one Slurm job per grouping with DESeq2 + apeglm. |
| 3 | `DESeq2_AllRegionNT/pipeline/3_Regional_Gaba_Glut_Coenrichment.Rmd` | Score pairwise co-enrichment of peaks between groupings (top 1% of peaks by shrunken fold change); writes `data/shrLFC_q99_coenrich_analysis.csv`. |
| 3 | `DESeq2_AllRegionNT_RNA/pipeline/3_Regional_Gaba_Glut_RNA_Coenrichment.Rmd` | The same for genes (top 2% of genes); writes `data/shrLFC_q98_coenrich_analysis.csv`. |
| 4 | `DESeq2_AllRegionNT_RNA/pipeline/4_Heatmap_Gaba_Glut_Enrichment_RNA_and_ATAC.Rmd` | Reads both co-enrichment tables and makes Figure S7E. |
| 5 | `RegionalPeaksBySpecificity/pipeline/1_Find_Peaks_By_Specificity.Rmd` | Defines regional marker peaks by specificity, and exports the per-cluster score of a midbrain peak set. |
| 6 | `RegionalPeaksBySpecificity/pipeline/2_RegionalMarkerPeaks.Rmd` | Defines the regional marker peaks again with the later annotation table and draws the heatmap in Figures 7C and S16G. |

Step 2 calls `RunDESeq2apeglm_bulked_h5ad.R` from `3_Analysis/Nonneuron/FigS5H_microglia_DEG/`. The jobs here ran an
earlier version of that script, which took a single numerator group and a single denominator group; the call
is the same.

`2_RegionalMarkerPeaks.Rmd` has no imports of its own. It was run in the same R session as
`../ATAC_Regulon_Analysis/combined_analysis/notebooks/1_combined_analysis.Rmd`, after that notebook's Env,
Paths, Utilities and Imports sections, and it writes its plots under that notebook's plot directory
(`combined_analysis/plots/RegionalPeaksBySpecificity/`).

## Figure panel to code

| Panel | Output file | Code |
| --- | --- | --- |
| S7E | `DESeq2_AllRegionNT_RNA/Gaba-Glut_CoEnrichment/enr_rna_vs_enr_atac_gaba_glut_only.pdf` | `DESeq2_AllRegionNT_RNA/pipeline/4_Heatmap_Gaba_Glut_Enrichment_RNA_and_ATAC.Rmd` line 536 |
| 7C | MB and THM blocks of `RegionalPeaksBySpecificity/heatmap_top40_by_combined_rank_dropOlfC.pdf` | `RegionalPeaksBySpecificity/pipeline/2_RegionalMarkerPeaks.Rmd` line 715 |
| S16G | `RegionalPeaksBySpecificity/heatmap_top40_by_combined_rank_dropOlfC.pdf` (all regions) | same |
| 7D, 7E | spatial plots made by a co-author | per-cluster PC1 score of a midbrain peak set, `RegionalPeaksBySpecificity/data/MB_peaks_in_THM_clusters_cl_eigen.csv`, written by `RegionalPeaksBySpecificity/pipeline/1_Find_Peaks_By_Specificity.Rmd` line 1436 |

## Small tables kept with the code

| File | What it is |
| --- | --- |
| `DESeq2_AllRegionNT[_RNA]/data/all_region_nt.txt` | The region+neurotransmitter groupings that got a DESeq2 job (written by step 2, script 1). |
| `DESeq2_AllRegionNT/data/shrLFC_q99_coenrich_analysis.csv` | Pairwise co-enrichment of peaks (step 3); input to Figure S7E. |
| `DESeq2_AllRegionNT_RNA/data/shrLFC_q98_coenrich_analysis.csv` | Pairwise co-enrichment of genes (step 3); input to Figure S7E. |
| `RegionalPeaksBySpecificity/data/MB_peaks_in_THM_clusters.txt` | The midbrain marker peaks with high signal in a subset of THM GABA clusters (one k-means cluster of the midbrain marker peaks). |
| `RegionalPeaksBySpecificity/data/THM_clusters_with_MB_peaks.txt` | The THM GABA clusters with high signal at those peaks. |
| `RegionalPeaksBySpecificity/data/MB_peaks_in_THM_clusters_cl_eigen.csv` | PC1 score of that peak set for every cluster. |

## Inputs

These are not in the repository. Paths are relative to `/jukebox/krienen`. A leading `...` stands for
`marm_hmba_integration/250602_reprocess_and_recluster`.

| Input | Used by |
| --- | --- |
| `.../data/atac/atac_count_matrix/Cluster_v4_by_462k_peaks.h5ad` | step 1 |
| `.../data/rna/metacell_v2/subcortex/Cluster_v4_metacell.h5ad` | step 1 |
| `.../metadata/subcortex_anno_table_archive/subcortex_anno_table_20251008.csv` | step 1 |
| `.../metadata/subcortex_anno_table_archive/subcortex_anno_table_cluster_region_20251008.csv` (dominant region per cluster) | step 1 |
| SCENIC+ `eRegulon_direct.tsv` and `eRegulons_extended.tsv` for neurons (output of `1_Preprocess/Multiomics/5_SCENIC_Plus`). `eRegulon_direct.tsv` is read from the copy in `3_Analysis/HPF_Gradient_Analysis/data/`. | step 4 |
| Objects loaded by `../ATAC_Regulon_Analysis/combined_analysis/notebooks/1_combined_analysis.Rmd` | step 6 |

## Annotation versions

Steps 1-5 use the annotation table of 2025-10-08. Step 6 and the eRegulon analyses use later tables. Later
annotation updates were mostly Group renames, so they matter here only where an older Group name appears as a
label. Figures S7E and the scores behind Figures 7D-E were not regenerated with the later tables.

Steps 5 and 6 define regional marker peaks with the same method. Figures 7C and S16G show the peaks from
step 6, and the score behind Figures 7D-E was computed in step 5.

## Not included

- The DESeq2 result tables, the exported count matrices and the plot files.

# eRegulon analysis

Analysis and plotting code for Figures 3 and 5 and Supplemental Figures S12B-D, S13 and S16E-F, plus the
numerical outputs handed to the spatial plots in Figures 3E and 5C.

The regional chromatin accessibility analyses (Figures 7C, 7E, S7E and S16G) are in
[`../ATAC_Regional_Analysis`](../ATAC_Regional_Analysis).

## Reading this code

- The code is published as it was run. It is not a pipeline that runs end to end: the scripts in steps 1-3
  below were run once each, and the notebook was run interactively, section by section.
- The notebook includes exploratory analyses that are not in the paper. The tables below point to the chunk
  that writes each figure panel.
- Paths are absolute. They start with `/jukebox/krienen` (the lab volume on the cluster) or `/Volumes/krienen`
  (the same volume mounted on a laptop). `lab_members/mike/251022_HMBA_Spatial_ATAC` in a path is this
  directory, so `common_data/`, `AUCell/data/` and `combined_analysis/plots/` resolve here.
- Final figure panels were edited in Illustrator (labels, legends, rotation, cropping, "MSN" relabelled as
  "SPN", "Glia" relabelled as "Nonneurons").

## Layout

```
common_pipeline/        step 1: export Cluster-level count matrices
AUCell/pipeline/        step 2: score eRegulons per cluster
DESeq2_RegionNT_RNA/    step 3: differential expression by region and neurotransmitter
DESeq2_RegionNT_ATAC/   step 3: differential accessibility by region and neurotransmitter
combined_analysis/
    notebooks/          step 4: the main notebook
    data/               small input tables
    plots/              small tables the notebook writes and reads back (no plot files)
```

## Run order

| Step | Script | What it does |
| --- | --- | --- |
| 1 | `common_pipeline/1_export_rna_pseudobulk_csv.py`, `2_export_atac_pseudobulk_csv.py` | Export the Cluster-level raw count matrices (gene-by-cluster and peak-by-cluster) to `common_data/*.csv`. |
| 2 | `AUCell/pipeline/1_get_metacell_AUCell.R`, `2_get_metacell_AUCell_nonneuron.R` | Score SCENIC+ eRegulons per cluster (gene AUC from RNA, peak AUC from ATAC); writes `AUCell/data/score_{rna,atac}_auc[_nonneuron].csv`. |
| 3 | `DESeq2_RegionNT_{RNA,ATAC}/pipeline/1_*.py`, `2_run_deseq2.sh`, `batchscript_region_nt_deseq2.sh` | Code each region+neurotransmitter grouping vs. neurons from other regions, then run one Slurm job per grouping with DESeq2 + apeglm. |
| 4 | `combined_analysis/notebooks/1_combined_analysis.Rmd` | Main notebook: imports everything above and makes the panels. |

Step 3 calls `RunDESeq2apeglm_bulked_h5ad.R`. A copy of that script is in this repository at
`3_Analysis/Nonneuron/FigS5H_microglia_DEG/RunDESeq2apeglm_bulked_h5ad.R`. The jobs here ran an earlier
version of it that took a single numerator group and a single denominator group.

## Figure panel to code

Line numbers refer to `combined_analysis/notebooks/1_combined_analysis.Rmd`. Output files are relative to
`combined_analysis/plots/`.

### Figure 3

| Panel | Output file | Code |
| --- | --- | --- |
| 3A | schematic | none |
| 3B | `heatmap_regulon_by_class_neurons3.pdf` | line 1378, section "(!) Neuron class by regulon heatmap > Using min-max scaling" |
| 3C | UMAPs of the OTX2 / cCRE / GRIK1 triplet | not in this directory (see "Not included") |
| 3D | browser track | made in IGV desktop from a Class-level bigWig; no code |
| 3E | spatial plot made by a co-author | input is the OTX2 gene AUC per cluster from `AUCell/data/score_rna_auc.csv` (step 2) |
| 3F | `OTX2/OTX2_exp_by_group.pdf` | line 1893, "OTX2 > Find groups enriching OTX2" |
| 3G | `OTX2/redo_OTX2_target_exp_by_selected_neur_v_glia[_BIG].pdf` | lines 1994-2001, "Scatterplot expression" |
| 3H | `OTX2/network_otx2_top_targets_other_TFs_TF2G_1p5_clipped.pdf` | line 3178, "New network: diff peaks/genes, same-peak TFs" |

### Figure 5

| Panel | Output file | Code |
| --- | --- | --- |
| 5A | `GanglionicEminence/heatmap_zscore_cor_ge_rgnt_by_nt/GABA_and_MSN.pdf` | line 6521, "(!) Heatmap correlations by NT" |
| 5B, 5H | schematics | none |
| 5C | spatial plot made by a co-author | cluster lists in `combined_analysis/data/pvalb_sst_clusters/` |
| 5D | `STR_v_CTX_MGE/scatter_pvalb_and_sst_str_v_ctx_tf_exp[_LARGE].pdf` (PVALB half) | lines 9382-9390, "Scatterplot all TFs" |
| 5E | `STR_v_CTX_MGE/nkx2_network_complete_with_activity_forMGE_STR_GABA2.pdf` | line 10874, final "Plot" section |
| 5F | `GanglionicEminence/heatmap_zscore_cor_ge_rgnt_by_nt/Glut.pdf` | line 6489 (same data without the "Glut" suffix in labels: `eregulon_plots/heatmap_zscore_cor_rgnt_by_nt/Glut.pdf`, line 1680) |
| 5G | the LGE and TEL Glut block of the Figure S13A heatmap | line 6535 |
| 5I | `TEL_Glut_vs_LGE/scatter_tel_glut_v_lge_coenrich_regulon_z[_big,_rev].pdf` | lines 3427-3447, "(!) Scatterplot regulon, enrichment vs. background" |
| 5J | `TEL_Glut_vs_LGE/dotplot_enr_selected_terms_signal[2].pdf` | lines 4173-4219, "DotPlot selected terms" |
| 5K | network images downloaded from string-db, annotated in Illustrator | gene list exported at line 3959; selected terms at line 4130 |
| 5L | `CRACDL/cracdl_regulon_z_sideways_topbar.pdf` (rotated) | line 5751, "Plot regulons & enrichment together". The CPM version at line 5851 was not used. |

### Supplemental figures

| Panel | Output file | Code |
| --- | --- | --- |
| S12A | non-neuron eRegulon heatmap | not in this directory (see "Not included") |
| S12B | `OTX2/redo_OTX2_target_lfc_gene_v_peak_allpeaks[_BIG].pdf` | lines 2450-2457 |
| S12C | `OTX2/redo_OTX2_target_lfc_gene_v_peak_bestpeaks_notInTarget[_BIG].pdf` | lines 2694-2701 |
| S12D | `OTX2/redo_OTX2_peak_TF_coenrichment_volcano[_BIG].pdf` | lines 2918-2925 |
| S13A | `GanglionicEminence/heatmap_ereg_ge_rgnt_z_cor.pdf` | line 6535, "Complete" |
| S13B | `STR_v_CTX_MGE/scatter_pvalb_and_sst_str_v_ctx_tf_exp[_LARGE].pdf` (SST half) | lines 9382-9390 |
| S13C | `TEL_Glut_vs_LGE/scatter_tf_vs_reg_z_by_ge_rgnt_for_fav2_tfs[_BIG].pdf` | lines 3638-3643 |
| S13D | `TEL_Glut_vs_LGE/upset_selected_coenr_terms.pdf` | line 4277 |
| S13E | `CRACDL/cracdl_network_side_by_side.pdf` | line 6403, "(!) CRACDL centric network plot" |
| S16E | `THM_vs_MB_GABA/barplot_ereg_rgnt_cor_with_thm_gaba_flip_top10.pdf` | line 7185 |
| S16F | `THM_vs_MB_GABA/scatter_ereg_z_mb_vs_thm_gaba_BIG.pdf` | line 7223 |

## Small tables kept with the code

| File | What it is |
| --- | --- |
| `DESeq2_RegionNT_{RNA,ATAC}/data/all_region_nt.txt` | The region+neurotransmitter groupings that got a DESeq2 job (written by step 3, script 1). |
| `combined_analysis/data/pvalb_sst_clusters/*.txt` | PVALB and SST cluster IDs by region (Figures 5C, 5D, S13B). |
| `combined_analysis/plots/TEL_Glut_vs_LGE/tel-glut_lge_coenr_genes_in_coenr_regulons_using_zdiff.csv` and `..._allGenes.txt` | The 279 co-enriched genes exported for string-db (Figure 5I-K). |
| `combined_analysis/plots/TEL_Glut_vs_LGE/tel_glut_lge_genes_stringdb_enr_results/enrichment.all.tsv` | Enrichment results downloaded from string-db; read back by the notebook at line 3979 for Figures 5J and S13D. Permalink in `tel_glut_lge_genes_stringdb_enr_results_permalink.txt`. |
| `combined_analysis/plots/TEL_Glut_vs_LGE/selected_string_terms.csv` | The manually selected terms shown in Figures 5J and S13D. |

The last three sit under `plots/` because that is where the notebook reads and writes them.

## Inputs

Paths are relative to `/jukebox/krienen`. A leading `...` stands for
`marm_hmba_integration/250602_reprocess_and_recluster`.

| Input | Used by |
| --- | --- |
| `.../data/rna/metacell_v2/subcortex/Cluster_v4_metacell.h5ad` (raw counts layer) | steps 1 and 3 |
| `.../data/atac/atac_count_matrix/Cluster_v4_by_462k_peaks.h5ad` | steps 1 and 3 |
| `.../metadata/preprint_freeze/subcortex_anno_table_preprint.csv` | step 3 |
| `.../metadata/revision/subcortex_anno_table_updated_clean_20260218.csv` | step 4 |
| SCENIC+ `eRegulon_direct.tsv` and `eRegulons_extended.tsv` for neurons and non-neurons (output of `1_Preprocess/Multiomics/5_SCENIC_Plus`) | steps 2 and 4 |
| `genomic_annotations/cellranger_ref/mCalJa1.2.pat.X_rs_mitos2/.../genes/genes.gtf.gz` | step 4 (TSS positions) |

The two annotation tables give the same cluster groupings. Both keep the same 1,642 clusters
(`preprint_analysis_include` in one, `!inconclusive_mapping` in the other), and Class, neurotransmitter and
dominant region are identical for every cluster. They differ in eight Group names (72 clusters), colors and
figure-panel assignments. One of the 1,642 clusters (a THM Glut cluster) has no ATAC data, so ATAC analyses
use 1,641.

Some panels were drawn before the later annotation table existed, by the same code reading the preprint
table. Because the cluster groupings are identical, this affects only the eight Group names where they
appear as labels.

## Not included

- Figure 3C (UMAPs) and Figure S12A (non-neuron eRegulon heatmap) were not made by the code in this directory.
- The two AUCell score tables (about 15 MB each) are not included; step 2 regenerates them.
- Plot files are not included.

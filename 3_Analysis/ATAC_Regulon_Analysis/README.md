# eRegulon and regional chromatin accessibility analysis

Analysis and plotting code for Figures 3, 5, 7C, S7E, S12B-D, S13 and S16E-G, plus the numerical outputs
handed to the spatial plots in Figures 3E, 5C and 7E.

The code is published as it was run. It was not refactored and it is not a pipeline that runs end to end:
the notebooks were run interactively, and many chunks are exploratory. The tables below point to the chunk
that writes each figure panel.

## Where the code came from

This directory mirrors the working directory `251022_HMBA_Spatial_ATAC`, so the relative paths in the code
(`common_data/`, `AUCell/data/`, `combined_analysis/plots/...`) resolve against this directory. Absolute
paths in the code start with `/jukebox/krienen` (HPC) or `/Volumes/krienen` (the same volume mounted on a
laptop).

`legacy_251007/` holds the parts of an earlier working directory (`251007_HMBA_Spatial_ATAC`) that still
produced something in the paper. Everything else in that earlier directory was superseded by the code here.

Final figure panels were edited in Illustrator (labels, legends, rotation, cropping, "MSN" relabelled as
"SPN", "Glia" relabelled as "Nonneurons").

## Run order

| Step | Script | What it does |
| --- | --- | --- |
| 1 | `common_pipeline/1_export_rna_pseudobulk_csv.py`, `2_export_atac_pseudobulk_csv.py` | Export the Cluster-level raw count matrices (gene-by-cluster and peak-by-cluster) to `common_data/*.csv`. |
| 2 | `AUCell/pipeline/1_get_metacell_AUCell.R`, `2_get_metacell_AUCell_nonneuron.R` | Score SCENIC+ eRegulons per cluster (gene AUC from RNA, peak AUC from ATAC); writes `AUCell/data/score_{rna,atac}_auc[_nonneuron].csv`. |
| 3 | `DESeq2_RegionNT_{RNA,ATAC}/pipeline/1_*.py`, `2_run_deseq2.sh`, `batchscript_region_nt_deseq2.sh` | Code each region+neurotransmitter grouping vs. neurons from other regions, then one Slurm job per grouping running DESeq2 + apeglm. |
| 4 | `combined_analysis/notebooks/1_combined_analysis.Rmd` | Main notebook: imports everything above and makes most panels. |
| 5 | `combined_analysis/notebooks/2_RegionalMarkerPeaks.Rmd` | Regional ATAC marker peaks. Has no imports of its own: it is run in the same R session after the import and setup sections of notebook 1 (it uses `obs`, `atac_by_cl_cpm`, `pdir_plots`). |

Step 3 calls `RunDESeq2apeglm_bulked_h5ad.R`. A cleaned-up copy of that script is in this repository at
`3_Analysis/Nonneuron/FigS5H_microglia_DEG/RunDESeq2apeglm_bulked_h5ad.R`; the jobs here ran the earlier
version of the same script (single numerator group and single denominator group).

`1_combined_analysis.Rmd` is the working notebook `1_combined_analysis_CleaningUp3.Rmd`, renamed. It is the
last of several cleanup passes over the original notebook. Plots dated before March 2026 were written by
the earlier passes, which used the same code with different variable names and the preprint annotation
table; see "Inputs" for why that does not change the cluster groupings.

## Figure panel to code

Line numbers refer to `combined_analysis/notebooks/1_combined_analysis.Rmd` unless another file is named.
Output files are relative to `combined_analysis/plots/`.

### Figure 3

| Panel | Output file | Code |
| --- | --- | --- |
| 3A | schematic | none |
| 3B | `heatmap_regulon_by_class_neurons3.pdf` | line 1544, section "(!) Neuron class by regulon heatmap > Using min-max scaling" |
| 3C | UMAPs of the OTX2 / cCRE / GRIK1 triplet | not in this directory (see "Not included") |
| 3D | browser track | made in IGV desktop from a Class-level bigWig; no code |
| 3E | spatial plot made by a co-author | input is the OTX2 gene AUC per cluster from `AUCell/data/score_rna_auc.csv` (step 2) |
| 3F | `OTX2/OTX2_exp_by_group.pdf` | line 2068, "OTX2 > Find groups enriching OTX2" |
| 3G | `OTX2/redo_OTX2_target_exp_by_selected_neur_v_glia[_BIG].pdf` | lines 2169-2176, "Scatterplot expression" |
| 3H | `OTX2/network_otx2_top_targets_other_TFs_TF2G_1p5_clipped.pdf` | line 3365, "New network: diff peaks/genes, same-peak TFs" |

### Figure 5

| Panel | Output file | Code |
| --- | --- | --- |
| 5A | `GanglionicEminence/heatmap_zscore_cor_ge_rgnt_by_nt/GABA_and_MSN.pdf` | line 7292, "(!) Heatmap correlations by NT" |
| 5B, 5H | schematics | none |
| 5C | spatial plot made by a co-author | cluster lists in `combined_analysis/data/pvalb_sst_clusters/` |
| 5D | `STR_v_CTX_MGE/scatter_pvalb_and_sst_str_v_ctx_tf_exp[_LARGE].pdf` (PVALB half) | lines 10179-10187, "Scatterplot all TFs" |
| 5E | `STR_v_CTX_MGE/nkx2_network_complete_with_activity_forMGE_STR_GABA2.pdf` | line 11711, final "Plot" section |
| 5F | `GanglionicEminence/heatmap_zscore_cor_ge_rgnt_by_nt/Glut.pdf` | line 7257 (same data without the "Glut" suffix in labels: `eregulon_plots/heatmap_zscore_cor_rgnt_by_nt/Glut.pdf`, line 1855) |
| 5G | the LGE and TEL Glut block of the Figure S13A heatmap | line 7306 |
| 5I | `TEL_Glut_vs_LGE/scatter_tel_glut_v_lge_coenrich_regulon_z[_big,_rev].pdf` | lines 3615-3635, "(!) Scatterplot regulon, enrichment vs. background" |
| 5J | `TEL_Glut_vs_LGE/dotplot_enr_selected_terms_signal[2].pdf` | lines 4888-4934, "DotPlot selected terms" |
| 5K | network images downloaded from string-db, annotated in Illustrator | gene list exported at line 4147; selected terms at line 4845 |
| 5L | `CRACDL/cracdl_exp_regulon_z_sideways_topbar.pdf` (rotated) | line 6607, "(!) Change to using exp (CPM), not LFC" |

### Figure 7

| Panel | Output file | Code |
| --- | --- | --- |
| 7C | MB and THM blocks of `RegionalPeaksBySpecificity/heatmap_top40_by_combined_rank_dropOlfC.pdf` | `2_RegionalMarkerPeaks.Rmd` line 755 |
| 7D, 7E | spatial plots made by a co-author | per-cluster PC1 score of the MB peak set: `legacy_251007/RegionalPeaksBySpecificity/data/MB_peaks_in_THM_clusters_cl_eigen.csv`, written by `legacy_251007/RegionalPeaksBySpecificity/pipeline/1_Find_Peaks_By_Specificity.Rmd` line 1560. An earlier version of the spatial plots is `2_spatial_plots_thm_vs_mb_peaks.py` in the same folder. |

### Supplemental figures

| Panel | Output file | Code |
| --- | --- | --- |
| S7E | `Gaba-Glut_CoEnrichment/enr_rna_vs_enr_atac_gaba_glut_only.pdf` | `legacy_251007/DESeq2_AllRegionNT_RNA/pipeline/4_Heatmap_Gaba_Glut_Enrichment_RNA_and_ATAC.Rmd` line 539 |
| S12A | non-neuron eRegulon heatmap | not in this directory (see "Not included") |
| S12B | `OTX2/redo_OTX2_target_lfc_gene_v_peak_allpeaks[_BIG].pdf` | lines 2628-2635 |
| S12C | `OTX2/redo_OTX2_target_lfc_gene_v_peak_bestpeaks_notInTarget[_BIG].pdf` | lines 2873-2880 |
| S12D | `OTX2/redo_OTX2_peak_TF_coenrichment_volcano[_BIG].pdf` | lines 3097-3104 |
| S13A | `GanglionicEminence/heatmap_ereg_ge_rgnt_z_cor.pdf` | line 7306, "Complete" |
| S13B | `STR_v_CTX_MGE/scatter_pvalb_and_sst_str_v_ctx_tf_exp[_LARGE].pdf` (SST half) | lines 10179-10187 |
| S13C | `TEL_Glut_vs_LGE/scatter_tf_vs_reg_z_by_ge_rgnt_for_fav2_tfs[_BIG].pdf` | lines 3826-3831 |
| S13D | `TEL_Glut_vs_LGE/upset_selected_coenr_terms.pdf` | line 4992 |
| S13E | `CRACDL/cracdl_network_side_by_side.pdf` | line 7168, "(!) CRACDL centric network plot" |
| S16E | `THM_vs_MB_GABA/barplot_ereg_rgnt_cor_with_thm_gaba_flip_top10.pdf` | line 7956 |
| S16F | `THM_vs_MB_GABA/scatter_ereg_z_mb_vs_thm_gaba_BIG.pdf` | line 8002 |
| S16G | `RegionalPeaksBySpecificity/heatmap_top40_by_combined_rank_dropOlfC.pdf` | `2_RegionalMarkerPeaks.Rmd` line 755 |

## Small tables kept with the code

| File | What it is |
| --- | --- |
| `DESeq2_RegionNT_{RNA,ATAC}/data/all_region_nt.txt` | The region+neurotransmitter groupings that got a DESeq2 job (written by step 3, script 1). |
| `combined_analysis/data/pvalb_sst_clusters/*.txt` | PVALB and SST cluster IDs by region (Figures 5C, 5D, S13B). |
| `combined_analysis/plots/TEL_Glut_vs_LGE/tel-glut_lge_coenr_genes_in_coenr_regulons_using_zdiff.csv` and `..._allGenes.txt` | The 279 co-enriched genes exported for string-db (Figure 5I-K). |
| `combined_analysis/plots/TEL_Glut_vs_LGE/tel_glut_lge_genes_stringdb_enr_results/enrichment.all.tsv` | Enrichment results downloaded from string-db; read back by the notebook at line 4682 for Figures 5J and S13D. Permalink in `tel_glut_lge_genes_stringdb_enr_results_permalink.txt`. |
| `combined_analysis/plots/TEL_Glut_vs_LGE/selected_string_terms.csv` | The manually selected terms shown in Figures 5J and S13D. |

Two of these sit under `plots/` because that is where the notebook reads and writes them.

## Inputs

| Input | Used by |
| --- | --- |
| `marm_hmba_integration/250602_reprocess_and_recluster/data/rna/metacell_v2/subcortex/Cluster_v4_metacell.h5ad` (raw counts layer) | steps 1 and 3 |
| `.../data/atac/atac_count_matrix/Cluster_v4_by_462k_peaks.h5ad` | steps 1 and 3 |
| `.../metadata/preprint_freeze/subcortex_anno_table_preprint.csv` | step 3 |
| `.../metadata/revision/subcortex_anno_table_updated_clean_20260218.csv` | notebook 1 |
| SCENIC+ `eRegulon_direct.tsv` and `eRegulons_extended.tsv` for neurons and non-neurons (output of `1_Preprocess/Multiomics/5_SCENIC_Plus`) | step 2, notebook 1 |
| `genomic_annotations/cellranger_ref/mCalJa1.2.pat.X_rs_mitos2/.../genes/genes.gtf.gz` | notebook 1 (TSS positions) |

The two annotation tables give the same cluster groupings. Both keep the same 1,642 clusters
(`preprint_analysis_include` in one, `!inconclusive_mapping` in the other), and Class, neurotransmitter and
dominant region are identical for every cluster. They differ in eight Group names (72 clusters), colors and
figure-panel assignments. One of the 1,642 clusters (a THM Glut cluster) has no ATAC data, so ATAC analyses
use 1,641.

## `legacy_251007/`

Earlier versions of steps 1-3 and of the regional peak analysis, built on the annotation freeze of 2025-10-08
(`metadata/subcortex_anno_table_archive/subcortex_anno_table_20251008.csv` and
`subcortex_anno_table_cluster_region_20251008.csv`; `common_pipeline/0b_*_redo.py` writes the cluster-level
dominant regions for that freeze). Two results in the paper come from here and were not regenerated with the
later annotation:

- Figure S7E. `DESeq2_AllRegionNT/` (ATAC) and `DESeq2_AllRegionNT_RNA/` run the same DESeq2 + apeglm design
  as step 3. Notebook `3_*.Rmd` in each folder scores pairwise co-enrichment and writes
  `data/shrLFC_q99_coenrich_analysis.csv` (ATAC, top 1% of peaks) and `data/shrLFC_q98_coenrich_analysis.csv`
  (RNA, top 2% of genes). `4_Heatmap_Gaba_Glut_Enrichment_RNA_and_ATAC.Rmd` reads both and makes the panel.
- Figure 7E. `RegionalPeaksBySpecificity/pipeline/1_Find_Peaks_By_Specificity.Rmd` defines the regional marker
  peaks with the same method as `2_RegionalMarkerPeaks.Rmd` and exports the PC1 score of the MB peak set per
  cluster. Figures 7C and S16G show the later rerun of the marker peaks, so the peak set behind 7E may differ
  slightly from the one drawn in 7C.

`common_pipeline/0a_tmp_update_anno_cluster_region.py` is the first attempt at those region assignments. It
reads a Google Sheets export (`251007_Marmoset_subcortex_annot_table.csv`) that is no longer on disk. The
assignments it carried are the `for_mike` column of the archived annotation table
`subcortex_anno_table_archive/subcortex_anno_table_20251008.csv`. They map one to one onto the
`Dominant_Region_Fine` column of `subcortex_anno_table_cluster_region_20251008.csv` (MSN and STR to STR, HPF to
HiF, Junk to Multi), which is what the later scripts read. Only the part of `0b_*_redo.py` above the "SCRAP"
marker runs.

Annotation updates between this freeze and the later ones were mostly Group renames, so they matter here only
where an old Group name appears as a label.

## Not included

- Figure 3C (UMAPs) and Figure S12A (non-neuron eRegulon heatmap) were not made by the code in this directory.
- The two AUCell score tables (about 15 MB each) are not included; step 2 regenerates them.
- Plot files are not included.
- The main notebook contains exploratory sections that did not go into the paper (for example the MSN Group
  eRegulon comparisons and the LHX6 networks). They were left in place.

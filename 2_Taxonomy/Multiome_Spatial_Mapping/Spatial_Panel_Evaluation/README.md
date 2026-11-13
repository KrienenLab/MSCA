# Spatial gene panel evaluation (Figure S2)

Analysis and plotting code for Figure S2: how well the 300-gene Xenium panel resolves the multiome taxonomy.
This directory mirrors the working directory `260216_Spatial_Panel_Eval`. The notebooks were run
interactively and are published as run; they contain exploratory sections that did not go into the paper.
Absolute paths start with `/Volumes/krienen`, which is `/jukebox/krienen` mounted on a laptop.

## Figure panel to code

Final panels were edited in Illustrator (labels and legends).

| Panel | Output file | Code |
| --- | --- | --- |
| S2A | `spatial_metacell_analysis/plots/panel_vs_hvg_gene_expression_percentile[2].pdf`, `panel_vs_hvg_gene_variance_percentile[2].pdf` | `spatial_metacell_analysis/pipeline/2_Verify_Spatial_Metacells.Rmd` lines 933-1012, "Check on panel gene variance" |
| S2B | `spatial_metacell_analysis/plots/scatter_genesets_cl_pairwise_gexp_dist_all_cl.pdf` | same notebook, line 329 |
| S2C | `mmc_multiome_selfmapping/plots/boxplot_region_precision_recall_f1.pdf` | `mmc_multiome_selfmapping/pipeline/1_Analyze_MMC_SelfMapping_Results.Rmd` line 1749 |
| S2D | `.../scatter_set_size_vs_avg_cluster_f1score.pdf` | same notebook, line 859 |
| S2E | `.../scatter_group_f1_4khvgs_vs_panel[2].pdf` | lines 1028-1076 |
| S2F | `random_forest_classifier/union_4k_hvgs_spatial_panel/plots/boxplot_rank_pct_importance_hvgs-vs-panel[2].pdf` | `random_forest_classifier/RF_FeatureImportancePlots.Rmd` lines 236-243 |
| S2G | `.../boxplot_cluster_and_group_precision_recall_f1[_flip].pdf` | `1_Analyze_MMC_SelfMapping_Results.Rmd` lines 1550-1574 |
| S2H | `.../boxplot_group_f1_by_class.pdf` | line 1693 |
| S2I | `.../scatter_precision_vs_palt1_confusion_mapping[2].pdf` | lines 2052-2088 |
| S2J | `.../rankscatter_avg_cluster_f1_by_group.pdf` | line 2930 |
| S2K | `spatial_metacell_analysis/plots/cell_abundance_comparisons_neuron_classes_regline.pdf` | `spatial_metacell_analysis/pipeline/3_Compare_Spatial_Multiome_CellCounts.Rmd` line 363 |
| S2L | `spatial_metacell_analysis/plots/cell_abundance_comparisons.pdf` | same notebook, line 266 |

## Upstream steps

| Step | Code | Output |
| --- | --- | --- |
| MapMyCells self-mapping of the multiome RNA with each gene set | not in this directory; run by a co-author | `aibs_fileshare/.../Subcortex/revisions/self_mapping_MMC_output/<gene set>/rna_raw_v4_freeze_noJunkClusters_mmc_results-cluster_v4.csv` |
| Random forest feature importance | `../RF_Feature_importance/rf_feature_importance.py` | `random_forest_classifier/union_4k_hvgs_spatial_panel/feature_importances_4104__1000_trees.csv` |
| Pseudobulk the Xenium data by cluster | `spatial_metacell_analysis/pipeline/1_run_spatial_pseudobulk.sh`, `pseudobulk_spatial_by_cluster.py` | `spatial_metacell_analysis/data/spatial_pseudobulk_Cluster_v4.csv` (not included) and `spatial_Cluster_v4_stats.csv` |

## Small tables kept with the code

| File | What it is |
| --- | --- |
| `mmc_multiome_selfmapping/data/MapMyCells_selfmapping_accuracy_by_{Cluster,Group,Region}.csv` | Precision, recall and F1 from the self-mapping, written by `1_Analyze_MMC_SelfMapping_Results.Rmd` lines 1353-1365. |
| `mmc_multiome_selfmapping/data/MapMyCells_selfmapping_truthtable_by_Group.csv` | True vs. predicted Group counts (line 1380). The Cluster-level table (6 MB) is not included. |
| `random_forest_classifier/union_4k_hvgs_spatial_panel/feature_importances_4104__1000_trees.csv` | Input to S2F. |
| `spatial_metacell_analysis/data/spatial_Cluster_v4_stats.csv` | Spatial cell counts per cluster; input to S2K-L. |

## Inputs

| Input | Used by |
| --- | --- |
| `lab_members/daisy/hmba/xenium_genes_2024.csv` (the 300-gene panel) | all notebooks |
| `marm_hmba_integration/250602_reprocess_and_recluster/data/rna/preprint_freeze/RNA_HVG_4K_preprint_freeze.csv` | all notebooks |
| `.../metadata/revision/subcortex_anno_table_updated_clean_20260218.csv`, `color_palette.csv` | all notebooks |
| `.../metadata/MSA_cluster_donor_table.csv` | S2K-L, self-mapping notebook |
| `.../data/rna/revision/rna_cell_metadata_with_preprint_umap_20260219.csv` | self-mapping notebook |
| Cluster-level RNA count matrix `common_data/Cluster_v4_by_genes_matrix.csv` | S2A-B; written by `3_Analysis/ATAC_Regulon_Analysis/common_pipeline/1_export_rna_pseudobulk_csv.py` |
| Xenium h5ad `Xenium-CJ23.56.004-marmoset-subcortex-cluster-v4_PREPRINT_FREEZE_20251022_subcortex.h5ad` | spatial pseudobulk |
| Gene set lists for S2D (NSForest markers, AnimalTFDB TFs, ABC atlas genes) | self-mapping notebook reads them from a laptop folder (`pdir_sets`); see "Gene sets" below |

## Gene sets

The self-mapping notebook reads five gene lists from `pdir_sets` (a laptop folder) to get the set sizes plotted
on the x-axis of S2D. A copy of that folder is in `mmc_multiome_selfmapping/data/gene_sets/`. The marker genes
MapMyCells actually used for each run are in
`self_mapping_MMC_output/<gene set>/query_markers_dispatch_*.json`.

| Gene set | File in `gene_sets/` | Size in notebook and Methods | Genes in the MapMyCells marker file | Where the list comes from |
| --- | --- | --- | --- | --- |
| Cluster NSForest markers | `Cluster_NSForest_Markers_Gene_List.txt` | 2,916 | 2,916 (same genes) | unique `nsforest_markers` in `metadata/preprint_freeze/supp_table/S2_Cluster_Info.csv` |
| Group NSForest markers | `Group_NSForest_Markers_Gene_List.txt` | 398 | 398 (same genes) | unique `nsforest_markers` in `metadata/preprint_freeze/taxonomy_group_class_preprint_nsforest.csv`; NSForest results in `data/rna/preprint_freeze/nsforest/Group/` |
| Marmoset TFs (AnimalTFDB v4.0) | `marm_TFs/Callithrix_jacchus_TF.txt` (not included; see `marm_TFs/readme.txt`) | 1,329 | 1,279 | AnimalTFDB download, kept at `genomic_annotations/AnimalTFDB4/callithrix_jacchus/Callithrix_jacchus_TF.txt` |
| Mouse ABC atlas MERFISH panel | `abc_wmb_tax_genes/abc_merfish_genes.txt` | 500 | 475 | written by `abc_wmb_tax_genes/260204_get_abc_tax_genes.Rmd` |
| Mouse ABC WMB taxonomy genes | `abc_wmb_tax_genes/abc_tax_genes.txt` | 194 | 190 | written by `abc_wmb_tax_genes/260204_get_abc_tax_genes.Rmd` |

`260204_get_abc_tax_genes.Rmd` reads `atlases/allen-brain-cell-atlas/metadata/WMB-taxonomy/20231215/cl.df_CCN202307220.xlsx`
and `.../metadata/WMB-10X/20231215/gene.csv`. The taxonomy genes are the words in the cluster, supertype,
subclass and class labels that are also mouse gene symbols.

For the last three sets, S2D plots the size of the full list, while the mapping used the smaller number of
genes in the marker file. The mouse gene names were converted to marmoset genes for the mapping, so the
marker files are not a simple subset of the upper-cased lists.

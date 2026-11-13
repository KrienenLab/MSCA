# MOXD1 eRegulons in astrocytes and neurons (Figures 4D and 4F)

Code for the z-score heatmaps in Figure 4D and the browser tracks in Figure 4F. This directory mirrors the
working directory `260303_MOXD1_Analysis`. The notebooks were run interactively and are published as run.
Absolute paths start with `/Volumes/krienen`, which is `/jukebox/krienen` mounted on a laptop.

| Panel | Output file | Code |
| --- | --- | --- |
| 4D (heatmaps) | `plots/heatmap_zscore_genes_regulons_peaks2.pdf` | `pipeline/1_Heatmaps_for_MOXD1_Regulons.Rmd` line 839, section "Heatmap w/ eregulons" |
| 4D (schematic) | drawn by hand | none |
| 4F | `browsershot/zoom_wholeRegion1*.pdf`, `zoom_Region1_zoom1_*.pdf`, `zoom_Region1_zoom2_enhancers*.pdf` | `browsershot/260305_Make_browsershot.Rmd` lines 390-632, section "Version 2 (with zooms)" |

Both panels were edited in Illustrator. In 4D, "Glia" was relabelled "Nonneurons" and the cCREs were renamed
from E1-E4 to cCRE1-cCRE4. Panel 4F was assembled from the separate whole-region and zoomed plots.

The first half of the browser shot notebook ("shot" to "shot4") is an earlier whole-locus layout that the
zoomed version replaced.

`browsershot/browser_shot_plotter/` is a copy of the plotting functions the notebook sources, from
<https://github.com/mdeber/browser_shot_plotter>.

## Inputs

| Input | Used by |
| --- | --- |
| Cluster-level raw count matrices `common_data/Cluster_v4_by_genes_matrix.csv` and `Cluster_v4_by_462k_peaks_matrix.csv` | heatmaps; written by `3_Analysis/ATAC_Regulon_Analysis/common_pipeline` |
| AUCell scores `AUCell/data/score_{rna,atac}_auc.csv` and `score_{rna,atac}_auc_nonneuron.csv` | heatmaps; written by `3_Analysis/ATAC_Regulon_Analysis/AUCell/pipeline` |
| SCENIC+ `eRegulon_direct.tsv` for neurons and non-neurons | heatmaps |
| `metadata/revision/subcortex_anno_table_updated_clean_20260218.csv` and `color_palette.csv` | heatmaps |
| CPM-normalized bigWigs by Group and Class, `data/atac/bw/Subcortex_{Group,Class}_v4/cpm/*.bw` (output of `1_Preprocess/Multiomics/3_ATAC_Preprocessing`) | browser tracks |
| `genomic_annotations/cellranger_ref/mCalJa1.2.pat.X_rs_mitos2/.../genes/genes.gtf.gz` | browser tracks |

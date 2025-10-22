# Marmoset Subcortical Cell Atlas
*for internal Krienen Lab users only*

## Project Summary
Understanding the cellular composition and gene regulatory landscape of the human brain remains a central challenge in neuroscience. Although whole-brain spatial and transcriptomic atlases exist for the mouse, the size and complexity of human and nonhuman primate brains have limited similarly comprehensive efforts in primates. Here, we present the most extensive, spatially resolved cellular atlas of the primate subcortex to date, encompassing the basal ganglia and all regions of the central nervous system outside of the neocortex and cerebellum. Leveraging the compact size of the common marmoset brain, we performed cellular-resolution spatial transcriptomics on intact hemispheres, imaging 5.6 million cells at ~200 μm-spaced intervals across the marmoset subcortex. In parallel, we generated paired RNA and ATAC single nucleus sequencing of ~700,000 nuclei sampled across densely tiled subcortical regions obtained from brains with atlas-registered structural and functional MRI. Cross-species alignment to an integrated basal ganglia taxonomy and existing mouse and human atlases enabled interoperable cell type annotations and highlighted divergent features between primates and mice. Comprehensive sampling in spatial and molecular modalities clarified cell type localization in border regions of the primate basal ganglia and nearby anatomically complex structures of the basal forebrain. Furthermore, transcriptomically defined cell types resolved discrete anatomical compartments, such as subdivisions of the primate-expanded medial pulvinar nucleus of the thalamus. We also find transcriptional similarities between cell types that often span anatomical transitions between areas, such as related cholinergic cell types in both basal ganglia and adjacent basal forebrain structures. Subcortical neurons exhibited strong regional identity but varied in heterogeneity across structures. For example, hypothalamic neurons shared greater transcriptional similarity with midbrain populations than with neighboring forebrain domains. We combined gene regulatory network inference with spatial mapping to identify regional patterns of both gene expression and cis-regulatory element activity, revealing anatomical- and lineage-associated determinants of cellular identity. In the basal ganglia, chromatin accessibility and gene expression hierarchies are largely aligned, but not perfectly concordant, suggesting that transcriptional identity is often (but not always) mirrored by cis-regulatory architecture. An interactive Cytosplore app enables multi-modal exploration of this resource.

## Data Avilability

### Marmoset Subcortical Cell Atlas

- **RNA**: `/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/data/rna/rna_clean_0715.h5ad`

  The most granular annotation level is stored as `Cluster_v4` and all annotation levels above are aggregated clusters. 
- **Cluster annotation**: `/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/metadata/subcortex_anno_table_20251019.csv`
  - **Oct 25, 2025 update**: preprint taxonomy frozen and stored as `/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/metadata/preprint_freeze/subcortex_anno_table_preprint.csv`

  This is work in progress and updated regularly. Relevant columns are `['Subcortex_Class_v4', 'Subcortex_Group_v4']`. 
- **ATAC**: `/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/data/atac/atac_count_matrix`
- **Spatial transcriptoimcs**: `/jukebox/krienen/data/MSCA/Spatial`
  - `...subcortex.h5ad` contains the hand drawn subcortical cells; other files without the subcortex name contain the full hemisphere； within the files, `Cluster_v4_name` or `CDM_Cluster_v4_label` corresponds to the mapped RNA `Cluster_v4` identity (through MapMyCells Flat)
    *Note: the spatial gene panel only contains 300 genes, designed to optimize for transcriptomic diversity in the basal ganglia. Clusters outside of the basal ganglia might be messy*
 
- *Interactively Viewing RNA and/or Spatial data*: Please navigate to this file `/jukebox/krienen/data/MSCA/readme.txt` for instruction on connecting to cirro servers hosted on scotty. 
  

### Cross Species Basal Ganglia Atlas

*Note: these files are all in synced AWS buckets with the Allen Brain Institute. Please do not overwrite these files, as none of the changes you make will be saved when the buckets get synced from cloud to local when we update.*

#### RNA
In each of the h5ad, you will find metadata columns named `['cluster_id', 'Neighborhood', 'Class', 'Subclass', 'Group']`. `cluster_id` is species-specific and the rest of the columns in the hierarchical order (most coarse to most granular) `Neighborhood, Class, Subclass, Group` are harmonized across species (effort led by Nelson Johansen and Yuanyuan Fu from AIBS). 

- Human:`/jukebox/krienen/aibs_fileshare/hmba-human-wg-802451596237-us-west-2/Aim1_Atlases/BasalGanglia/Human_basalganglia_HMBA_AIT19-5_anno_latest.h5ad`
- Macaque: `/jukebox/krienen/aibs_fileshare/hmba-macaque-wg-802451596237-us-west-2/Aim1_Atlases/BasalGanglia/Macaque_basalganglia_HMBA_AIT11-9_anno_latest.h5ad`
- Marmoset:  `/jukebox/krienen/aibs_fileshare/hmba-marmoset-wg-802451596237-us-west-2/Aim1_Atlases/BasalGanglia/Marmoset_basalganglia_anno_latest.h5ad`
- Xspecies: `/jukebox/krienen/aibs_fileshare/hmba-cross-species-wg-802451596237-us-west-2/HMBA_BG_Human_Macaque_Marmoset_BG_alignment.h5ad`

#### ATAC
- master folder: `/jukebox/krienen/aibs_fileshare/hmba-bican-sharing-802451596237-us-west-2/BasalGanglia_pre-print_ATAC`

#### Spatial Transcriptomics (Xenium)
- master folder: `/jukebox/krienen/aibs_fileshare/hmba-marmoset-wg-802451596237-us-west-2/spatial/for_data_and_tech_handoff/20250617_for_release/20250617_h5ad` contains the data and a .md file that explains all the metadata columns.
*This is the same data as the subcortical cell atlas, but subsetted to manually selected basal ganglia and then mapped to the marmoset BG RNA reference with the consensus taxonomy

## Code Availability
- Basic RNA data exploration: `MSCA/tutorial/explore_subcortical_data.ipynb`
- Cross-species integration example: `/jukebox/krienen/aibs_fileshare/241018_RE_AWS_bucket_for_HMBA_attachment/HMBA/align_hmba_marmoset_macaque_human.py`

### Utils
- Make metacells (pseudobulk): `/jukebox/krienen/marm_hmba_integration/250414_preprocess_br_pxr_250220_1/analysis/code/util/metacell`
- Mapmycells Example Pipeline: 


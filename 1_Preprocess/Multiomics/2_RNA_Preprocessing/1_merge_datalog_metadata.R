# This script will reformat the metadata in the krienen lab datalog such that,
# instead of having separate lines for both the RNA and ATAC data, there will be
# one line per sample and different columns associated with the RNA or ATAC.
#
# Some metadata fields are left behind here, but a large number are kept.
# [dropped columns are shown in a section near the end].
# ---
# Note: column name "studies" was changed to "study" between 250205 and 250318
# ---
suppressPackageStartupMessages(library(dplyr))
pdir <- '/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster'
path_datalog <- file.path(pdir, '250529_krienen_data_log_hmba.csv')
path_out <- file.path(pdir, 'metadata/merged_datalog_metadata.csv')

if (!dir.exists(file.path(pdir, 'metadata')))
    dir.create(file.path(pdir, 'metadata'))

# Paths to alignments -----------------------------------------------------

# Want to keep the metadata aligned with the actual aligned samples
# (Or should I include the QC fail in this as well?)
align_dirs <- c(
    '/jukebox/krienen/marm_hmba/2312xx/alignments/231224_cellranger-arc',
    '/jukebox/krienen/marm_hmba/240103/alignments/240110_cellranger-arc',
    '/jukebox/krienen/marm_hmba/240215/alignments/240307_cellranger-arc',
    '/jukebox/krienen/marm_hmba/240313/alignments/240323_cellranger-arc',
    '/jukebox/krienen/marm_hmba/240625/alignments/240807_cellranger-arc',
    '/jukebox/krienen/marm_hmba/240822_align_unmatched',
    '/jukebox/krienen/marm_hmba/241115/alignments/241212_cellranger-arc',
    '/jukebox/krienen/marm_hmba/250106/alignments/250114_cellranger-arc',
    '/jukebox/krienen/marm_hmba/250407/alignments/250529_cellranger-arc'
)
sample_names <- sapply(align_dirs, function(align_dir) {
    readLines(file.path(align_dir, "sample_names.txt"))
})

# keep the alignment date as the "alignment_job_id" (and join at end)
align_dates <- data.frame(
    sample_name = unlist(unname(sample_names)),
    align_dir = rep(names(sample_names), lengths(sample_names))
)
align_dates$alignment_job_id <- sub("_.*", "", sub(".*/", "", align_dates$align_dir))

# keep just the vector of sample_names for everything else
sample_names <- unlist(unname(sample_names))


# Join metadata -----------------------------------------------------------

# Note, these are the samples that probably have swapped sorted/unsorted
# 231207_HMBA_cjNutmeg_Slab6_Tile6_pooled_RNA1
# 231207_HMBA_cjNutmeg_Slab6_Tile6_pooled_RNA2
# 231207_HMBA_cjNutmeg_Slab6_Tile6_unsorted_RNA3
# 231207_HMBA_cjNutmeg_Slab6_Tile6_unsorted_RNA4

kdatalog <- read.csv(path_datalog)

# some sample(s) was submitted with a dot "." in the name, which had to be 
# modified
kdatalog$krienen_lab_identifier <- sub(
    "\\.", "p", kdatalog$krienen_lab_identifier
)

# unique(kdatalog$library_method)
# [1] "10xMultiome-RSeq" "10xMultiome-ASeq" "10xV3.1"          ""
#
# all of these samples have one of the "Multiome" library_method names
rna_datalog  <- subset(kdatalog, library_method == "10xMultiome-RSeq")
atac_datalog <- subset(kdatalog, library_method == "10xMultiome-ASeq")

# this is what is in find_samples_and_files.R; it doesn't work b/c the first
# tranche of data wasn't processed with that script, it was entered manually
#   rna_datalog$sample_name  <- sub("_RNA*", "", rna_datalog$krienen_lab_identifier)
#   atac_datalog$sample_name <- sub("_ATAC*", "", atac_datalog$krienen_lab_identifier)
#
# So this is updated
rna_datalog$sample_name  <- sub("_RNA(Seq|seq)*", "", rna_datalog$krienen_lab_identifier)
atac_datalog$sample_name <- sub("_ATAC(Seq|seq)*", "", atac_datalog$krienen_lab_identifier)

rna_datalog  <- subset(rna_datalog,  sample_name %in% sample_names)
atac_datalog <- subset(atac_datalog, sample_name %in% sample_names)

stopifnot(all(sample_names %in% rna_datalog$sample_name))
stopifnot(all(sample_names %in% atac_datalog$sample_name))

# Now, how to merge metadata from the two modalities

# Found in both (to keep, both should be identical)
col_identical <- c(
    "experiment_start_date",
    "mit_name",
    "donor_name",
    "tissue_name",
    "dissociated_cell_sample_name",
    "facs_population_plan",
    "cell_prep_type",
    "study",
    "expc_cell_capture",
    "enriched_cell_sample_quantity_count",
    "port_well",
    "barcoded_cell_sample_name"
    # "BG_containing_prob",
    # "BG_structures",
    # "Minor_BG",
    # "Other_Structures",
    # "All_structures"
)
stopifnot(
    all(col_identical %in% intersect(colnames(rna_datalog), colnames(atac_datalog)))
)
stopifnot(all(
    sapply(col_identical, function(x) all(rna_datalog[,x] == atac_datalog[,x]))
))

df_meta <- left_join(
    data.frame(sample_name = sample_names), 
    rna_datalog[, c(col_identical, "sample_name")], 
    by = "sample_name"
)

# Found in both, both should be kept but prefixed with RNA_ or ATAC_
col_both_prefix <- c(
    "library_creation_date",
    "library_prep_set",
    "library_name",
    "tapestation_avg_size_bp",
    "library_num_cycles",
    "lib_quantification_ng",
    "library_pool_name"
)

df_meta <- rna_datalog[, col_both_prefix] %>% 
    setNames(., paste0("RNA_", names(.))) %>% 
    cbind(sample_name = rna_datalog$sample_name) %>% 
    left_join(df_meta, ., by = "sample_name")

df_meta <- atac_datalog[, col_both_prefix] %>% 
    setNames(., paste0("ATAC_", names(.))) %>% 
    cbind(sample_name = atac_datalog$sample_name) %>% 
    left_join(df_meta, ., by = "sample_name")

# Found in RNA only (to be kept, no renaming)
col_rna_no_rename <- c(
    "cDNA_amplification_method",
    "cDNA_amplification_date",
    "amplified_cdna_name",
    "cDNA_pcr_cycles",
    "percent_cdna_longer_than_400bp",
    "cdna_amplified_quantity_ng",
    "cDNA_library_input_ng"
)

df_meta <- left_join(
    df_meta, 
    rna_datalog[, c(col_rna_no_rename, "sample_name")], 
    by = "sample_name"
)

# Only in RNA, but should still be prefixed with RNA_
col_rna_prefix <- c(
    "r1_index",
    "r2_index"
    # "r1_index_sequence",
    # "r2_index_sequence"
)

df_meta <- rna_datalog[, col_rna_prefix] %>% 
    setNames(., paste0("RNA_", names(.))) %>% 
    cbind(sample_name = rna_datalog$sample_name) %>% 
    left_join(df_meta, ., by = "sample_name")

# Found in ATAC only (to be kept, no renaming)
col_atac_no_rename <- c(
    "ATAC_index"
    # "ATAC_index_sequence_0",
    # "ATAC_index_sequence_1",
    # "ATAC_index_sequence_2",
    # "ATAC_index_sequence_3"
)

df_meta <- left_join(
    df_meta, 
    atac_datalog[, c(col_atac_no_rename, "sample_name")], 
    by = "sample_name"
)

stopifnot(all(sample_names == df_meta$sample_name))


# Identify dropped columns ------------------------------------------------

# names(kdatalog)[!names(kdatalog) %in% c(
#     col_identical,
#     col_both_prefix,
#     col_rna_no_rename,
#     col_rna_prefix,
#     col_atac_no_rename
# )]

# [1] "krienen_lab_identifier"             
# [2] "seq_portal"                         
# [3] "elab_link"                          
# [4] "tissue_name_old"                    
# [5] "enriched_cell_sample_container_name"
# [6] "enriched_cell_sample_name"          
# [7] "library_method"                     
# [8] "rna_amplification_pass_fail"        
# [9] "library_prep_pass_fail"             
# [10] "BG_containing_prob"                 
# [11] "BG_structures"                      
# [12] "Minor_BG"                           
# [13] "Other_Structures"                   
# [14] "All_structures"


# Annotate if FACS sorted -------------------------------------------------

# Was sample FACS sorted?
df_meta$facs_sorted <- df_meta$facs_population_plan != "no_FACS"

# > unique(df_meta$facs_population_plan)
# [1] "no_FACS"     "sorted_NeuN" "70/20/10"    "90/0/10"     "90/10/0"    
# [6] "70/0/30"     "DAPI"        "70/10/20"    "85/5/10"     "10/10/80" 

# 
# Addressing the swapped samples!
#
# Note: Here I'm only addressing two variables, not correcting all of the 
#       metadata between these. However, I'm not sure if there are any other
#       metadata entries that would be relevant, and the samples aren't being
#       renamed just yet.
#
# These are the samples that apparently reversed sorted vs. unsorted
# 
swapped_samples <- c(
    "231207_HMBA_cjNutmeg_Slab6_Tile6_pooled1",
    "231207_HMBA_cjNutmeg_Slab6_Tile6_pooled2",
    "231207_HMBA_cjNutmeg_Slab6_Tile6_unsorted3",
    "231207_HMBA_cjNutmeg_Slab6_Tile6_unsorted4"
)
idx_swapped <- df_meta$sample_name %in% swapped_samples
df_meta$facs_sorted[idx_swapped] <- !df_meta$facs_sorted[idx_swapped]

# Now will also edit the FACS population plan as well
swapped_entries <- unique(df_meta$facs_population_plan[idx_swapped])
df_meta$facs_population_plan[idx_swapped] <- ifelse(
    df_meta$facs_population_plan[idx_swapped] == swapped_entries[1],
    swapped_entries[2],
    swapped_entries[1]
)

# > df_meta[idx_swapped, c("sample_name", "facs_sorted", "facs_population_plan")]
#                                   sample_name facs_sorted facs_population_plan
# 16   231207_HMBA_cjNutmeg_Slab6_Tile6_pooled1       FALSE              no_FACS
# 17   231207_HMBA_cjNutmeg_Slab6_Tile6_pooled2       FALSE              no_FACS
# 18 231207_HMBA_cjNutmeg_Slab6_Tile6_unsorted3        TRUE              90/0/10
# 19 231207_HMBA_cjNutmeg_Slab6_Tile6_unsorted4        TRUE              90/0/10



# Change FACS slashes to dashes -------------------------------------------

# Lydia was concerned about date conversion

# > df_meta$facs_population_plan %>% unique
# [1] "no_FACS"     "sorted_NeuN" "70/20/10"    "90/0/10"     "90/10/0"     "70/0/30"     "DAPI"        "70/10/20"   
# [9] "85/5/10"     "10/10/80"    "12/13/75"    "11/18/71"    "84/6/10"     "84/7/9"      "80/10/10"    "100/0/0"    
# [17] "76/14/10"    "54/23/23"    "68/19/13"    "50/8/42"     "75/15/10"    "66/8/26"     "70/14/16"    "70/13/17"   
# [25] "50/0/50"     "70/19/11"

# > gsub("/", "-", df_meta$facs_population_plan) %>% unique
# [1] "no_FACS"     "sorted_NeuN" "70-20-10"    "90-0-10"     "90-10-0"     "70-0-30"     "DAPI"        "70-10-20"   
# [9] "85-5-10"     "10-10-80"    "12-13-75"    "11-18-71"    "84-6-10"     "84-7-9"      "80-10-10"    "100-0-0"    
# [17] "76-14-10"    "54-23-23"    "68-19-13"    "50-8-42"     "75-15-10"    "66-8-26"     "70-14-16"    "70-13-17"   
# [25] "50-0-50"     "70-19-11" 

df_meta$facs_population_plan <- gsub("/", "-", df_meta$facs_population_plan)


# Add donor_slab_tile + hemisphere coding ---------------------------------

# note that donor_slab_tile does not correspond 1:1 with tissue_name because
# they can come from different hemispheres
df_meta$donor_slab_tile <- mutate(
    df_meta,
    slab = strsplit(sample_name, "_") %>% unlist %>% grep("Slab", ., value = TRUE),
    # make tile whatever comes after slab
    tile = sapply(sample_name, function(x) {
        splits <- unlist(strsplit(x, "_"))
        idx_slab <- grep("Slab", splits)
        splits[idx_slab + 1]
    }),
    slab_tile = ifelse(
        tile == "", slab, paste0(slab, "_", tile)
    ),
    donor_slab_tile = paste0(mit_name, "_", slab_tile)
)$donor_slab_tile


### highlight cases where multiple hemispheres are present for a given tile
#
# > length(unique(df_meta$tissue_name))
# [1] 73
# > length(unique(df_meta$donor_slab_tile))
# [1] 69
#
# tst <- df_meta[, c("donor_slab_tile", "tissue_name")] %>% 
#     unique %>% 
#     .[order(.$donor_slab_tile), ]
# dupes <- tst$donor_slab_tile[duplicated(tst$donor_slab_tile)]
#
# > dupes
# [1] "cjNutmeg_Slab4_Tile3" "cjNutmeg_Slab5_Tile3" "cjNutmeg_Slab6_Tile3"
# [4] "cjNutmeg_Slab6_Tile6"
#
# > subset(tst, donor_slab_tile %in% dupes) %>% print.data.frame(row.names=F)
#       donor_slab_tile          tissue_name
#  cjNutmeg_Slab4_Tile3 CJ23.56.003.CX.44.03
#  cjNutmeg_Slab4_Tile3 CJ23.56.003.CX.04.03
#  cjNutmeg_Slab5_Tile3 CJ23.56.003.CX.05.03
#  cjNutmeg_Slab5_Tile3 CJ23.56.003.CX.45.03
#  cjNutmeg_Slab6_Tile3 CJ23.56.003.CX.06.03
#  cjNutmeg_Slab6_Tile3 CJ23.56.003.CX.46.03
#  cjNutmeg_Slab6_Tile6 CJ23.56.003.CX.06.06
#  cjNutmeg_Slab6_Tile6 CJ23.56.003.CX.46.06
#
###

### code hemisphere: lh <40, rh >=40, and bi >=90
df_meta$hemisphere <- as.numeric(
    sapply(strsplit(df_meta$tissue_name, "\\."), "[", 5)
)
df_meta$hemisphere <- ifelse(
    df_meta$hemisphere < 40,
    "lh",
    ifelse(df_meta$hemisphere < 90, "rh", "bi")
)


# Re-order some columns ---------------------------------------------------

stopifnot(all(
    tail(colnames(df_meta), 3) == c("facs_sorted", "donor_slab_tile", "hemisphere")
))
nc <- ncol(df_meta)
df_meta <- df_meta[, c(1:2, (nc-1):nc, 3:7, nc-2, 8:(nc-3))]


# Join in "alignment_job_id" ----------------------------------------------

df_meta <- left_join(
    df_meta, 
    align_dates[, c("sample_name", "alignment_job_id")],
    by = "sample_name"
)

# Export metadata ---------------------------------------------------------

# first sort
df_meta <- df_meta[order(df_meta$sample_name), ]
rownames(df_meta) <- NULL

cat("Exporting merged datalog metadata...\n")
write.csv(df_meta, path_out, row.names = FALSE)

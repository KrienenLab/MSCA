# July 24, 2025
#
# For making the data landing pages for the subcortical paper's data, Lauren
# (maybe with Suvvi) sent this "publication landing page template.xlsx" file
# which has this at the bottom:
#
# ---
# Please email us a separate lists of 1) raw filenames, associated aliquot local
# names and "tmp" GCP bucket names and 2) analysis filenames with associated
# aliquot local names and/or barcoded cell sample names as applicable along with
# "tmp" GCP bucket names. Please note, only data that’s released through
# Continuous Release will be included in the BDBags. If the data hasn’t been
# released or is only partly released by the time the manuscript is submitted,
# NeMO will update the BDBags once the release happens.
# ---
#
# I'm just going to make a single table out of this, for each sample, the
# raw (fastq) library names for each modality, the bucket it's from, and
# the barcoded_cell_sample_name (which is also the name of the h5ad file).
# 
# so start with sample_name, and then get:
# - barcoded_cell_sample_name
# - bucket name
# - fastq names
# 
# will use the libraries.csv files, which have already done most of the 
# heavy lifting (except for the barcoded_cell_sample_name).
# 
# and do only for samples actually in the subcortical data (passing QC)
suppressPackageStartupMessages(library(dplyr))
pdir <- '/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster/'
merged_datalog <- read.csv(file.path(pdir, 'metadata/merged_datalog_metadata.csv'))
qc_table <- read.csv(file.path(pdir, 'qc/library_qc_pass_fail_table.csv'))
qc_table <- subset(qc_table, RNA_all_qc_pass == "PASS")
# also from file "additional_qc_fails.csv" (RNA fail)
qc_table <- subset(qc_table, sample_name != "240207_HMBA_cjNutmeg_Slab4_Tile4_pooled1")

# not really the output object, just joining in barcoded_cell_sample_name
out <- qc_table[, c("sample_name", "align_dir")]
out <- left_join(
    out, 
    merged_datalog[, c("sample_name", "barcoded_cell_sample_name")], 
    by = "sample_name"
)

# now make the output dataframe
out_df <- do.call(rbind, lapply(1:nrow(out), function(i) {
    suppressWarnings(
        lib_tbl <- read.csv(
            file.path(out$align_dir[i], out$sample_name[i], "libraries.csv"), 
            header=TRUE
        )
    )
    # verify structure of libraries.csv file (one RNA, one ATAC raw name, in that
    # order)
    stopifnot(nrow(lib_tbl) == 2)
    stopifnot(identical(
        lib_tbl$library_type, 
        c("Gene Expression", "Chromatin Accessibility")
    ))
    
    # RNA, ATAC
    buckets <- strsplit(lib_tbl$fastqs, "/") %>% 
        lapply(function(x) grep("^nemo-tmp", x, value=TRUE)) %>% 
        unlist
    
    data.frame(
        raw_name   = lib_tbl$sample,
        tmp_bucket = buckets, 
        modality   = c("RNA", "ATAC"),
        barcoded_cell_sample_name = rep(out$barcoded_cell_sample_name[i], 2),
        internal_sample_name      = rep(out$sample_name[i], 2)
    )
}))

# and get rid of the ATAC files from the ATAC QC fail samples
ATAC_fails <- subset(qc_table, ATAC_all_qc_pass != "PASS")$sample_name
out_df <- subset(out_df, !((internal_sample_name %in% ATAC_fails) & modality == "ATAC"))
# > sum(qc_table$RNA_all_qc_pass == "PASS")
# [1] 146
# > sum(qc_table$ATAC_all_qc_pass == "PASS")
# [1] 133
stopifnot(nrow(out_df) == sum(qc_table$RNA_all_qc_pass == "PASS") + sum(qc_table$ATAC_all_qc_pass == "PASS"))

write.csv(
    out_df, file.path(pdir, 'metadata/MSA_library_sample_table.csv'), row.names = FALSE
)


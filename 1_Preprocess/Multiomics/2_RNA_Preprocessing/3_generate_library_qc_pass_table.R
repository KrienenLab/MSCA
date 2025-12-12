# For our library-level metrics, generate a table with all samples, showing 
# for each QC metric if the library is a pass or fail
# 
# Note that the sample will not be in the cellranger summary if it failed to
# align, but those samples will be apparent in this table
# 
# RNA fail if: 
# - Reads mapped confidently to transcriptome <50% (all pass here)
# - Reads in cells <60%
# 
# ATAC fail if:
# - ATAC unmapped read pairs >0.1
# - Fraction of high quality fragments in cells <0.5
# - Fraction of transposition events in peaks in cells <0.2
#
# Will not add a column for post-clustering RNA QC fail (which is currently 
# vaguely defined but based on a large fraction of cells being found in clusters
# that are exlcuded for being low quality).
# -> This should be done somewhere else
# 
# suppressPackageStartupMessages(library(dplyr))
pdir <- '/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster'
out_path <- file.path(pdir, 'qc/library_qc_pass_fail_table.csv')
path_cellr_summary <- file.path(pdir, 'metadata/cellranger-arc_summary.csv')
path_library_meta <- file.path(pdir, 'metadata/merged_datalog_metadata.csv')
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

cellr_summary <- read.csv(path_cellr_summary)

# make the qc folder if it doesn't already exist
if (!dir.exists(dirname(out_path)))
    dir.create(dirname(out_path))


# Initialize table --------------------------------------------------------

samples_by_align <- lapply(align_dirs, function(align_dir) {
    readLines(file.path(align_dir, "sample_names.txt"))
})
names(samples_by_align) <- align_dirs

out_table <- data.frame(
    sample_name = unname(unlist(samples_by_align)),
    align_dir = rep(names(samples_by_align), lengths(samples_by_align))
)

# This wouldn't necessarily be a problem for the code, but it would be 
# unexpected
stopifnot(all(cellr_summary$Sample.ID %in% out_table$sample_name))


# Find QC fails without alignment -----------------------------------------

early_fails <- unlist(lapply(align_dirs, function(align_dir) {
    sample_names <- readLines(file.path(align_dir, "sample_names.txt"))
    sample_names[file.exists(file.path(align_dir, sample_names, "QC_FAIL"))]
}))

out_table$succeed_alignment <- ifelse(
    out_table$sample_name %in% early_fails,
    "FAIL",
    "PASS"
)

# Find RNA align failures -------------------------------------------------

# GEX.Reads.mapped.confidently.to.transcriptome < 0.5
rna_reads_transcriptome_fails <- cellr_summary$Sample.ID[
    cellr_summary$GEX.Reads.mapped.confidently.to.transcriptome < 0.5
]

out_table$RNA_reads_map_confidently_to_transcriptome_0.5 <- ifelse(
    out_table$sample_name %in% rna_reads_transcriptome_fails,
    "FAIL",
    "PASS"
)

# GEX.Fraction.of.transcriptomic.reads.in.cells < 0.6
rna_reads_in_cells_fails <- cellr_summary$Sample.ID[
    cellr_summary$GEX.Fraction.of.transcriptomic.reads.in.cells < 0.6
]

out_table$RNA_transcriptomic_reads_in_cells_0.6 <- ifelse(
    out_table$sample_name %in% rna_reads_in_cells_fails,
    "FAIL",
    "PASS"
)

# GEX.Reads.with.TSO > 0.55
rna_reads_tso_fails <- cellr_summary$Sample.ID[
    cellr_summary$GEX.Reads.with.TSO > 0.55
]

out_table$RNA_transcriptomic_reads_tso_0.55 <- ifelse(
    out_table$sample_name %in% rna_reads_tso_fails,
    "FAIL",
    "PASS"
)


# Find ATAC align failures ------------------------------------------------

# ATAC.Unmapped.read.pairs > 0.1
atac_reads_unmapped_fails <- cellr_summary$Sample.ID[
    cellr_summary$ATAC.Unmapped.read.pairs > 0.1
]

out_table$ATAC_reads_unmapped_0.1 <- ifelse(
    out_table$sample_name %in% atac_reads_unmapped_fails,
    "FAIL",
    "PASS"
)

# ATAC.Fraction.of.high.quality.fragments.in.cells < 0.5
atac_frags_in_cells_fails <- cellr_summary$Sample.ID[
    cellr_summary$ATAC.Fraction.of.high.quality.fragments.in.cells < 0.5
]

out_table$ATAC_quality_fragments_in_cells_0.5 <- ifelse(
    out_table$sample_name %in% atac_frags_in_cells_fails,
    "FAIL",
    "PASS"
)

# ATAC.Fraction.of.transposition.events.in.peaks.in.cells < 0.2
atac_events_in_peaks_fails <- cellr_summary$Sample.ID[
    cellr_summary$ATAC.Fraction.of.transposition.events.in.peaks.in.cells < 0.2
]

out_table$ATAC_events_in_peaks_in_cells_0.2 <- ifelse(
    out_table$sample_name %in% atac_events_in_peaks_fails,
    "FAIL",
    "PASS"
)


# Add final columns for overall PASS/FAIL status --------------------------

rna_cols <- grep("RNA_", colnames(out_table))
out_table$RNA_all_qc_pass <- ifelse(
    apply(out_table[, rna_cols], 1, function(x) "FAIL" %in% x),
    "FAIL",
    "PASS"
)

atac_cols <- grep("ATAC_", colnames(out_table))
out_table$ATAC_all_qc_pass <- ifelse(
    apply(out_table[, atac_cols], 1, function(x) "FAIL" %in% x),
    "FAIL",
    "PASS"
)

out_table$all_qc_pass <- ifelse(
    apply(out_table, 1, function(x) "FAIL" %in% x),
    "FAIL",
    "PASS"
)


# Fix early fail columns --------------------------------------------------

# Currently samples that failed to align get a PASS for the other QC metrics.
# Will set them to NA, and then for the final QC columns, they will be fails.
for (cname in c(rna_cols, atac_cols)) {
    out_table[[cname]][out_table$succeed_alignment == "FAIL"] <- NA
}
for (cname in c("RNA_all_qc_pass", "ATAC_all_qc_pass", "all_qc_pass")) {
    out_table[[cname]][out_table$succeed_alignment == "FAIL"] <- "FAIL"
}


# Add barcoded_cell_sample_name -------------------------------------------

library_meta <- read.csv(path_library_meta)

bc_map <- setNames(
    library_meta$barcoded_cell_sample_name,
    library_meta$sample_name
)
out_table$barcoded_cell_sample_name <- unname(bc_map[out_table$sample_name])

# make it the first column
n_col <- ncol(out_table)
out_table <- out_table[, c(n_col, 1:(n_col-1))]


# Export table ------------------------------------------------------------

write.csv(out_table, out_path, row.names = FALSE)

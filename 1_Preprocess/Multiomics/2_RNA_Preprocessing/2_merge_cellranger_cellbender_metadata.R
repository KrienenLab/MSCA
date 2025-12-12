# This script will combine multiple sets of metadata that were generated and
# output by the pipeline.
suppressPackageStartupMessages(library(dplyr))
pdir <- '/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster'
path_out_cellr_summary <- file.path(pdir, "metadata/cellranger-arc_summary.csv")
path_out_cellb_metrics <- file.path(pdir, "metadata/cellbender_metrics.csv")

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
    sample_names <- readLines(file.path(align_dir, "sample_names.txt"))
    sample_names <- sample_names[!file.exists(
        file.path(align_dir, sample_names, "QC_FAIL")
    )]
    return(sample_names)
})
sample_names <- unlist(unname(sample_names))


# cellranger-arc summaries ------------------------------------------------

cellr_summary <- lapply(
    align_dirs, function(x) read.csv(file.path(x, "cellranger-arc_summaries.csv"))
)
cellr_summary <- do.call(rbind, unname(cellr_summary))
stopifnot(all(cellr_summary$Sample.ID == sample_names))

cat("Exporting merged cellranger-arc summaries...\n")
cellr_summary <- cellr_summary[order(cellr_summary$Sample.ID), ]
write.csv(cellr_summary, path_out_cellr_summary, row.names = FALSE)


# cellbender metrics ------------------------------------------------------

cellb_paths <- sapply(align_dirs, function(align_dir) {
    sample_names <- readLines(file.path(align_dir, "sample_names.txt"))
    sample_names <- sample_names[!file.exists(
        file.path(align_dir, sample_names, "QC_FAIL")
    )]
    cellb_paths <- sapply(sample_names, function(sample_name) {
        file.path(
            align_dir, sample_name, "cellbender_out/cellbender_output_metrics.csv"
        )
    })
    stopifnot(all(file.exists(cellb_paths)))
    cellb_paths
})
cellb_paths <- do.call(c, unname(cellb_paths))
cellb_metrics <- lapply(cellb_paths, read.csv, header = FALSE)
stopifnot(all(names(cellb_metrics) == sample_names))

# metrics file has row for each metric, 2 columns (name, value), and no header;
# want to turn those rows into columns, and add a column for sample_name, and
# then combine into a single dataframe
cellb_metrics <- Map(
    function(df, nm) {
        # turn column 1 into column names, and turn column 2 into a single row
        out <- as.data.frame(setNames(as.list(df[[2]]), df[[1]]))
        cbind(data.frame(sample_name = nm), out)
    }, 
    cellb_metrics, 
    names(cellb_metrics)
) %>% 
    unname %>% 
    do.call(rbind, .)

cat("Exporting merged cellbender metrics...\n")
cellb_metrics <- cellb_metrics[order(cellb_metrics$sample_name), ]
write.csv(cellb_metrics, path_out_cellb_metrics, row.names = FALSE)

suppressPackageStartupMessages({
    library(Seurat)
    library(scDblFinder)
})

pdir <- '/jukebox/krienen/marm_hmba_integration/250602_reprocess_and_recluster'
path_qc_table <- file.path(pdir, 'qc/library_qc_pass_fail_table.csv')
path_added_qc_fails <- file.path(pdir, 'qc/additional_qc_fails.csv')

qc_table <- read.csv(path_qc_table)
added_fails <- read.csv(path_added_qc_fails)
added_fails <- subset(added_fails, modality_fail == "RNA")

out_dir_curated_rna <- file.path(pdir, "data/rna")
out_path <- file.path(out_dir_curated_rna, "rna_merged.rds") # basically not important
out_path_curation_stats   <- file.path(pdir, "metadata/cell_curation_stats.csv")
out_path_cell_cycle_stats <- file.path(pdir, "metadata/cell_cycle_stats.csv")

# Cell curation scheme
nCount_RNA_min   <- 2e3
nCount_RNA_max   <- 1e5
nFeature_RNA_min <- 1e3
nFeature_RNA_max <- 1.3e4
percent_mt_max   <- 3

if (!dir.exists(file.path(pdir, "data")))
    dir.create(file.path(pdir, "data"))

if (!dir.exists(out_dir_curated_rna))
    dir.create(out_dir_curated_rna)


# Make mapping for barcoded_cell_sample_name ------------------------------

#
# Need to use the  barcoded_cell_sample_name instead of the sample_name for 
# cell barcodes
# 
sample_bc_map <- setNames(
    qc_table$barcoded_cell_sample_name,
    qc_table$sample_name
)

# Paths to cellbender h5 outputs ------------------------------------------

cat("Finding cellbender output h5 files...\n")
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

h5_paths <- sapply(align_dirs, function(align_dir) {
    sample_names <- readLines(file.path(align_dir, "sample_names.txt"))
    sample_names <- sample_names[!file.exists(
        file.path(align_dir, sample_names, "QC_FAIL")
    )]
    h5_paths <- sapply(sample_names, function(sample_name) {
        file.path(
            align_dir, sample_name, "cellbender_out",
            "cellbender_output_filtered_seurat.h5"
        )
    })
    stopifnot(all(file.exists(h5_paths)))
    h5_paths
})
h5_paths <- do.call(c, unname(h5_paths))
h5_paths <- data.frame(
    sample_name = names(h5_paths),
    h5_path = unname(h5_paths)
)

# will sort these now
h5_paths <- h5_paths[order(h5_paths$sample_name), ]


# Create Seurat objects from h5 files -------------------------------------

cat("Making Seurat objects from h5 files...\n")
seur_list <- lapply(h5_paths$h5_path, function(x) {
    # Multiome data will be a list of 2 matrices
    suppressMessages(suppressWarnings({
        SeuratObject::CreateSeuratObject(
            counts = Read10X_h5(x)[["Gene Expression"]]
        )
    }))
})
names(seur_list) <- h5_paths$sample_name


# Some basic Seurat cell annotations --------------------------------------

cat("Adding some basic metadata to Seurat objects...\n")
seur_list <- Map(
    function(seur, sample_name) {
        seur$orig.ident <- unname(sample_bc_map[sample_name])
        seur[["barcoded_cell_sample_name"]] <- unname(sample_bc_map[sample_name])
        seur[["sample_name"]] <- sample_name
        seur[["percent.mt"]]  <- PercentageFeatureSet(seur, pattern = "^(mt|MT)-")
        seur[["percent.Rpl"]] <- PercentageFeatureSet(seur, pattern = "^(Rpl|RPL)")
        seur[["percent.Hbb"]] <- PercentageFeatureSet(seur, pattern = "^(Hbb|HBB)")
        seur
    },
    seur_list,
    names(seur_list)
)


# Initialize cell curation stats ------------------------------------------

# Unlike before, will do these "stepwise", so numbers reflect progressive
# curation.
#
# Doublet detection is run before most other cell curation, except very 
# obviously empty droplets (<200 UMIs) are removed before doublet detection.
#
# Later will add on:
# - cells_over_min_features
# - cells_under_max_pct_mt
# - singlets 
cell_curation_table <- data.frame(
    sample_name = names(seur_list),
    cells_before_curation = sapply(seur_list, function(x) dim(x)[2]),
    cells_over_200_umi = sapply(
        seur_list, function(x) sum(x$nCount_RNA >= 200)
    ),
    cells_over_min_umi = sapply(
        seur_list, function(x) sum(x$nCount_RNA >= nCount_RNA_min)
    ),
    cells_under_max_umi = sapply(
        seur_list, function(x) sum(
            (x$nCount_RNA >= nCount_RNA_min) & 
                (x$nCount_RNA <= nCount_RNA_max)
        )
    )
)
rownames(cell_curation_table) <- NULL


# Doublet annotation ------------------------------------------------------

cat("Pre-filtering cells & running scDblFinder to annotate doublets...\n")

# Pre-filter empties
seur_list <- lapply(seur_list, subset, nCount_RNA >= 200)

# Run scDblFinder and add to Seurat metadata:
# - scDblFinder.class
# - scDblFinder.cxds_score
# - scDblFinder.score
# - scDblFinder.weighted
## Increasing from default dbr.per1k=0.008 (for 1% doublet rate per 1k cells),
seur_list <- lapply(seur_list, function(x) {
    suppressMessages(
        sce <- scDblFinder(x[["RNA"]]$counts, dbr.per1k = 0.01, verbose = FALSE)
    )
    AddMetaData(x, as.data.frame(colData(sce)))
})


# Cell curation -----------------------------------------------------------

cat("Filtering cells (including removing doublets...\n")
seur_list <- lapply(seur_list, subset, nCount_RNA >= nCount_RNA_min)
seur_list <- lapply(seur_list, subset, nCount_RNA <= nCount_RNA_max)

seur_list <- lapply(seur_list, subset, nFeature_RNA >= nFeature_RNA_min)
cell_curation_table$cells_over_min_features <- sapply(
    seur_list, function(x) dim(x)[2]
)

seur_list <- lapply(seur_list, subset, nFeature_RNA <= nFeature_RNA_max)
cell_curation_table$cells_under_max_features <- sapply(
    seur_list, function(x) dim(x)[2]
)

seur_list <- lapply(seur_list, subset, percent.mt < percent_mt_max)
cell_curation_table$cells_under_max_pct_mt <- sapply(
    seur_list, function(x) dim(x)[2]
)

seur_list <- lapply(seur_list, subset, scDblFinder.class == "singlet")
cell_curation_table$singlets <- sapply(
    seur_list, function(x) dim(x)[2]
)


# Rename columns in cell curation table -----------------------------------

idx <- which(names(cell_curation_table) == "cells_over_min_umi")
names(cell_curation_table)[idx] <- sub(
    "min", nCount_RNA_min, names(cell_curation_table)[idx]
)
idx <- which(names(cell_curation_table) == "cells_under_max_umi")
names(cell_curation_table)[idx] <- sub(
    "max", nCount_RNA_max, names(cell_curation_table)[idx]
)

idx <- which(names(cell_curation_table) == "cells_over_min_features")
names(cell_curation_table)[idx] <- sub(
    "min", nFeature_RNA_min, names(cell_curation_table)[idx]
)
idx <- which(names(cell_curation_table) == "cells_under_max_features")
names(cell_curation_table)[idx] <- sub(
    "max", nFeature_RNA_max, names(cell_curation_table)[idx]
)

idx <- which(names(cell_curation_table) == "cells_under_max_pct_mt")
names(cell_curation_table)[idx] <- sub(
    "max", percent_mt_max, names(cell_curation_table)[idx]
)


# Export cell curation table ----------------------------------------------

cat("Exporting cell curation statistics...\n")
write.csv(cell_curation_table, out_path_curation_stats, row.names = FALSE)


# Merge Seurat objects ----------------------------------------------------

cat('Modifying cell barcodes to have suffix of "`barcoded_cell_sample_name`"...\n')
seur_list <- lapply(
    seur_list, function(x) RenameCells(x, new.names = paste0(
        sub("1$", "", Cells(x)), x$barcoded_cell_sample_name
    ))
)

cat("Merging Seurat objects...\n")
suppressWarnings(
    seur <- merge(seur_list[[1]], seur_list[-1])
)

cat("Exporting (temporary) Seurat before joining layers...\n")
SeuratObject::SaveSeuratRds(
    seur, sub("rna_merged", "rna_merged_nojoin", out_path)
)

# If integrating, can re-split these, but cannot do cell cycle scoring when
# multiple layers present
cat("Joining layers (into chunks of max 80 libraries)...\n")

# Shockingly small limit on non-zero dgCMatrix elements (2^31-1 nz elements);
# so we'll have to keep the matrices separate in R, merge later in python
nlayers <- length(Layers(seur))
max_chunk <- 80
if (nlayers > max_chunk) {
    layer_chunks <- findInterval(seq_len(nlayers), seq(0, nlayers, max_chunk))
    seur <- lapply(unique(layer_chunks), function(i) {
        lyrs <- Layers(seur)[layer_chunks == i]
        cls <- Cells(seur, layer = lyrs)
        JoinLayers(subset(seur, cells = cls))
    })
    invisible(gc())
} else {
    seur <- list(JoinLayers(seur))
}

cat("Exporting (temporary) chunked Seurat with joined layers...\n")
saveRDS(
    seur, sub("rna_merged", "rna_merged_chunkjoin", out_path)
)


# Embedding ---------------------------------------------------------------

# cat("Running NormalizeData, FindVariableFeatures, ScaleData, RunPCA...\n")
cat("Running NormalizeData (for cell cycle annotations)...\n")
seur <- lapply(seur, NormalizeData, verbose = FALSE)


# Cell cycle scoring ------------------------------------------------------

cat("Generating cell cycle annotations...\n")

# Check whether to use human style gene names (as provided by Seurat) or 
# mouse style (cap-case)
toCapCase <- function(x) {
    .cap <- function(s) paste0(
        toupper(substr(s, 1, 1)),
        tolower(substr(s, 2, nchar(s)))
    )
    unname(sapply(x, .cap))
}
cc_genes <- lapply(cc.genes.updated.2019, toCapCase)

# check which genes to use
if (sum(unlist(cc_genes) %in% rownames(seur[[1]])) < 
    sum(unlist(cc.genes.updated.2019) %in% rownames(seur[[1]]))) {
    cc_genes <- cc.genes.updated.2019
}

cc_genes <- lapply(cc_genes, function(x) x[x %in% rownames(seur[[1]])])

seur <- lapply(
    seur,
    CellCycleScoring,
    s.features = cc_genes$s.genes, 
    g2m.features = cc_genes$g2m.genes, 
    set.ident = FALSE
)

cell_cycle_results <- do.call(rbind, lapply(seur, function(x) {
    aggregate(Phase ~ sample_name, x@meta.data, FUN = table)
}))


cat("Exporting cell cycle annotations by sample...\n")
write.csv(cell_cycle_results, out_path_cell_cycle_stats, row.names = FALSE)


# Export merged RNA Seurat  -----------------------------------------------

cat("Exporting Seurat RNA (`rna_merged.rds`, which contains RNA libraries failing QC)...\n")
seur <- lapply(seur, function(x) {
    x[["RNA"]]$data <- NULL
    x
})
saveRDS(seur, out_path)

## from when I messed up barcodes
# seur <- lapply(seur, function(x) {
#     RenameCells(x, new.names = sub("_[A-Z]*-1$", "", Cells(x)))
# })
# saveRDS(seur, out_path)


# Subset to only RNA QC Pass libraries ------------------------------------

cat("Subsetting to only libraries passing RNA library QC...\n")
rna_pass <- qc_table$sample_name[qc_table$RNA_all_qc_pass == "PASS"]
rna_pass <- setdiff(rna_pass, added_fails$sample_name)
seur <- lapply(seur, function(x) x[, x$sample_name %in% rna_pass])


# Export cell metadata & count matrices -----------------------------------

cat("Exporting cell metadata and gene metadata csv files......\n")
cell_meta <- do.call(rbind, lapply(seur, function(x) x@meta.data))
gene_meta <- as.data.frame(seur[[1]][["RNA"]]@meta.data)
rownames(gene_meta) <- rownames(seur[[1]])

write.csv(cell_meta, file.path(out_dir_curated_rna, "rna_cell_metadata.csv"))
write.csv(gene_meta, file.path(out_dir_curated_rna, "rna_gene_metadata.csv"))

cat("Transposing and exporting matrix market files (to be in cell-by-gene format)...\n")
for (i in seq_along(seur)) {
    out <- file.path(out_dir_curated_rna, paste0("rna_counts_", i, ".mtx"))
    Matrix::writeMM(t(seur[[i]][["RNA"]]["counts"]), out)
}

cat("Script finished.\n")

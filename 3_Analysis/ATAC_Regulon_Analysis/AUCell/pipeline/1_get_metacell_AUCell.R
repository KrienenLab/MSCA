# Score each eRegulon by enrichment of its target genes in RNA and 
# its cistrome peaks in ATAC within each cluster, using rank-based AUC
#
# main outputs:
# - score_rna_auc.csv  (eRegulon by cluster)
# - score_atac_auc.csv (eRegulon by cluster)
#
# also the R AUC objects are exported:
# - auc_rna.rds
# - auc_atac.rds
#
# I'm not sure what else is in those files but they should be useable in 
# downstream AUCell pipeline
suppressPackageStartupMessages({
    library(dplyr)
    library(AUCell)
})

repo <- here::here()  # root of this repository (run from anywhere inside it)
pdir   <- file.path(repo, "3_Analysis/ATAC_Regulon_Analysis")
outdir <- file.path(pdir, "AUCell/data")

if (!dir.exists(outdir))
    dir.create(outdir)


# Paths -------------------------------------------------------------------

# RNA & ATAC cluster metacell matrices
path_mat_rna <- file.path(pdir, "common_data/Cluster_v4_by_genes_matrix.csv")
path_mat_atac <- file.path(pdir, "common_data/Cluster_v4_by_462k_peaks_matrix.csv")


# Imports -----------------------------------------------------------------

# gene-by-cluster
mat_rna <- data.table::fread(path_mat_rna, header=TRUE)
rnames <- mat_rna[[1]]
mat_rna <- as.matrix(mat_rna[,-1])
rownames(mat_rna) <- rnames

# > dim(mat_rna)
# [1] 35787  1811

# peak-by-cluster
mat_atac <- data.table::fread(path_mat_atac, header=TRUE)
rnames <- mat_atac[[1]]
mat_atac <- as.matrix(mat_atac[,-1])
rownames(mat_atac) <- rnames

# > dim(mat_atac)
# [1] 462698   1810


# Import eRegulons --------------------------------------------------------

ereg_dir <- read.delim(
    file.path(repo, "3_Analysis/HPF_Gradient_Analysis/data/eRegulon_direct.tsv")
)
ereg_ext <- read.delim(
    "/jukebox/krienen/lab_members/victor/marm_hmba_atac_50k/scplus_neuron/outs/eRegulons_extended.tsv"
)
ereg <- rbind(ereg_dir, ereg_ext)



# Group genes/peaks by eRegulons ------------------------------------------

# Target genes by regulon
ereg_genes <- ereg %>%
    group_by(eRegulon_name) %>%
    summarise(genes = list(unique(Gene)), .groups = "drop")

# Peaks by Regulon (if column name differs, change 'region')
ereg_peaks <- ereg %>%
    group_by(eRegulon_name) %>%
    summarise(peaks = list(unique(Region)), .groups = "drop")

# And now join (eRegulon_name, genes, peaks)
ereg_both <- inner_join(ereg_genes, ereg_peaks, by = "eRegulon_name")
stopifnot(length(unique(ereg$eRegulon_name)) == nrow(ereg_both))


# RNA eRegulon activity via AUCell (rank-AUC) -----------------------------

stopifnot(all(ereg$Gene %in% rownames(mat_rna)))
geneSets <- setNames(ereg_both$genes, ereg_both$eRegulon_name)

# Build rankings per cluster
rna_rankings <- AUCell_buildRankings(
    mat_rna,
    plotStats = FALSE,
    nCores = 1
)

# Compute regulon-by-cluster AUC martrix for RNA
# - aucMaxRank : how "top-weighted"; top 5% by default
auc_rna <- AUCell_calcAUC(
    geneSets,
    rna_rankings,
    aucMaxRank = ceiling(0.05 * nrow(rna_rankings)) 
)
score_rna_auc <- as.matrix(getAUC(auc_rna))


### EXPORT
write.csv(score_rna_auc, file.path(outdir, "score_rna_auc.csv"))
saveRDS(auc_rna, file.path(outdir, "auc_rna.rds"))



# ATAC eRegulon activity via rank-AUC on peaks ----------------------------

stopifnot(all(ereg$Region %in% rownames(mat_atac)))
peakSets <- setNames(ereg_both$peaks, ereg_both$eRegulon_name)

atac_rankings <- AUCell_buildRankings(
    mat_atac,
    plotStats = FALSE,
    nCores = 1
)
auc_atac <- AUCell_calcAUC(
    peakSets,
    atac_rankings,
    aucMaxRank = ceiling(0.05 * nrow(atac_rankings)) 
)
score_atac_auc <- as.matrix(getAUC(auc_atac))

### EXPORT
write.csv(score_atac_auc, file.path(outdir, "score_atac_auc.csv"))
saveRDS(auc_atac, file.path(outdir, "auc_atac.rds"))

stopifnot(all(rownames(score_rna_auc) == rownames(score_atac_auc)))


# Script end --------------------------------------------------------------

cat("Script finished.\n")
quit(save="no")


# SCRAP: ATAC chromVAR ----------------------------------------------------

# a good idea for when we have time

# BiocManager::install(c("chromVAR","BSgenome", "Biostrings)) 
library(SummarizedExperiment)
library(chromVAR)
library(BiocParallel)
register(MulticoreParam(1))

peak_gr <- GRanges(rownames(mat_atac))

# SummarizedExperiment with *raw counts* recommended
se <- SummarizedExperiment(
    assays = list(counts = as.matrix(mat_atac)),
    rowRanges = peak_gr
)

# Add GC bias 
# se <- addGCBias(se, genome = <BSgenome object>)
## just need to see how the output is actually formatted and do with getSeq

# Build a binary peak × regulon incidence matrix
all_peaks <- rownames(mat_atac)
motif_ix <- sparseMatrix(
    i = unlist(lapply(seq_along(peakSets), function(k) match(peakSets[[k]], all_peaks))),
    j = rep(seq_along(peakSets), lengths(peakSets)),
    x = 1,
    dims = c(length(all_peaks), length(peakSets)),
    dimnames = list(all_peaks, names(peakSets))
)

# Deviations
dev <- computeDeviations(object = se, annotations = motif_ix)
dev_mat <- deviations(dev)              # regulon × cluster (z-like)
score_atac_dev <- dev_mat[rownames(dev_mat) %in% rownames(score_rna_auc), , drop = FALSE]

# Combine with RNA z-scores (rename to align if necessary)
common_regs2 <- intersect(rownames(score_rna_auc), rownames(score_atac_dev))
rna_z2  <- scale_rows(score_rna_auc[common_regs2, , drop = FALSE])
atac_z2 <- scale_rows(score_atac_dev[common_regs2, , drop = FALSE])
score_combined_chromVAR <- (rna_z2 + atac_z2) / 2



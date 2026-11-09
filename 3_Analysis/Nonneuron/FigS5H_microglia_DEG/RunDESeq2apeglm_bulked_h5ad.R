#!/usr/bin/env Rscript
# DESeq2 differential expression between two sets of groups of a pseudobulk / metacell h5ad, with apeglm
# log fold change shrinkage.
#
# Usage
#   Rscript RunDESeq2apeglm_bulked_h5ad.R <pseudobulk.h5ad> <group_column_in_obs> "<numerator_groups>" \
#       "<denominator_groups>" <out_dir> [counts_layer]
#
#   <numerator_groups>, <denominator_groups>  comma-separated values of obs[[group_column_in_obs]]; samples in
#                                             neither list are removed. Positive log2FoldChange = higher in numerator.
#   [counts_layer]                            layer with raw integer counts (default: X)
#
# Requires R packages DESeq2, apeglm, reticulate, and a Python environment with anndata (reticulate uses the active
# environment or RETICULATE_PYTHON).
#
# Outputs (out_dir)
#   DESeq2_apeglm_results.csv   DESeq2 results + log2FoldChange_shr / lfcSE_shr (apeglm-shrunken)
#   DESeqDataSet.rds            fitted DESeqDataSet
#   sample_group_mapping.csv    sample -> original group -> numerator / denominator

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5 || length(args) > 6)
  stop("Usage: Rscript RunDESeq2apeglm_bulked_h5ad.R <pseudobulk.h5ad> <group_column_in_obs> ",
       "\"<numerator_groups>\" \"<denominator_groups>\" <out_dir> [counts_layer]")

path_h5ad     <- args[1]
group_col     <- args[2]
numer_groups  <- trimws(unlist(strsplit(args[3], ",")))
denom_groups  <- trimws(unlist(strsplit(args[4], ",")))
out_dir       <- args[5]
counts_layer  <- if (length(args) == 6) args[6] else "X"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

### ---------------------------------------------------------------------------
### Read the h5ad
### ---------------------------------------------------------------------------
cat("Importing h5ad...\n")
adata <- reticulate::import("anndata", convert = FALSE)$read_h5ad(path_h5ad)

obs_df    <- reticulate::py_to_r(adata$obs)
obs_names <- as.character(reticulate::py_to_r(adata$obs_names$tolist()))
rownames(obs_df) <- obs_names
if (!(group_col %in% colnames(obs_df)))
  stop("Group column '", group_col, "' not found in .obs.")

var_names <- as.character(reticulate::py_to_r(adata$var_names$tolist()))
if (any(is.na(var_names) | var_names == "")) stop("Feature names (var_names) contain empty or NA values.")
if (any(duplicated(var_names))) stop("Feature names (var_names) contain duplicates.")

### ---------------------------------------------------------------------------
### Sample table: numerator vs denominator, other samples removed
### ---------------------------------------------------------------------------
coldat <- data.frame(sample_name = obs_names,
                     original_group = as.character(obs_df[[group_col]]),
                     row.names = obs_names, stringsAsFactors = FALSE)

for (side in list(list("Numerator", numer_groups), list("Denominator", denom_groups))) {
  missing <- setdiff(side[[2]], coldat$original_group)
  if (length(missing) > 0)
    stop(side[[1]], " groups not found in obs[['", group_col, "']]: ", paste(missing, collapse = ", "))
}

coldat$contrast_group <- ifelse(coldat$original_group %in% numer_groups, "numer",
                                ifelse(coldat$original_group %in% denom_groups, "denom", NA))
coldat <- coldat[!is.na(coldat$contrast_group), , drop = FALSE]
coldat$contrast_group <- factor(coldat$contrast_group, levels = c("denom", "numer"))  # denominator = reference

### ---------------------------------------------------------------------------
### Count matrix (genes x samples)
### ---------------------------------------------------------------------------
cat("Getting count matrix...\n")
counts_py <- if (counts_layer == "X") adata$X else adata$layers$get(counts_layer)
if (is.null(counts_py) || reticulate::py_is_none(counts_py))
  stop("Layer '", counts_layer, "' not found in adata$layers.")
if (reticulate::py_has_attr(counts_py, "toarray")) counts_py <- counts_py$toarray()

mat <- t(as.matrix(reticulate::py_to_r(counts_py)))
rownames(mat) <- var_names
colnames(mat) <- obs_names
mat <- mat[, rownames(coldat), drop = FALSE]

if (!is.integer(mat)) {
  if (any(abs(mat - round(mat)) > .Machine$double.eps^0.5))
    stop("The matrix looks normalized (non-integer); DESeq2 needs raw counts.")
  mat <- round(mat)
  storage.mode(mat) <- "integer"
}

### ---------------------------------------------------------------------------
### DESeq2 + apeglm
### ---------------------------------------------------------------------------
cat("Running DESeq2...\n")
dds <- DESeqDataSetFromMatrix(countData = mat, colData = coldat, design = ~ contrast_group)
dds <- DESeq(dds)
saveRDS(dds, file.path(out_dir, "DESeqDataSet.rds"))

coef <- "contrast_group_numer_vs_denom"
res  <- results(dds, name = coef)

cat("Running apeglm...\n")
shr <- lfcShrink(dds, coef = coef, res = res, type = "apeglm")

out <- as.data.frame(res)
out$log2FoldChange_shr <- shr$log2FoldChange
out$lfcSE_shr          <- shr$lfcSE
write.csv(out, file.path(out_dir, "DESeq2_apeglm_results.csv"), row.names = TRUE)

write.csv(data.frame(sample = rownames(coldat), original_group = coldat$original_group,
                     contrast_group = as.character(coldat$contrast_group)),
          file.path(out_dir, "sample_group_mapping.csv"), row.names = FALSE)

cat("Finished.\n")

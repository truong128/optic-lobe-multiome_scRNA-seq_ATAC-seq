#!/usr/bin/env Rscript
# Within-library sign check.
# Cell-level Wilcoxon p-values on the atlas are not written. Nuclei in one
# genotype are not replicates, so those p-values are not a result.
user_lib <- path.expand("~/R/library")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
root <- Sys.getenv("PROJECT_ROOT", unset = "")
if (!nzchar(root)) {
  root <- if (length(file_arg)) normalizePath(file.path(dirname(file_arg[1]), "..")) else getwd()
}
out <- file.path(root, "results")
genes <- c("nAChRalpha6", "CG44153", "Tl", "Sulf1", "Drl-2", "kek6")

seu <- readRDS(file.path(root, "results", "dm8_subset.rds"))
md <- seu@meta.data
rna <- GetAssayData(seu, assay = "RNA", layer = "counts")
atac <- GetAssayData(seu, assay = "ATAC", layer = "counts")
beat <- as.numeric(rna["beat-IIIc", ])
dip <- as.numeric(rna["CG14521", ])
subtype <- ifelse(beat > 0 & dip == 0, "pDm8", ifelse(dip > 0 & beat == 0, "yDm8", "other"))
meta <- data.frame(cell = colnames(seu), stage = as.character(md$stage),
                   library = as.character(md$orig.ident), subtype = subtype,
                   stringsAsFactors = FALSE)

focus_rna <- c("beat-IIIc", "CG14521", genes)
peaks <- c("chrX-5543757-5544256", "chrX-5544476-5544975", "chrX-5546437-5546936",
           "chrX-5548482-5548981", "chrX-5549283-5549782", "chrX-5551716-5552215")

lib_rows <- list()
for (assay_name in c("RNA", "ATAC")) {
  counts <- if (assay_name == "RNA") rna else atac
  feats <- if (assay_name == "RNA") focus_rna else peaks
  feats <- feats[feats %in% rownames(counts)]
  for (st in c("P24", "P48")) {
    libs <- unique(meta$library[meta$stage == st])
    for (lib in libs) {
      cells_p <- meta$cell[meta$stage == st & meta$library == lib & meta$subtype == "pDm8"]
      cells_y <- meta$cell[meta$stage == st & meta$library == lib & meta$subtype == "yDm8"]
      if (length(cells_p) < 5 || length(cells_y) < 5) next
      mp <- Matrix::rowMeans(counts[feats, cells_p, drop = FALSE])
      my <- Matrix::rowMeans(counts[feats, cells_y, drop = FALSE])
      logfc <- log2(mp + 0.1) - log2(my + 0.1)
      lib_rows[[paste(assay_name, st, lib)]] <- data.frame(
        assay = assay_name, stage = st, library = lib,
        n_pDm8 = length(cells_p), n_yDm8 = length(cells_y),
        feature = feats, log2FC = as.numeric(logfc),
        stringsAsFactors = FALSE
      )
    }
  }
}
lib_tab <- rbindlist(lib_rows)
fwrite(lib_tab, file.path(out, "library_sign.csv"))
lib_sum <- lib_tab[, .(
  n_libraries = .N,
  n_same_sign = sum(sign(log2FC) == sign(log2FC[1]) & log2FC[1] != 0),
  n_positive = sum(log2FC > 0)
), by = .(assay, stage, feature)]
fwrite(lib_sum, file.path(out, "library_sign_summary.csv"))
message("Library sign summary")
print(lib_sum)
message("DONE")

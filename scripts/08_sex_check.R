#!/usr/bin/env Rscript
# Where does the sex label come from? Compare it with sex-chromosome genes.
user_lib <- path.expand("~/R/library")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
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
seu <- readRDS(file.path(root, "results", "dm8_subset.rds"))
md <- seu@meta.data
rna <- GetAssayData(seu, assay = "RNA", layer = "counts")
labels <- fread(file.path(root, "results", "cell_labels.csv"))
labels <- labels[cell %in% colnames(rna)]

sex_genes <- c("Sxl", "tra", "dsx", "fru", "roX1", "roX2", "msl-2", "lncRNA:roX1", "lncRNA:roX2")
sex_genes <- sex_genes[sex_genes %in% rownames(rna)]
message("Sex genes present: ", paste(sex_genes, collapse = ", "))

genes <- fread(file.path(root, "ref", "genes.tsv"))
xgenes <- unique(genes[chrom == "X", gene_name])
xgenes <- intersect(xgenes, rownames(rna))
agenes <- unique(genes[chrom %in% c("2L", "2R", "3L", "3R"), gene_name])
agenes <- intersect(agenes, rownames(rna))

sub <- labels[stage %in% c("P24", "P48") & subtype %in% c("pDm8", "yDm8")]
rows <- list()
for (st in c("P24", "P48")) {
  for (sx in unique(sub[stage == st, sex])) {
    cells <- sub[stage == st & sex == sx, cell]
    if (length(cells) < 5) next
    totals <- Matrix::colSums(rna[, cells, drop = FALSE])
    xshare <- Matrix::colSums(rna[xgenes, cells, drop = FALSE]) / pmax(totals, 1)
    ashare <- Matrix::colSums(rna[agenes, cells, drop = FALSE]) / pmax(totals, 1)
    rec <- data.frame(
      stage = st, sex = sx, n_cells = length(cells),
      median_x_fraction = median(xshare),
      median_autosome_fraction = median(ashare),
      stringsAsFactors = FALSE
    )
    for (g in sex_genes) {
      v <- as.numeric(rna[g, cells])
      rec[[paste0(g, "_pct")]] <- mean(v > 0)
      rec[[paste0(g, "_mean")]] <- mean(v)
    }
    rows[[paste(st, sx)]] <- rec
  }
}
tab <- rbindlist(rows, fill = TRUE)
fwrite(tab, file.path(out, "sex_gene_check.csv"))
print(tab)
message("DONE")

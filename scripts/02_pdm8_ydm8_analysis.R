#!/usr/bin/env Rscript
# pDm8 vs yDm8 analysis inside author-annotated Dm8 nuclei (GSE305940).
user_lib <- path.expand("~/R/library")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(limma)
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(readxl)
  library(data.table)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
root <- Sys.getenv("PROJECT_ROOT", unset = "")
if (!nzchar(root)) {
  root <- if (length(file_arg)) normalizePath(file.path(dirname(file_arg[1]), "..")) else getwd()
}
fig <- file.path(root, "figures")
res <- file.path(root, "results")
dir.create(fig, showWarnings = FALSE, recursive = TRUE)
dir.create(res, showWarnings = FALSE, recursive = TRUE)

log_con <- file(file.path(res, "02_analysis_log.txt"), open = "wt")
sink(log_con, split = TRUE)
on.exit({ sink(); close(log_con) }, add = TRUE)
message("START ", Sys.time())

seu <- readRDS(file.path(res, "dm8_subset.rds"))
md <- seu@meta.data
stopifnot(all(md$AnnotatedJun26 == "Dm8"))

rna_counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
atac_counts <- GetAssayData(seu, assay = "ATAC", layer = "counts")
message("RNA ", nrow(rna_counts), " x ", ncol(rna_counts),
        " ATAC ", nrow(atac_counts), " x ", ncol(atac_counts))

beat <- as.numeric(rna_counts["beat-IIIc", ])
dip <- as.numeric(rna_counts["CG14521", ])
names(beat) <- colnames(seu)
names(dip) <- colnames(seu)

subtype <- ifelse(beat > 0 & dip == 0, "pDm8",
           ifelse(dip > 0 & beat == 0, "yDm8",
           ifelse(beat > 0 & dip > 0, "double_positive", "double_negative")))
subtype <- factor(subtype, levels = c("pDm8", "yDm8", "double_positive", "double_negative"))

labels <- data.frame(
  cell = colnames(seu),
  stage = md$stage,
  sex = md$sex,
  dgrp = md$dgrp,
  dgrp_new = md$dgrp.new,
  orig.ident = md$orig.ident,
  cluster_numeric = md$Clusters,
  subtype = as.character(subtype),
  beat_IIIc_count = beat,
  DIP_gamma_count = dip,
  nCount_RNA = md$nCount_RNA,
  nFeature_RNA = md$nFeature_RNA,
  percent.mt = md$percent.mt,
  nCount_ATAC = md$nCount_ATAC,
  nFeature_ATAC = md$nFeature_ATAC,
  TSS.enrichment = md$TSS.enrichment,
  nucleosome_signal = md$nucleosome_signal,
  stringsAsFactors = FALSE
)
fwrite(labels, file.path(res, "cell_labels.csv"))

assign_tab <- as.data.frame(table(stage = labels$stage, subtype = labels$subtype))
fwrite(assign_tab, file.path(res, "assignment_by_stage.csv"))
message("Assignment by stage:")
print(xtabs(~ stage + subtype, labels))

# mutual exclusivity among cells that express at least one marker
excl <- labels %>%
  group_by(stage) %>%
  summarise(
    n = n(),
    n_any_marker = sum(subtype %in% c("pDm8", "yDm8", "double_positive")),
    n_single = sum(subtype %in% c("pDm8", "yDm8")),
    frac_exclusive_among_marker_positive = ifelse(n_any_marker == 0, NA, n_single / n_any_marker),
    frac_p = mean(subtype == "pDm8"),
    frac_y = mean(subtype == "yDm8"),
    frac_double = mean(subtype == "double_positive"),
    frac_none = mean(subtype == "double_negative"),
    .groups = "drop"
  )
fwrite(excl, file.path(res, "marker_exclusivity_by_stage.csv"))
print(excl)

# QC by subtype x stage
qc_vars <- c("nCount_RNA", "nFeature_RNA", "percent.mt", "nCount_ATAC", "nFeature_ATAC", "TSS.enrichment", "nucleosome_signal")
qc_rows <- list()
for (st in unique(labels$stage)) {
  for (v in qc_vars) {
    for (g in c("pDm8", "yDm8", "double_positive", "double_negative")) {
      x <- labels[[v]][labels$stage == st & labels$subtype == g]
      qc_rows[[length(qc_rows) + 1]] <- data.frame(
        stage = st, subtype = g, metric = v, n = length(x),
        median = if (length(x)) median(x, na.rm = TRUE) else NA,
        q25 = if (length(x)) quantile(x, 0.25, na.rm = TRUE) else NA,
        q75 = if (length(x)) quantile(x, 0.75, na.rm = TRUE) else NA
      )
    }
    xp <- labels[[v]][labels$stage == st & labels$subtype == "pDm8"]
    xy <- labels[[v]][labels$stage == st & labels$subtype == "yDm8"]
    pval <- if (length(xp) >= 5 && length(xy) >= 5) suppressWarnings(wilcox.test(xp, xy)$p.value) else NA
    qc_rows[[length(qc_rows) + 1]] <- data.frame(
      stage = st, subtype = "wilcox_p_vs_y", metric = v, n = length(xp) + length(xy),
      median = pval, q25 = NA, q75 = NA
    )
  }
}
qc_df <- bind_rows(qc_rows)
fwrite(qc_df, file.path(res, "qc_by_subtype.csv"))

# ---------- plots: coexpression ----------
plot_df <- labels
plot_df$log_beat <- log1p(plot_df$beat_IIIc_count)
plot_df$log_dip <- log1p(plot_df$DIP_gamma_count)
plot_df$stage <- factor(plot_df$stage, levels = c("P0", "P24", "P48", "Adult"))
cols <- c(pDm8 = "#0072B2", yDm8 = "#E69F00", double_positive = "#009E73", double_negative = "#BBBBBB")

p_sc <- ggplot(plot_df, aes(log_beat, log_dip, color = subtype)) +
  geom_jitter(width = 0.08, height = 0.08, size = 0.7, alpha = 0.75) +
  facet_wrap(~ stage, nrow = 1) +
  scale_color_manual(values = cols, name = "Label") +
  labs(
    title = "beat-IIIc and DIP-gamma counts inside author-annotated Dm8 nuclei",
    subtitle = "GSE305940 multiome. Labels use raw counts: single-positive cells are confident subtypes.",
    x = "log1p(beat-IIIc count)", y = "log1p(DIP-gamma / CG14521 count)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(fig, "01_marker_coexpression.png"), p_sc, width = 11, height = 4.2, dpi = 160)

p_bar <- ggplot(assign_tab, aes(stage, Freq, fill = subtype)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = cols) +
  scale_x_discrete(limits = c("P0", "P24", "P48", "Adult")) +
  labs(
    title = "Dm8 subtype labels by developmental stage",
    subtitle = "pDm8 = beat-IIIc only; yDm8 = DIP-gamma only. Double-positive and double-negative cells are not used in differential tests.",
    x = "Stage", y = "Nuclei", fill = "Label"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(fig, "02_assignment_counts.png"), p_bar, width = 7.2, height = 4.4, dpi = 160)

# ---------- per-stage UMAP on RNA ----------
set.seed(1)
umap_list <- list()
for (st in c("P0", "P24", "P48", "Adult")) {
  cells <- rownames(md)[md$stage == st]
  if (length(cells) < 50) next
  sub <- subset(seu, cells = cells)
  DefaultAssay(sub) <- "RNA"
  sub <- NormalizeData(sub, verbose = FALSE)
  sub <- FindVariableFeatures(sub, nfeatures = 2000, verbose = FALSE)
  sub <- ScaleData(sub, verbose = FALSE)
  sub <- RunPCA(sub, npcs = 20, verbose = FALSE)
  sub <- RunUMAP(sub, dims = 1:15, verbose = FALSE, seed.use = 1)
  emb <- as.data.frame(Embeddings(sub, "umap"))
  emb$cell <- rownames(emb)
  emb$stage <- st
  emb$subtype <- labels$subtype[match(emb$cell, labels$cell)]
  emb$beat <- log1p(labels$beat_IIIc_count[match(emb$cell, labels$cell)])
  emb$dip <- log1p(labels$DIP_gamma_count[match(emb$cell, labels$cell)])
  umap_list[[st]] <- emb
  rm(sub); gc()
}
umap_df <- bind_rows(umap_list)
fwrite(umap_df, file.path(res, "dm8_umap_coordinates.csv"))

p_umap <- ggplot(umap_df, aes(umap_1, umap_2, color = subtype)) +
  geom_point(size = 0.7, alpha = 0.8) +
  facet_wrap(~ stage, nrow = 1) +
  scale_color_manual(values = cols) +
  labs(
    title = "RNA UMAP of Dm8 nuclei, computed separately at each stage",
    subtitle = "Dimensionality reduction is descriptive. Subtype color is the marker label, not a cluster cut.",
    x = "UMAP 1", y = "UMAP 2", color = "Label"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(fig, "03_umap_subtype.png"), p_umap, width = 11, height = 4.2, dpi = 160)

umap_long <- rbind(
  data.frame(umap_df[, c("umap_1", "umap_2", "stage")], gene = "beat-IIIc", value = umap_df$beat),
  data.frame(umap_df[, c("umap_1", "umap_2", "stage")], gene = "DIP-gamma", value = umap_df$dip)
)
p_feat <- ggplot(umap_long, aes(umap_1, umap_2, color = value)) +
  geom_point(size = 0.55) +
  facet_grid(gene ~ stage) +
  scale_color_viridis_c(option = "magma", name = "log1p(count)") +
  labs(
    title = "Marker expression on the within-stage Dm8 RNA UMAP",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme_bw(base_size = 10)
ggsave(file.path(fig, "04_umap_markers.png"), p_feat, width = 11, height = 5.6, dpi = 160)

# ---------- pseudobulk limma ----------
pseudobulk_limma <- function(counts, cell_ids, stage_labels, subtype_labels, genotype_labels,
                             min_cells = 5, min_genotypes = 4, min_pct = 0.10) {
  use <- subtype_labels %in% c("pDm8", "yDm8")
  cell_ids <- cell_ids[use]
  stage_labels <- stage_labels[use]
  subtype_labels <- subtype_labels[use]
  genotype_labels <- genotype_labels[use]
  keep_g <- !(genotype_labels %in% c("0", "NA", "", "unknown")) & !is.na(genotype_labels)
  cell_ids <- cell_ids[keep_g]
  stage_labels <- stage_labels[keep_g]
  subtype_labels <- subtype_labels[keep_g]
  genotype_labels <- genotype_labels[keep_g]

  out_tabs <- list()
  out_samples <- list()
  for (st in unique(stage_labels)) {
    message("Pseudobulk ", deparse(substitute(counts)), " ", st)
    ii <- which(stage_labels == st)
    if (length(ii) < 20) {
      message("  skip, too few labeled cells")
      next
    }
    genos <- sort(unique(genotype_labels[ii]))
    good_genos <- c()
    for (g in genos) {
      n_p <- sum(genotype_labels[ii] == g & subtype_labels[ii] == "pDm8")
      n_y <- sum(genotype_labels[ii] == g & subtype_labels[ii] == "yDm8")
      if (n_p >= min_cells && n_y >= min_cells) good_genos <- c(good_genos, g)
    }
    message("  genotypes passing >=", min_cells, " cells each subtype: ", length(good_genos),
            " / ", length(genos))
    if (length(good_genos) < min_genotypes) {
      message("  skip stage, fewer than ", min_genotypes, " paired genotypes")
      next
    }
    # build pseudobulk
    pb_cols <- list()
    sample_rows <- list()
    for (g in good_genos) {
      for (sbt in c("pDm8", "yDm8")) {
        cells <- cell_ids[ii][genotype_labels[ii] == g & subtype_labels[ii] == sbt]
        pb_cols[[paste(st, g, sbt, sep = "|")]] <- Matrix::rowSums(counts[, cells, drop = FALSE])
        sample_rows[[length(sample_rows) + 1]] <- data.frame(
          sample = paste(st, g, sbt, sep = "|"),
          stage = st, genotype = g, subtype = sbt, n_cells = length(cells),
          stringsAsFactors = FALSE
        )
      }
    }
    pb <- do.call(cbind, pb_cols)
    samples <- bind_rows(sample_rows)
    samples$libsize <- colSums(pb)
    # detection filter on the single cells that entered pseudobulks
    cells_p <- cell_ids[ii][genotype_labels[ii] %in% good_genos & subtype_labels[ii] == "pDm8"]
    cells_y <- cell_ids[ii][genotype_labels[ii] %in% good_genos & subtype_labels[ii] == "yDm8"]
    pct_p <- Matrix::rowMeans(counts[, cells_p, drop = FALSE] > 0)
    pct_y <- Matrix::rowMeans(counts[, cells_y, drop = FALSE] > 0)
    keep_feat <- (pct_p >= min_pct) | (pct_y >= min_pct)
    message("  features tested ", sum(keep_feat), " of ", length(keep_feat),
            " | cells p ", length(cells_p), " y ", length(cells_y))
    pb_f <- pb[keep_feat, , drop = FALSE]
    y <- DGEList(counts = pb_f, samples = samples)
    y <- calcNormFactors(y, method = "TMM")
    logcpm <- cpm(y, log = TRUE, prior.count = 1)
    subtype_f <- factor(samples$subtype, levels = c("yDm8", "pDm8"))
    genotype_f <- factor(samples$genotype)
    design <- model.matrix(~ genotype_f + subtype_f)
    fit <- lmFit(logcpm, design)
    fit <- eBayes(fit, trend = TRUE, robust = TRUE)
    coef <- grep("^subtype_f", colnames(design), value = TRUE)
    tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
    tt$feature <- rownames(tt)
    tt$stage <- st
    tt$pct_pDm8 <- pct_p[keep_feat]
    tt$pct_yDm8 <- pct_y[keep_feat]
    tt$n_pDm8 <- length(cells_p)
    tt$n_yDm8 <- length(cells_y)
    tt$n_genotypes <- length(good_genos)
    # direction consistency across genotypes
    # logcpm columns follow samples order
    p_idx <- which(samples$subtype == "pDm8")
    y_idx <- which(samples$subtype == "yDm8")
    # match genotypes
    p_gen <- samples$genotype[p_idx]
    y_gen <- samples$genotype[y_idx]
    ord <- match(p_gen, y_gen)
    delta <- logcpm[, p_idx, drop = FALSE] - logcpm[, y_idx[ord], drop = FALSE]
    mean_delta <- rowMeans(delta)
    tt$frac_genotypes_same_sign <- rowMeans(sign(delta) == sign(mean_delta))
    tt$mean_paired_log2FC <- mean_delta
    tt <- tt[order(tt$adj.P.Val, -abs(tt$logFC)), ]
    out_tabs[[st]] <- tt
    samples$assay_stage <- st
    out_samples[[st]] <- samples
    message("  FDR<0.05 ", sum(tt$adj.P.Val < 0.05),
            " and |logFC|>=0.5 ", sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.5))
  }
  list(tables = out_tabs, samples = out_samples)
}

message("\n===== RNA pseudobulk =====")
rna_res <- pseudobulk_limma(
  counts = rna_counts,
  cell_ids = labels$cell,
  stage_labels = labels$stage,
  subtype_labels = labels$subtype,
  genotype_labels = labels$dgrp_new,
  min_cells = 5,
  min_genotypes = 4,
  min_pct = 0.10
)

message("\n===== ATAC pseudobulk =====")
atac_res <- pseudobulk_limma(
  counts = atac_counts,
  cell_ids = labels$cell,
  stage_labels = labels$stage,
  subtype_labels = labels$subtype,
  genotype_labels = labels$dgrp_new,
  min_cells = 5,
  min_genotypes = 4,
  min_pct = 0.05
)

write_de <- function(tabs, path) {
  if (!length(tabs)) return(invisible(NULL))
  all <- bind_rows(lapply(tabs, function(tt) {
    tt$feature <- rownames(tt)
    tt
  }))
  # rownames were already copied to feature before bind; topTable rownames survive in each tt
  fwrite(all, path)
  all
}

# topTable data.frames already have feature column; bind_rows keeps it
rna_all <- if (length(rna_res$tables)) bind_rows(rna_res$tables) else NULL
atac_all <- if (length(atac_res$tables)) bind_rows(atac_res$tables) else NULL
if (!is.null(rna_all)) fwrite(rna_all, file.path(res, "de_rna_pseudobulk.csv"))
if (!is.null(atac_all)) fwrite(atac_all, file.path(res, "de_atac_pseudobulk.csv"))
if (length(rna_res$samples)) fwrite(bind_rows(rna_res$samples), file.path(res, "pseudobulk_samples_rna.csv"))
if (length(atac_res$samples)) fwrite(bind_rows(atac_res$samples), file.path(res, "pseudobulk_samples_atac.csv"))

# sanity: classifier genes
if (!is.null(rna_all)) {
  message("Classifier genes in RNA results:")
  print(rna_all[rna_all$feature %in% c("beat-IIIc", "CG14521"),
                c("feature", "stage", "logFC", "P.Value", "adj.P.Val", "pct_pDm8", "pct_yDm8")])
}

# volcano plots
volcano_plot <- function(tt, title, label_n = 12, outfile) {
  if (is.null(tt) || !nrow(tt)) return(invisible(NULL))
  tt$neglog10 <- -log10(pmax(tt$P.Value, 1e-300))
  tt$sig <- tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.5
  # labels: controls plus top by p
  controls <- tt$feature %in% c("beat-IIIc", "CG14521")
  top <- tt[order(tt$adj.P.Val, -abs(tt$logFC)), ]
  top <- top[!(top$feature %in% c("beat-IIIc", "CG14521")), ]
  lab_feat <- unique(c(tt$feature[controls], head(top$feature, label_n)))
  tt$label <- ifelse(tt$feature %in% lab_feat, tt$feature, "")
  p <- ggplot(tt, aes(logFC, neglog10, color = sig)) +
    geom_point(size = 0.7, alpha = 0.7) +
    geom_text_repel(aes(label = label), size = 2.6, max.overlaps = 30, color = "black", show.legend = FALSE) +
    scale_color_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#B0B0B0"),
                       labels = c(`TRUE` = "FDR<0.05 and |log2FC|>=0.5", `FALSE` = "not significant"),
                       name = NULL) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, linewidth = 0.3) +
    labs(
      title = title,
      subtitle = "Positive log2 fold-change means higher in pDm8 than yDm8. Limma on TMM-normalized pseudobulk, paired by DGRP group.",
      x = "log2 fold-change (pDm8 / yDm8)",
      y = "-log10(p-value)"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(outfile, p, width = 7.4, height = 5.4, dpi = 160)
}

for (st in names(rna_res$tables)) {
  volcano_plot(rna_res$tables[[st]],
               paste0("RNA pseudobulk: pDm8 vs yDm8 at ", st),
               outfile = file.path(fig, paste0("05_volcano_rna_", st, ".png")))
}
for (st in names(atac_res$tables)) {
  volcano_plot(atac_res$tables[[st]],
               paste0("ATAC pseudobulk: pDm8 vs yDm8 at ", st),
               label_n = 8,
               outfile = file.path(fig, paste0("06_volcano_atac_", st, ".png")))
}

# heatmap of top RNA DEGs (exclude classifier genes from the discovered set, but show them)
if (!is.null(rna_all) && nrow(rna_all)) {
  heat_genes <- c()
  for (st in names(rna_res$tables)) {
    tt <- rna_res$tables[[st]]
    tt <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) >= 0.5 & !(tt$feature %in% c("beat-IIIc", "CG14521")), ]
    tt_p <- head(tt[order(-tt$logFC), ], 10)
    tt_y <- head(tt[order(tt$logFC), ], 10)
    heat_genes <- unique(c(heat_genes, tt_p$feature, tt_y$feature))
  }
  heat_genes <- unique(c("beat-IIIc", "CG14521", heat_genes))
  heat_genes <- heat_genes[heat_genes %in% rownames(rna_counts)]
  if (length(heat_genes) >= 2 && length(rna_res$samples)) {
    # average log1p CPM-like per subtype x stage using single cells in the tested sets
    heat_rows <- list()
    for (st in names(rna_res$tables)) {
      for (sbt in c("pDm8", "yDm8")) {
        cells <- labels$cell[labels$stage == st & labels$subtype == sbt]
        if (length(cells) < 5) next
        m <- Matrix::rowMeans(rna_counts[heat_genes, cells, drop = FALSE])
        heat_rows[[paste(st, sbt)]] <- log1p(m * 1e4 / mean(Matrix::colSums(rna_counts[, cells, drop = FALSE])))
      }
    }
    heat_mat <- do.call(cbind, heat_rows)
    heat_mat <- heat_mat[apply(heat_mat, 1, sd) > 0, , drop = FALSE]
    heat_z <- t(scale(t(heat_mat)))
    heat_z[heat_z > 2.5] <- 2.5
    heat_z[heat_z < -2.5] <- -2.5
    hm <- as.data.frame(as.table(heat_z))
    colnames(hm) <- c("gene", "sample", "z")
    hm$gene <- factor(hm$gene, levels = rev(rownames(heat_z)))
    p_hm <- ggplot(hm, aes(sample, gene, fill = z)) +
      geom_tile() +
      scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#E69F00", name = "row z") +
      labs(
        title = "Top subtype-differential genes in Dm8",
        subtitle = "Mean expression per subtype and stage, row-scaled. beat-IIIc and CG14521 are the labels used to define the groups.",
        x = NULL, y = NULL
      ) +
      theme_bw(base_size = 9) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(fig, "07_heatmap_top_genes.png"), p_hm, width = 8.2, height = max(5, 0.22 * nrow(heat_z) + 2), dpi = 160)
    fwrite(as.data.frame(heat_mat), file.path(res, "heatmap_mean_expression.csv"), row.names = TRUE)
  }
}

# ---------- gene coordinates and peak annotation ----------
genes <- fread(file.path(root, "ref", "genes.tsv"))
genes[, start := as.integer(start)]
genes[, end := as.integer(end)]
genes[, chrom := as.character(chrom)]
genes[, chrom := ifelse(startsWith(chrom, "chr"), chrom, paste0("chr", chrom))]
genes[, tss := ifelse(strand == "+", start, end)]
# NCBI aliases: map CG ids used in the RNA matrix onto symbols
gi_path <- file.path(root, "ref", "Drosophila_melanogaster.gene_info.gz")
alias <- NULL
if (file.exists(gi_path)) {
  gi <- fread(cmd = paste("zcat", shQuote(gi_path)), header = TRUE)
  # fread may keep the hash
  if (colnames(gi)[1] != "tax_id" && grepl("tax", colnames(gi)[1])) {
    setnames(gi, 1, "tax_id")
  }
  # rebuild map
  rows <- list()
  syms <- gi$Symbol
  loc <- gi$LocusTag
  syns <- gi$Synonyms
  for (i in seq_len(nrow(gi))) {
    toks <- unique(c(syms[i], gsub("^Dmel_", "", loc[i]), strsplit(syns[i], "|", fixed = TRUE)[[1]]))
    toks <- toks[toks != "" & toks != "-"]
    if (length(toks)) rows[[length(rows) + 1]] <- data.frame(alias = toks, symbol = syms[i], stringsAsFactors = FALSE)
  }
  alias <- bind_rows(rows)
  alias <- alias[!duplicated(alias$alias), ]
  fwrite(alias, file.path(res, "gene_alias_map.csv"))
  message("Alias map ", nrow(alias))
}

map_symbol <- function(ids) {
  if (is.null(alias)) return(ids)
  hit <- alias$symbol[match(ids, alias$alias)]
  ifelse(is.na(hit), ids, hit)
}

parse_peak <- function(peak) {
  # chr2L-5464-5963
  m <- regexec("^(chr[0-9A-Za-z]+)-([0-9]+)-([0-9]+)$", peak)
  r <- regmatches(peak, m)
  chrom <- vapply(r, function(z) if (length(z) == 4) z[2] else NA_character_, "")
  start <- vapply(r, function(z) if (length(z) == 4) as.integer(z[3]) else NA_integer_, 1L)
  end <- vapply(r, function(z) if (length(z) == 4) as.integer(z[4]) else NA_integer_, 1L)
  data.frame(peak = peak, chrom = chrom, start = start, end = end, stringsAsFactors = FALSE)
}

nearest_genes <- function(peaks_df, window = 50000) {
  # Distance is from the peak interval to the TSS, matching the LinkPeaks window.
  out <- vector("list", nrow(peaks_df))
  for (i in seq_len(nrow(peaks_df))) {
    ch <- peaks_df$chrom[i]
    if (is.na(ch)) next
    g <- genes[chrom == ch]
    if (!nrow(g)) next
    tss <- g$tss
    dist <- ifelse(tss < peaks_df$start[i], peaks_df$start[i] - tss,
            ifelse(tss > peaks_df$end[i], tss - peaks_df$end[i], 0L))
    near <- which(dist <= window)
    j <- if (length(near)) near[order(dist[near])] else which.min(dist)
    sym <- map_symbol(g$gene_name)
    out[[i]] <- data.frame(
      peak = peaks_df$peak[i],
      nearest_gene = sym[j[1]],
      nearest_dist = as.integer(dist[j[1]]),
      genes_within_50kb = if (length(near)) paste(sym[j], collapse = ";") else NA_character_,
      stringsAsFactors = FALSE
    )
  }
  bind_rows(out)
}

if (!is.null(atac_all) && nrow(atac_all)) {
  sig <- atac_all[atac_all$adj.P.Val < 0.05 & abs(atac_all$logFC) >= 0.5, ]
  message("Significant DARs ", nrow(sig))
  if (nrow(sig)) {
    pk <- parse_peak(sig$feature)
    message("Annotating ", nrow(pk), " peaks")
    ann <- nearest_genes(pk[!duplicated(pk$peak), , drop = FALSE], 50000)
    sig$chrom <- pk$chrom
    sig$start <- pk$start
    sig$end <- pk$end
    sig$nearest_gene <- ann$nearest_gene[match(sig$feature, ann$peak)]
    sig$nearest_dist <- ann$nearest_dist[match(sig$feature, ann$peak)]
    sig$genes_within_50kb <- ann$genes_within_50kb[match(sig$feature, ann$peak)]
    sig$nearest_symbol <- sig$nearest_gene
    fwrite(sig, file.path(res, "dar_significant.csv"))
  } else {
    sig <- sig
    fwrite(sig, file.path(res, "dar_significant.csv"))
  }
} else {
  sig <- NULL
}

# published links
link_sheets <- c(P24 = "P24_links", P48 = "P48_links", Adult = "Adult_links")
link_all <- list()
for (st in names(link_sheets)) {
  d <- as.data.frame(read_excel(file.path(root, "data", "GSE305940_Links.xlsx"), sheet = link_sheets[[st]]))
  d$stage <- st
  d$pvalue <- as.numeric(d$pvalue)
  d$zscore <- as.numeric(d$zscore)
  d$score <- as.numeric(d$score)
  d$zscore.abs <- as.numeric(d$zscore.abs)
  link_all[[st]] <- d
}
links <- bind_rows(link_all)
links$gene_symbol <- map_symbol(links$gene)
fwrite(links[links$gene %in% c("beat-IIIc", "CG14521", "DIP-gamma", "CG15138") |
             links$gene_symbol %in% c("beat-IIIc", "DIP-gamma"), ],
       file.path(res, "published_links_marker_genes.csv"))

cand_rows <- list()
if (!is.null(sig) && nrow(sig) && !is.null(rna_all)) {
  de_sig <- rna_all[rna_all$adj.P.Val < 0.05 & abs(rna_all$logFC) >= 0.25, c("feature", "stage", "logFC", "adj.P.Val", "pct_pDm8", "pct_yDm8")]
  de_sig$symbol <- map_symbol(de_sig$feature)
  for (st in unique(sig$stage)) {
    dar <- sig[sig$stage == st, ]
    de <- de_sig[de_sig$stage == st, ]
    lk <- links[links$stage == st & links$pvalue < 0.05, ]
    if (!nrow(dar) || !nrow(de)) next
    de_by_symbol <- split(de, de$symbol)
    for (i in seq_len(nrow(dar))) {
      near <- unlist(strsplit(ifelse(is.na(dar$genes_within_50kb[i]), "", dar$genes_within_50kb[i]), ";", fixed = TRUE))
      near <- near[near != ""]
      linked <- lk$gene_symbol[lk$peak == dar$feature[i]]
      linked <- linked[!is.na(linked)]
      both <- unique(c(near, linked))
      if (!length(both)) next
      for (g in both) {
        if (!g %in% names(de_by_symbol) && !g %in% de$symbol) next
        hits <- de[de$symbol == g, ]
        if (!nrow(hits)) next
        # concordant: same sign of logFC
        conc <- hits[sign(hits$logFC) == sign(dar$logFC[i]) & hits$logFC != 0, ]
        if (!nrow(conc)) next
        how <- c()
        if (g %in% near) how <- c(how, "TSS_50kb")
        if (g %in% linked) how <- c(how, "published_LinkPeaks")
        cand_rows[[length(cand_rows) + 1]] <- data.frame(
          stage = st,
          peak = dar$feature[i],
          chrom = dar$chrom[i],
          start = dar$start[i],
          end = dar$end[i],
          peak_log2FC_p_over_y = dar$logFC[i],
          peak_FDR = dar$adj.P.Val[i],
          peak_pct_pDm8 = dar$pct_pDm8[i],
          peak_pct_yDm8 = dar$pct_yDm8[i],
          gene = g,
          gene_rna_id = paste(unique(conc$feature), collapse = ";"),
          gene_log2FC_p_over_y = conc$logFC[1],
          gene_FDR = conc$adj.P.Val[1],
          link = paste(how, collapse = "+"),
          higher_in = ifelse(dar$logFC[i] > 0, "pDm8", "yDm8"),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
if (length(cand_rows)) {
  cand <- bind_rows(cand_rows)
  cand <- cand[order(cand$stage, cand$peak_FDR), ]
  fwrite(cand, file.path(res, "candidate_enhancers.csv"))
  message("Candidate enhancer links ", nrow(cand))
} else {
  cand <- NULL
  message("No concordant candidate enhancers")
}

# locus-focused table: all tested ATAC peaks within 50kb of the two markers, whether or not significant
locus_genes <- genes[gene_name %in% c("beat-IIIc", "DIP-gamma")]
if (!is.null(atac_all)) {
  pk_all <- parse_peak(unique(atac_all$feature))
  loc_rows <- list()
  for (i in seq_len(nrow(locus_genes))) {
    g <- locus_genes[i]
    win_a <- g$start - 50000
    win_b <- g$end + 50000
    hit <- pk_all$chrom == g$chrom & pk_all$end >= win_a & pk_all$start <= win_b
    if (!any(hit)) next
    subpk <- pk_all[hit, ]
    sub <- atac_all[atac_all$feature %in% subpk$peak, , drop = FALSE]
    sub <- cbind(sub, subpk[match(sub$feature, subpk$peak), c("chrom", "start", "end")])
    sub$locus_gene <- g$gene_name
    loc_rows[[length(loc_rows) + 1]] <- sub
  }
  if (length(loc_rows)) {
    locus <- bind_rows(loc_rows)
    fwrite(locus, file.path(res, "locus_peaks_markers.csv"))
    # plot significant or top peaks as bars of logFC
    locus$mid <- (locus$start + locus$end) / 2
    locus$sig <- locus$adj.P.Val < 0.05
    p_loc <- ggplot(locus, aes(mid / 1e6, logFC, color = sig)) +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      geom_point(size = 1.4) +
      facet_grid(locus_gene ~ stage, scales = "free_x") +
      scale_color_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#B0B0B0"),
                         name = "FDR<0.05") +
      labs(
        title = "Differential accessibility within 50 kb of beat-IIIc and DIP-gamma",
        subtitle = "Each point is one ATAC peak. Positive log2 fold-change is higher in pDm8.",
        x = "Genomic position (Mb)",
        y = "log2 fold-change (pDm8 / yDm8)"
      ) +
      theme_bw(base_size = 10)
    ggsave(file.path(fig, "08_locus_accessibility.png"), p_loc, width = 10, height = 6.2, dpi = 160)
  }
}

# TF-like genes among RNA DE results, using NCBI description
if (!is.null(rna_all) && file.exists(gi_path)) {
  gi <- fread(cmd = paste("zcat", shQuote(gi_path)))
  if (!"Symbol" %in% names(gi)) setnames(gi, 2, sub("^#", "", names(gi)[2]))
  # columns may include #tax_id
  desc_col <- intersect(c("description", "Full_name_from_nomenclature_authority"), names(gi))[1]
  gi$is_tf <- grepl("transcription factor|homeobox|nuclear receptor|forkhead|bHLH|helix-loop-helix|T-box|zinc finger protein",
                    gi[[desc_col]], ignore.case = TRUE)
  tf_symbols <- unique(gi$Symbol[gi$is_tf])
  # also CG aliases of those TFs
  tf_alias <- alias$alias[alias$symbol %in% tf_symbols]
  rna_all$symbol <- map_symbol(rna_all$feature)
  rna_all$is_tf <- rna_all$symbol %in% tf_symbols | rna_all$feature %in% tf_alias
  tf_de <- rna_all[rna_all$is_tf & rna_all$adj.P.Val < 0.1, ]
  tf_de <- tf_de[order(tf_de$stage, tf_de$adj.P.Val), ]
  fwrite(tf_de, file.path(res, "tf_differential_expression.csv"))
  message("TF-like genes at FDR<0.1 ", nrow(tf_de))
  fwrite(rna_all, file.path(res, "de_rna_pseudobulk.csv"))
}

# compact summary json-like tsv
summ <- data.frame(
  item = c("n_Dm8", "n_P0", "n_P24", "n_P48", "n_Adult"),
  value = c(nrow(labels),
            sum(labels$stage == "P0"),
            sum(labels$stage == "P24"),
            sum(labels$stage == "P48"),
            sum(labels$stage == "Adult"))
)
fwrite(summ, file.path(res, "summary_counts.csv"))
message("DONE ", Sys.time())
saveRDS(list(rna = rna_res$tables, atac = atac_res$tables, exclusivity = excl),
        file.path(res, "analysis_objects.rds"))

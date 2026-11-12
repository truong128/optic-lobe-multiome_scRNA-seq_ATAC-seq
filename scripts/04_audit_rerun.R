#!/usr/bin/env Rscript
# Audit revision of the pDm8 vs yDm8 tests.
# Writes the individual-line tables into results/.
# Changes requested by the audit:
#   1. LinkPeaks uses gene_symbol, not the missing column symbol.
#   2. A gene is near a peak when its TSS, not its whole body, is within 50 kb.
#   3. The paired unit is the individual DGRP line (dgrp), not the pooled dgrp.new group.
#   4. Sensitivity: voom, edgeR quasi-likelihood, female-only lines, leave-one-line-out.
user_lib <- path.expand("~/R/library")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
  library(limma)
  library(edgeR)
  library(ggplot2)
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
out <- file.path(root, "results")
fig <- file.path(root, "figures")
dir.create(out, showWarnings = FALSE, recursive = TRUE)
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

log_con <- file(file.path(out, "04_audit_log.txt"), open = "wt")
sink(log_con, split = TRUE)
on.exit({
  sink()
  close(log_con)
}, add = TRUE)
message("START ", Sys.time())

seu <- readRDS(file.path(root, "results", "dm8_subset.rds"))
md <- seu@meta.data
stopifnot(all(md$AnnotatedJun26 == "Dm8"))
rna_counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
atac_counts <- GetAssayData(seu, assay = "ATAC", layer = "counts")
message("RNA ", nrow(rna_counts), " x ", ncol(rna_counts),
        " ATAC ", nrow(atac_counts), " x ", ncol(atac_counts))

beat <- as.numeric(rna_counts["beat-IIIc", ])
dip <- as.numeric(rna_counts["CG14521", ])
subtype <- ifelse(beat > 0 & dip == 0, "pDm8",
           ifelse(dip > 0 & beat == 0, "yDm8",
           ifelse(beat > 0 & dip > 0, "double_positive", "double_negative")))

meta <- data.frame(
  cell = colnames(seu),
  stage = as.character(md$stage),
  sex = as.character(md$sex),
  dgrp = as.character(md$dgrp),
  dgrp_new = as.character(md$dgrp.new),
  orig.ident = as.character(md$orig.ident),
  subtype = subtype,
  stringsAsFactors = FALSE
)

tested <- meta[meta$subtype %in% c("pDm8", "yDm8") & meta$stage %in% c("P24", "P48"), ]
sex_line <- tested %>%
  group_by(stage, dgrp, sex) %>%
  summarise(n = n(), .groups = "drop")
fwrite(sex_line, file.path(out, "sex_by_line.csv"))
lib_line <- tested %>%
  group_by(stage, dgrp, orig.ident) %>%
  summarise(n = n(), .groups = "drop")
fwrite(lib_line, file.path(out, "library_by_line.csv"))
nsex <- sex_line %>% group_by(stage, dgrp) %>% summarise(n_sex = n_distinct(sex), .groups = "drop")
message("Lines with more than one sex label: ", sum(nsex$n_sex > 1), " of ", nrow(nsex))

prepare_pb <- function(counts, meta, stage, genotype_col,
                       min_cells = 5, min_genotypes = 4, min_pct = 0.10,
                       sex_keep = NULL, drop_genos = character()) {
  md <- meta[meta$stage == stage & meta$subtype %in% c("pDm8", "yDm8"), ]
  if (!is.null(sex_keep)) md <- md[md$sex %in% sex_keep, ]
  if (length(drop_genos)) md <- md[!md[[genotype_col]] %in% drop_genos, ]
  md <- md[!is.na(md[[genotype_col]]) & nzchar(md[[genotype_col]]) & md[[genotype_col]] != "0", ]
  md <- md[md$cell %in% colnames(counts), ]
  if (!nrow(md)) return(NULL)
  tab <- table(md[[genotype_col]], md$subtype)
  if (!all(c("pDm8", "yDm8") %in% colnames(tab))) return(NULL)
  good <- rownames(tab)[tab[, "pDm8"] >= min_cells & tab[, "yDm8"] >= min_cells]
  if (length(good) < min_genotypes) return(NULL)
  md <- md[md[[genotype_col]] %in% good, ]
  samples <- expand.grid(genotype = good, subtype = c("yDm8", "pDm8"),
                         KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  samples$sample_id <- paste(samples$genotype, samples$subtype, sep = "__")
  pb_cols <- vector("list", nrow(samples))
  samples$n_cells <- NA_integer_
  for (i in seq_len(nrow(samples))) {
    cells <- md$cell[md[[genotype_col]] == samples$genotype[i] & md$subtype == samples$subtype[i]]
    pb_cols[[i]] <- Matrix::rowSums(counts[, cells, drop = FALSE])
    samples$n_cells[i] <- length(cells)
  }
  pb <- do.call(cbind, pb_cols)
  colnames(pb) <- samples$sample_id
  cells_p <- md$cell[md$subtype == "pDm8"]
  cells_y <- md$cell[md$subtype == "yDm8"]
  pct_p <- Matrix::rowMeans(counts[, cells_p, drop = FALSE] > 0)
  pct_y <- Matrix::rowMeans(counts[, cells_y, drop = FALSE] > 0)
  keep <- which((pct_p >= min_pct) | (pct_y >= min_pct))
  pb <- pb[keep, , drop = FALSE]
  lib <- colSums(pb)
  if (any(lib <= 0)) {
    samples <- samples[lib > 0, ]
    pb <- pb[, lib > 0, drop = FALSE]
    still <- names(which(table(samples$genotype) == 2))
    if (length(still) < min_genotypes) return(NULL)
    samples <- samples[samples$genotype %in% still, ]
    pb <- pb[, samples$sample_id, drop = FALSE]
    good <- still
  }
  list(
    pb = pb,
    samples = samples,
    pct_p = pct_p[rownames(pb)],
    pct_y = pct_y[rownames(pb)],
    n_p = length(cells_p),
    n_y = length(cells_y),
    n_genotypes = length(good),
    good = good
  )
}

test_prep <- function(prep, method = c("trend", "voom", "edger")) {
  method <- match.arg(method)
  samples <- prep$samples
  y <- DGEList(counts = prep$pb)
  y <- calcNormFactors(y, method = "TMM")
  subtype_f <- factor(samples$subtype, levels = c("yDm8", "pDm8"))
  genotype_f <- factor(samples$genotype)
  design <- model.matrix(~ genotype_f + subtype_f)
  coef <- grep("^subtype_f", colnames(design), value = TRUE)
  if (length(coef) != 1) stop("Expected one subtype coefficient, got: ", paste(coef, collapse = ", "))
  if (method == "trend") {
    logcpm <- cpm(y, log = TRUE, prior.count = 1)
    fit <- eBayes(lmFit(logcpm, design), trend = TRUE, robust = TRUE)
    tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  } else if (method == "voom") {
    fit <- eBayes(lmFit(voom(y, design), design), robust = TRUE)
    tt <- topTable(fit, coef = coef, number = Inf, sort.by = "none")
  } else {
    y <- estimateDisp(y, design, robust = TRUE)
    qt <- glmQLFTest(glmQLFit(y, design, robust = TRUE), coef = coef)
    raw <- topTags(qt, n = Inf, sort.by = "none")$table
    tt <- data.frame(
      logFC = raw$logFC,
      P.Value = raw$PValue,
      adj.P.Val = raw$FDR,
      feature = rownames(raw),
      row.names = rownames(raw),
      stringsAsFactors = FALSE
    )
  }
  tt$feature <- rownames(tt)
  tt$FDR <- tt$adj.P.Val
  logcpm <- cpm(y, log = TRUE, prior.count = 1)
  p_idx <- which(samples$subtype == "pDm8")
  y_idx <- which(samples$subtype == "yDm8")
  ord <- match(samples$genotype[p_idx], samples$genotype[y_idx])
  delta <- logcpm[, p_idx, drop = FALSE] - logcpm[, y_idx[ord], drop = FALSE]
  delta <- delta[tt$feature, , drop = FALSE]
  mean_delta <- rowMeans(delta)
  tt$frac_genotypes_same_sign <- rowMeans(sign(delta) == sign(mean_delta))
  tt$mean_paired_log2FC <- as.numeric(mean_delta)
  tt$pct_pDm8 <- prep$pct_p[tt$feature]
  tt$pct_yDm8 <- prep$pct_y[tt$feature]
  tt$n_pDm8 <- prep$n_p
  tt$n_yDm8 <- prep$n_y
  tt$n_genotypes <- prep$n_genotypes
  tt$method <- method
  tt[, c("feature", "logFC", "P.Value", "FDR", "pct_pDm8", "pct_yDm8",
         "n_pDm8", "n_yDm8", "n_genotypes", "frac_genotypes_same_sign",
         "mean_paired_log2FC", "method")]
}

run_models <- function(counts, assay, min_pct, genotype_col, methods, sex_keep = NULL, tag = "line") {
  tabs <- list()
  for (st in c("P24", "P48")) {
    prep <- prepare_pb(counts, meta, st, genotype_col, min_pct = min_pct, sex_keep = sex_keep)
    if (is.null(prep)) {
      message("SKIP ", assay, " ", st, " ", tag)
      next
    }
    message(assay, " ", st, " ", tag, " | lines ", prep$n_genotypes,
            " | cells p/y ", prep$n_p, "/", prep$n_y,
            " | features ", nrow(prep$pb))
    for (m in methods) {
      tt <- test_prep(prep, m)
      tt$stage <- st
      tt$assay <- assay
      tt$model <- tag
      tabs[[paste(st, m, sep = "__")]] <- tt
      message("  ", m, " sig ", sum(tt$FDR < 0.05 & abs(tt$logFC) >= 0.5, na.rm = TRUE))
    }
  }
  if (!length(tabs)) return(NULL)
  rbindlist(tabs)
}

message("\n===== individual DGRP line models =====")
rna_line <- run_models(rna_counts, "RNA", 0.10, "dgrp", c("trend", "voom", "edger"), tag = "dgrp_line")
atac_line <- run_models(atac_counts, "ATAC", 0.05, "dgrp", c("trend", "voom", "edger"), tag = "dgrp_line")
message("\n===== female-only, individual lines, limma-trend =====")
rna_fem <- run_models(rna_counts, "RNA", 0.10, "dgrp", "trend", sex_keep = "female", tag = "female_dgrp_line")
atac_fem <- run_models(atac_counts, "ATAC", 0.05, "dgrp", "trend", sex_keep = "female", tag = "female_dgrp_line")

fwrite(rna_line, file.path(out, "de_rna_line_models.csv"))
fwrite(atac_line, file.path(out, "de_atac_line_models.csv"))
fwrite(rna_fem, file.path(out, "de_rna_female_trend.csv"))
fwrite(atac_fem, file.path(out, "de_atac_female_trend.csv"))
message("Wrote full model tables ", Sys.time())

base_rna <- fread(file.path(root, "results", "de_rna_pseudobulk.csv"))
base_atac <- fread(file.path(root, "results", "de_atac_pseudobulk.csv"))
base_rna[, assay := "RNA"]
base_atac[, assay := "ATAC"]

parse_peak <- function(peak) {
  m <- regexec("^(chr[0-9A-Za-z]+)-([0-9]+)-([0-9]+)$", peak)
  r <- regmatches(peak, m)
  data.frame(
    peak = peak,
    chrom = vapply(r, function(z) if (length(z) == 4) z[2] else NA_character_, ""),
    start = vapply(r, function(z) if (length(z) == 4) as.integer(z[3]) else NA_integer_, 1L),
    end = vapply(r, function(z) if (length(z) == 4) as.integer(z[4]) else NA_integer_, 1L),
    stringsAsFactors = FALSE
  )
}

vsx_chr <- "chrX"
vsx_start <- 5528927L
vsx_end <- 5560183L
pk <- parse_peak(unique(base_atac$feature))
in_vsx <- pk$chrom == vsx_chr & pk$end >= vsx_start & pk$start <= vsx_end
vsx_peaks <- pk$peak[in_vsx]
base_vsx <- base_atac[feature %in% vsx_peaks & adj.P.Val < 0.05 & abs(logFC) >= 0.5]
vsx_both <- base_vsx[, .N, by = feature][N == 2, feature]
message("Vsx2 peaks significant at both stages in the baseline: ", length(vsx_both))

focus_rna <- c("beat-IIIc", "CG14521", "nAChRalpha6", "CG44153", "Tl", "Sulf1", "Drl-2", "kek6")
focus_atac <- vsx_both

loo_one <- function(counts, assay, min_pct, focus) {
  rows <- list()
  for (st in c("P24", "P48")) {
    full <- prepare_pb(counts, meta, st, "dgrp", min_pct = min_pct)
    if (is.null(full)) next
    for (g in full$good) {
      prep <- prepare_pb(counts, meta, st, "dgrp", min_pct = min_pct, drop_genos = g)
      if (is.null(prep)) next
      tt <- test_prep(prep, "trend")
      got <- match(focus, tt$feature)
      rows[[paste(st, g)]] <- data.frame(
        assay = assay,
        stage = st,
        dropped_line = g,
        feature = focus,
        logFC = tt$logFC[got],
        FDR = tt$FDR[got],
        stringsAsFactors = FALSE
      )
      message("  LOO ", assay, " ", st, " drop ", g, " features kept ", sum(!is.na(got)))
    }
  }
  if (!length(rows)) return(NULL)
  rbindlist(rows)
}

message("\n===== leave-one-line-out =====")
loo_rna <- loo_one(rna_counts, "RNA", 0.10, focus_rna)
loo_atac <- loo_one(atac_counts, "ATAC", 0.05, focus_atac)
loo <- rbindlist(list(loo_rna, loo_atac), fill = TRUE)
fwrite(loo, file.path(out, "leave_one_line_out.csv"))

pick <- function(tab, feat, st, method = NULL) {
  # Use a data.frame. data.table would treat the argument names as columns.
  if (is.null(tab) || !nrow(tab)) return(list(logFC = NA_real_, FDR = NA_real_, frac = NA_real_))
  hit <- as.data.frame(tab)
  hit <- hit[hit$feature == feat & hit$stage == st, , drop = FALSE]
  if (!is.null(method)) hit <- hit[hit$method == method, , drop = FALSE]
  if (!nrow(hit)) return(list(logFC = NA_real_, FDR = NA_real_, frac = NA_real_))
  list(logFC = hit$logFC[1], FDR = hit$FDR[1],
       frac = if ("frac_genotypes_same_sign" %in% names(hit)) hit$frac_genotypes_same_sign[1] else NA_real_)
}

cmp_rows <- list()
for (assay in c("RNA", "ATAC")) {
  feats <- if (assay == "RNA") focus_rna else focus_atac
  base <- if (assay == "RNA") base_rna else base_atac
  line <- if (assay == "RNA") rna_line else atac_line
  fem <- if (assay == "RNA") rna_fem else atac_fem
  for (st in c("P24", "P48")) {
    for (f in feats) {
      b <- as.data.frame(base)
      b <- b[b$feature == f & b$stage == st, , drop = FALSE]
      tr <- pick(line, f, st, "trend")
      vo <- pick(line, f, st, "voom")
      ed <- pick(line, f, st, "edger")
      fe <- pick(fem, f, st, "trend")
      cmp_rows[[paste(assay, st, f)]] <- data.frame(
        assay = assay, stage = st, feature = f,
        base_logFC = if (nrow(b)) b$logFC[1] else NA_real_,
        base_FDR = if (nrow(b)) b$adj.P.Val[1] else NA_real_,
        base_frac = if (nrow(b)) b$frac_genotypes_same_sign[1] else NA_real_,
        line_logFC = tr$logFC, line_FDR = tr$FDR, line_frac = tr$frac,
        voom_logFC = vo$logFC, voom_FDR = vo$FDR,
        edger_logFC = ed$logFC, edger_FDR = ed$FDR,
        female_logFC = fe$logFC, female_FDR = fe$FDR,
        stringsAsFactors = FALSE
      )
    }
  }
}
cmp <- rbindlist(cmp_rows)
loo_ref <- cmp[, .(assay, stage, feature, line_logFC)]
loo2 <- merge(loo, loo_ref, by = c("assay", "stage", "feature"), all.x = TRUE)
loo_sum <- loo2[!is.na(logFC), .(
  loo_n = .N,
  loo_same_sign = sum(sign(logFC) == sign(line_logFC) & line_logFC != 0),
  loo_fdr05 = sum(FDR < 0.05),
  loo_min_abs_logFC = min(abs(logFC)),
  loo_max_FDR = max(FDR)
), by = .(assay, stage, feature)]
cmp <- merge(cmp, loo_sum, by = c("assay", "stage", "feature"), all.x = TRUE)
cmp[, survives_models := line_FDR < 0.05 & voom_FDR < 0.05 & edger_FDR < 0.05 &
      female_FDR < 0.05 & sign(line_logFC) == sign(base_logFC) &
      sign(voom_logFC) == sign(base_logFC) & sign(edger_logFC) == sign(base_logFC) &
      sign(female_logFC) == sign(base_logFC)]
cmp[, survives_all_LOO := !is.na(loo_n) & loo_n > 0 & loo_fdr05 == loo_n & loo_same_sign == loo_n]
cmp[, survives := survives_models & survives_all_LOO]
setorder(cmp, assay, stage, feature)
fwrite(cmp, file.path(out, "focus_comparison.csv"))
message("Focus features that keep FDR<0.05 and the same sign in every model: ",
        sum(cmp$survives, na.rm = TRUE), " of ", nrow(cmp))

count_sig <- function(tab, assay, model) {
  if (is.null(tab) || !nrow(as.data.frame(tab))) return(NULL)
  x <- as.data.frame(tab)
  x <- x[x$FDR < 0.05 & abs(x$logFC) >= 0.5, ]
  if (!nrow(x)) return(data.frame(assay = assay, model = model, stage = c("P24", "P48"), n_sig = 0))
  agg <- aggregate(feature ~ stage, x, function(z) length(unique(z)))
  names(agg)[2] <- "n_sig"
  agg$assay <- assay
  agg$model <- model
  agg
}
sig_summary <- rbind(
  data.frame(assay = "RNA", model = "baseline_dgrp.new", stage = c("P24", "P48"),
             n_sig = c(sum(base_rna$stage == "P24" & base_rna$adj.P.Val < 0.05 & abs(base_rna$logFC) >= 0.5),
                       sum(base_rna$stage == "P48" & base_rna$adj.P.Val < 0.05 & abs(base_rna$logFC) >= 0.5))),
  data.frame(assay = "ATAC", model = "baseline_dgrp.new", stage = c("P24", "P48"),
             n_sig = c(sum(base_atac$stage == "P24" & base_atac$adj.P.Val < 0.05 & abs(base_atac$logFC) >= 0.5),
                       sum(base_atac$stage == "P48" & base_atac$adj.P.Val < 0.05 & abs(base_atac$logFC) >= 0.5))),
  count_sig(rna_line[method == "trend"], "RNA", "line_trend"),
  count_sig(rna_line[method == "voom"], "RNA", "line_voom"),
  count_sig(rna_line[method == "edger"], "RNA", "line_edger"),
  count_sig(rna_fem, "RNA", "female_trend"),
  count_sig(atac_line[method == "trend"], "ATAC", "line_trend"),
  count_sig(atac_line[method == "voom"], "ATAC", "line_voom"),
  count_sig(atac_line[method == "edger"], "ATAC", "line_edger"),
  count_sig(atac_fem, "ATAC", "female_trend")
)
fwrite(sig_summary, file.path(out, "sig_counts.csv"))
print(sig_summary)

# ---------- revised candidate enhancers: TSS +/- 50 kb and LinkPeaks gene_symbol ----------
alias <- fread(file.path(root, "results", "gene_alias_map.csv"))
map_symbol <- function(ids) {
  hit <- alias$symbol[match(ids, alias$alias)]
  ifelse(is.na(hit) | hit == "", ids, hit)
}
genes <- fread(file.path(root, "ref", "genes.tsv"))
genes[, start := as.integer(start)]
genes[, end := as.integer(end)]
genes[, chrom := ifelse(startsWith(as.character(chrom), "chr"), as.character(chrom), paste0("chr", chrom))]
genes[, tss := ifelse(strand == "+", start, end)]
genes[, symbol := map_symbol(gene_name)]

link_sheets <- c(P24 = "P24_links", P48 = "P48_links")
link_all <- rbindlist(lapply(names(link_sheets), function(st) {
  d <- as.data.table(read_excel(file.path(root, "data", "GSE305940_Links.xlsx"), sheet = link_sheets[[st]]))
  d[, stage := st]
  d[, pvalue := as.numeric(pvalue)]
  d[, gene_symbol := map_symbol(gene)]
  d[pvalue < 0.05, .(stage, peak, gene, gene_symbol, pvalue)]
}))

rna_trend <- rna_line[method == "trend"]
atac_trend <- atac_line[method == "trend"]
sig <- atac_trend[FDR < 0.05 & abs(logFC) >= 0.5]
de <- rna_trend[FDR < 0.05 & abs(logFC) >= 0.25]
de[, symbol := map_symbol(feature)]

cand <- list()
for (i in seq_len(nrow(sig))) {
  st <- sig$stage[i]
  loc <- parse_peak(sig$feature[i])
  ch <- loc$chrom
  g <- genes[chrom == ch & tss >= (loc$start - 50000L) & tss <= (loc$end + 50000L)]
  near <- unique(g$symbol)
  linked <- unique(link_all[stage == st & peak == sig$feature[i], gene_symbol])
  both <- unique(c(near, linked))
  both <- both[!is.na(both) & nzchar(both)]
  if (!length(both)) next
  de_st <- de[stage == st]
  for (gs in both) {
    hits <- de_st[symbol == gs]
    if (!nrow(hits)) next
    hits <- hits[sign(logFC) == sign(sig$logFC[i]) & logFC != 0]
    if (!nrow(hits)) next
    how <- c()
    if (gs %in% near) how <- c(how, "TSS_50kb")
    if (gs %in% linked) how <- c(how, "published_LinkPeaks")
    cand[[length(cand) + 1]] <- data.frame(
      stage = st,
      peak = sig$feature[i],
      chrom = loc$chrom,
      start = loc$start,
      end = loc$end,
      peak_log2FC_p_over_y = sig$logFC[i],
      peak_FDR = sig$FDR[i],
      peak_pct_pDm8 = sig$pct_pDm8[i],
      peak_pct_yDm8 = sig$pct_yDm8[i],
      gene = gs,
      gene_rna_id = paste(unique(hits$feature), collapse = ";"),
      gene_log2FC_p_over_y = hits$logFC[1],
      gene_FDR = hits$FDR[1],
      link = paste(how, collapse = "+"),
      higher_in = ifelse(sig$logFC[i] > 0, "pDm8", "yDm8"),
      stringsAsFactors = FALSE
    )
  }
}
cand <- if (length(cand)) rbindlist(cand) else data.table()
if (nrow(cand)) setorder(cand, stage, gene, peak_FDR)
fwrite(cand, file.path(out, "candidate_enhancers_revised.csv"))
message("Revised candidate links: ", nrow(cand),
        " | with published_LinkPeaks: ", if (nrow(cand)) sum(grepl("published_LinkPeaks", cand$link)) else 0)

# comparison plot of the focus estimates
long <- rbindlist(list(
  cmp[, .(assay, stage, feature, model = "baseline", logFC = base_logFC, FDR = base_FDR)],
  cmp[, .(assay, stage, feature, model = "line trend", logFC = line_logFC, FDR = line_FDR)],
  cmp[, .(assay, stage, feature, model = "voom", logFC = voom_logFC, FDR = voom_FDR)],
  cmp[, .(assay, stage, feature, model = "edgeR QL", logFC = edger_logFC, FDR = edger_FDR)],
  cmp[, .(assay, stage, feature, model = "female only", logFC = female_logFC, FDR = female_FDR)]
))
long[, model := factor(model, levels = c("baseline", "line trend", "voom", "edgeR QL", "female only"))]
long[, sig := FDR < 0.05]
long[, label := ifelse(assay == "ATAC",
                       paste0("Vsx2 ", sub("^chrX-", "", feature)),
                       feature)]
p <- ggplot(long, aes(logFC, label, color = sig, shape = model)) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_point(position = position_dodge(width = 0.6), size = 2) +
  scale_color_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#B0B0B0"),
                     labels = c(`TRUE` = "FDR < 0.05", `FALSE` = "FDR >= 0.05")) +
  facet_grid(assay ~ stage, scales = "free_y") +
  labs(x = "log2 fold-change, pDm8 / yDm8", y = NULL,
       title = "Focus genes and Vsx2 peaks across models",
       color = NULL, shape = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")
ggsave(file.path(fig, "10_sensitivity.png"), p, width = 11, height = 7, dpi = 150)
message("DONE ", Sys.time())

#!/usr/bin/env Rscript
# Replace the pooled-model figures with the individual-line limma-trend results.
user_lib <- path.expand("~/R/library")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
root <- Sys.getenv("PROJECT_ROOT", unset = "")
if (!nzchar(root)) {
  root <- if (length(file_arg)) normalizePath(file.path(dirname(file_arg[1]), "..")) else getwd()
}
fig <- file.path(root, "figures")
rna <- fread(file.path(root, "results", "de_rna_line_models.csv"))
atac <- fread(file.path(root, "results", "de_atac_line_models.csv"))
rna <- rna[method == "trend"]
atac <- atac[method == "trend"]

volcano_plot <- function(tt, title, label_n, outfile) {
  tt <- as.data.frame(tt)
  tt$neglog10 <- -log10(pmax(tt$P.Value, 1e-300))
  tt$sig <- tt$FDR < 0.05 & abs(tt$logFC) >= 0.5
  controls <- tt$feature %in% c("beat-IIIc", "CG14521", "DIP-gamma")
  top <- tt[order(tt$FDR, -abs(tt$logFC)), ]
  top <- top[!(top$feature %in% c("beat-IIIc", "CG14521", "DIP-gamma")), ]
  lab_feat <- unique(c(tt$feature[controls], head(top$feature, label_n)))
  tt$label <- ifelse(tt$feature %in% lab_feat, tt$feature, "")
  p <- ggplot(tt, aes(logFC, neglog10, color = sig)) +
    geom_point(size = 0.7, alpha = 0.7) +
    geom_text_repel(aes(label = label), size = 2.6, max.overlaps = 30,
                    color = "black", show.legend = FALSE) +
    scale_color_manual(
      values = c(`TRUE` = "#0072B2", `FALSE` = "#B0B0B0"),
      labels = c(`TRUE` = "FDR<0.05 and |log2FC|>=0.5", `FALSE` = "not significant"),
      name = NULL
    ) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, linewidth = 0.3) +
    labs(
      title = title,
      subtitle = "Individual DGRP lines, limma-trend. Positive log2 fold-change is higher in pDm8.",
      x = "log2 fold-change (pDm8 / yDm8)",
      y = "-log10(p-value)"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")
  ggsave(outfile, p, width = 7.4, height = 5.4, dpi = 160)
}

for (st in c("P24", "P48")) {
  volcano_plot(
    rna[stage == st],
    paste0("RNA: pDm8 vs yDm8 at ", st),
    12,
    file.path(fig, paste0("05_volcano_rna_", st, ".png"))
  )
  volcano_plot(
    atac[stage == st],
    paste0("ATAC: pDm8 vs yDm8 at ", st),
    8,
    file.path(fig, paste0("06_volcano_atac_", st, ".png"))
  )
}

sig <- rna[FDR < 0.05 & abs(logFC) >= 0.5]
wide <- dcast(sig, feature ~ stage, value.var = "logFC")
for (st in c("P24", "P48")) if (!st %in% names(wide)) wide[[st]] <- NA_real_
wide <- as.data.frame(wide)
wide$mean_fc <- rowMeans(wide[, c("P24", "P48")], na.rm = TRUE)
wide <- wide[order(wide$mean_fc), ]
long <- rbind(
  data.frame(gene = wide$feature, stage = "P24", logFC = wide$P24),
  data.frame(gene = wide$feature, stage = "P48", logFC = wide$P48)
)
long$gene <- factor(long$gene, levels = wide$feature)
p_hm <- ggplot(long, aes(stage, gene, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#E69F00", mid = "white", high = "#0072B2",
    midpoint = 0, na.value = "grey90", name = "log2FC"
  ) +
  labs(
    title = "Genes differential between pDm8 and yDm8",
    subtitle = "Individual-line limma-trend. FDR < 0.05 and |log2FC| >= 0.5 in at least one stage. Blue is higher in pDm8.",
    x = NULL, y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(axis.text.y = element_text(size = 7))
ggsave(
  file.path(fig, "07_heatmap_top_genes.png"),
  p_hm, width = 6.2, height = max(5, 0.18 * nrow(wide) + 1.5), dpi = 160
)

genes <- fread(file.path(root, "ref", "genes.tsv"))
want <- genes[gene_name %in% c("beat-IIIc", "DIP-gamma", "Vsx2")]
want[, chrom := ifelse(startsWith(chrom, "chr"), chrom, paste0("chr", chrom))]
want[, tss := ifelse(strand == "+", start, end)]
peaks <- atac[, .(
  stage, logFC, FDR,
  chrom = sub("-.*", "", feature),
  start = as.integer(sub(".*-(\\d+)-\\d+$", "\\1", feature)),
  end = as.integer(sub(".*-(\\d+)$", "\\1", feature))
)]
loc_rows <- list()
for (i in seq_len(nrow(want))) {
  g <- want[i]
  hit <- peaks[chrom == g$chrom & end >= g$tss - 50000 & start <= g$tss + 50000]
  if (!nrow(hit)) next
  hit[, locus := g$gene_name]
  loc_rows[[length(loc_rows) + 1L]] <- hit
}
loc <- rbindlist(loc_rows)
loc[, mid := (start + end) / 2]
loc[, sig := FDR < 0.05 & abs(logFC) >= 0.5]
loc[, locus := factor(locus, levels = c("Vsx2", "DIP-gamma", "beat-IIIc"))]
p_loc <- ggplot(loc, aes(mid / 1e6, logFC, color = sig)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_point(size = 1.4) +
  facet_grid(locus ~ stage, scales = "free_x") +
  scale_color_manual(
    values = c(`TRUE` = "#0072B2", `FALSE` = "#B0B0B0"),
    name = "FDR<0.05 and |log2FC|>=0.5"
  ) +
  labs(
    title = "Differential accessibility within 50 kb of the TSS",
    subtitle = "Individual-line limma-trend. Positive log2 fold-change is higher in pDm8.",
    x = "Genomic position (Mb)",
    y = "log2 fold-change (pDm8 / yDm8)"
  ) +
  theme_bw(base_size = 10)
ggsave(file.path(fig, "08_locus_accessibility.png"), p_loc, width = 10, height = 7.2, dpi = 160)
message("Replaced figures 05-08 from the line-trend tables")

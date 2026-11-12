#!/usr/bin/env Rscript
# Coverage tracks at three loci from the fragment file, for labeled Dm8 nuclei.
# The fragment file is queried only at these loci. It is not a genome-wide reanalysis.
user_lib <- path.expand("~/R/library")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicRanges)
  library(IRanges)
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
root <- Sys.getenv("PROJECT_ROOT", unset = "")
if (!nzchar(root)) {
  root <- if (length(file_arg)) normalizePath(file.path(dirname(file_arg[1]), "..")) else getwd()
}
frag <- file.path(root, "data", "GSE305940_atac_fragments.tsv.gz")
labels <- fread(file.path(root, "results", "cell_labels.csv"))
map <- fread(file.path(root, "results", "barcode_map.tsv"))
labels <- merge(labels, map[, .(cell, atac)], by = "cell", all.x = TRUE)

use <- labels[subtype %in% c("pDm8", "yDm8") & stage %in% c("P24", "P48")]
use[, cell_suf := sub(".*-", "", cell)]
use[, atac_key := paste0(sub("-.*", "", atac), "-", cell_suf)]
use[, gex_key := cell]

message("Cells in tracks: ", nrow(use))
print(use[, .N, by = .(stage, subtype)])

# Gene coordinates are dm6 / Ensembl, chromosome names without the chr prefix,
# matching the fragment file.
loci <- data.table(
  locus = c("Vsx2", "DIP-gamma", "beat-IIIc"),
  chrom = c("X", "3R", "2L"),
  start = c(5523927L, 28920004L, 17182637L),
  end = c(5565183L, 28975372L, 17265734L),
  gene_start = c(5528927L, 28925004L, 17187637L),
  gene_end = c(5560183L, 28970372L, 17260734L),
  strand = c("-", "-", "-")
)
bin_size <- 100L

tb <- TabixFile(frag)
gr <- GRanges(loci$chrom, IRanges(loci$start, loci$end))
names(gr) <- loci$locus
message("Querying fragments...")
raw <- scanTabix(tb, param = gr)
message("Returned regions: ", paste(names(raw), collapse = " | "))
message("Fragment lines: ", paste(vapply(raw, length, 1L), collapse = ", "))

choose_key <- function(lines, gex_keys, atac_keys) {
  if (!length(lines)) return("gex_key")
  take <- lines[seq_len(min(length(lines), 200000L))]
  bc <- vapply(strsplit(take, "\t", fixed = TRUE), function(z) if (length(z) >= 4) z[4] else NA_character_, "")
  g <- sum(bc %in% gex_keys)
  a <- sum(bc %in% atac_keys)
  message("Barcode match in sample: cell-name ", g, "  atac+library ", a, "  of ", length(bc))
  if (a > g) "atac_key" else "gex_key"
}

key_col <- choose_key(raw[[1]], use$gex_key, use$atac_key)
barcodes <- use[[key_col]]
meta <- use[, .(stage, subtype, nCount_ATAC)]
meta[, barcode := barcodes]
meta <- meta[!duplicated(barcode)]
message("Using barcode column ", key_col, " n barcodes ", nrow(meta))

cov_rows <- list()
gene_rows <- list()
support_rows <- list()
for (i in seq_len(nrow(loci))) {
  loc <- loci[i]
  lines <- raw[[i]]
  message(loc$locus, " fragments in window (all cells) ", length(lines))
  if (!length(lines)) next
  # 10x fragments.tsv: chrom, 0-based start, half-open end, barcode, readSupport.
  # readSupport counts read pairs, including PCR duplicates. Each record is one
  # unique fragment, so coverage uses weight 1. IRanges is 1-based and inclusive,
  # so the 1-based interval is (start + 1) through end.
  parts <- tstrsplit(lines, "\t", fixed = TRUE, keep = 1:5)
  dt <- data.table(
    start = as.integer(parts[[2]]),
    end = as.integer(parts[[3]]),
    barcode = parts[[4]],
    read_support = as.integer(parts[[5]])
  )
  dt <- dt[!is.na(barcode)]
  dt <- meta[dt, on = "barcode", nomatch = 0]
  n_multi <- dt[read_support > 1, .N]
  message(
    "  kept Dm8 p/y unique fragments ", nrow(dt),
    "; readSupport>1 ", n_multi,
    " (", round(100 * n_multi / max(nrow(dt), 1), 2), "%)",
    "; mean readSupport ", round(mean(dt$read_support), 3)
  )
  support_rows[[length(support_rows) + 1]] <- data.table(
    locus = loc$locus,
    n_unique_fragments = nrow(dt),
    n_read_support_gt1 = n_multi,
    frac_read_support_gt1 = n_multi / max(nrow(dt), 1),
    mean_read_support = mean(dt$read_support),
    max_read_support = max(dt$read_support)
  )
  bins <- seq(loc$start, loc$end, by = bin_size)
  # The last bin stops at the window end, so it can be shorter than 100 bp.
  bin_width <- pmin(bin_size, loc$end - bins + 1L)
  for (st in c("P24", "P48")) {
    for (sbt in c("pDm8", "yDm8")) {
      sub <- dt[stage == st & subtype == sbt]
      n_cells <- meta[stage == st & subtype == sbt, .N]
      # nCount_ATAC is the Seurat peak-count total, not genomic fragment depth.
      depth <- meta[stage == st & subtype == sbt, sum(nCount_ATAC)]
      cov <- numeric(length(bins))
      if (nrow(sub)) {
        start1 <- pmax(sub$start + 1L, loc$start)
        end1 <- pmin(sub$end, loc$end)
        ok <- !is.na(start1) & !is.na(end1) & end1 >= start1
        width <- loc$end - loc$start + 1L
        if (any(ok)) {
        rle <- coverage(
          IRanges(
            start = start1[ok] - loc$start + 1L,
            end = end1[ok] - loc$start + 1L
          ),
          weight = rep.int(1, sum(ok)),
          width = width
        )
        vals <- as.numeric(rle)
        bin_id <- ((seq_along(vals) - 1L) %/% bin_size) + 1L
        bin_id <- pmin(bin_id, length(bins))
        agg <- tapply(vals, bin_id, sum)
        cov[as.integer(names(agg))] <- as.numeric(agg)
        }
      }
      cov_rows[[length(cov_rows) + 1]] <- data.table(
        locus = loc$locus,
        stage = st,
        subtype = sbt,
        bin_start = bins,
        bin_width = bin_width,
        bin_end = bins + bin_width - 1L,
        bin_mid = bins + (bin_width - 1L) / 2,
        fragment_base_sum = cov,
        # cov sums unique-fragment coverage over the bases actually in the bin.
        # atac_depth is the sum of nCount_ATAC (peak counts), not genomic fragments.
        mean_per_base = cov / bin_width,
        per_million_atac = if (depth > 0) (cov / bin_width) / depth * 1e6 else NA_real_,
        n_cells = n_cells,
        atac_depth = depth
      )
    }
  }
  gene_rows[[length(gene_rows) + 1]] <- data.table(
    locus = loc$locus,
    xmin = loc$gene_start,
    xmax = loc$gene_end,
    strand = loc$strand
  )
}

cov <- rbindlist(cov_rows)
genes <- rbindlist(gene_rows)
fwrite(cov, file.path(root, "results", "locus_coverage_bins.csv"))
fwrite(rbindlist(support_rows), file.path(root, "results", "fragment_read_support.csv"))

# scale the gene bar to sit just under each panel's maximum
ymax <- cov[, .(ymax = max(per_million_atac, na.rm = TRUE)), by = .(locus, stage)]
genes <- genes[rep(seq_len(nrow(genes)), each = 2)]
genes[, stage := rep(c("P24", "P48"), length.out = .N)]
genes <- merge(genes, ymax, by = c("locus", "stage"))
genes[, y := -0.04 * ymax]
genes[, yend := y]

cov[, panel := factor(paste(locus, stage), levels = c(
  "Vsx2 P24", "Vsx2 P48", "DIP-gamma P24", "DIP-gamma P48", "beat-IIIc P24", "beat-IIIc P48"
))]
genes[, panel := factor(paste(locus, stage), levels = levels(cov$panel))]

p <- ggplot() +
  geom_segment(
    data = genes,
    aes(x = xmin / 1e6, xend = xmax / 1e6, y = y, yend = yend),
    linewidth = 1.6, color = "black", lineend = "butt"
  ) +
  geom_line(
    data = cov,
    aes(bin_mid / 1e6, per_million_atac, color = subtype),
    linewidth = 0.45
  ) +
  facet_wrap(~ panel, scales = "free", ncol = 2) +
  scale_color_manual(
    values = c(pDm8 = "#0072B2", yDm8 = "#E69F00"),
    name = "Subtype"
  ) +
  labs(
    title = "ATAC coverage at three Dm8 loci",
    subtitle = "Each unique fragment counted once. Mean per base, per million nCount_ATAC peak counts. Black bar is the gene body.",
    x = "Genomic position (Mb)",
    y = "Unique-fragment coverage per million peak counts"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.subtitle = element_text(size = 9),
    plot.margin = margin(8, 16, 8, 8)
  )

ggsave(
  file.path(root, "figures", "09_coverage_tracks.png"),
  p, width = 12, height = 8.2, dpi = 140
)
message("Wrote figures/09_coverage_tracks.png")

# IGV bedGraph is 0-based and half-open. bin_end is the last 1-based base,
# so the interval is bin_start-1 .. bin_end and does not run past the window.
chrom_of <- c("Vsx2" = "chrX", "DIP-gamma" = "chr3R", "beat-IIIc" = "chr2L")
color_of <- c("pDm8" = "0,114,178", "yDm8" = "230,159,0")
igv <- file.path(root, "igv")
dir.create(igv, showWarnings = FALSE)
for (st in c("P24", "P48")) {
  for (sbt in c("pDm8", "yDm8")) {
    sub <- cov[stage == st & subtype == sbt]
    sub[, chrom := chrom_of[locus]]
    setorder(sub, chrom, bin_start)
    path <- file.path(igv, paste0(st, "_", sbt, ".bedGraph"))
    con <- file(path, "w")
    writeLines(sprintf(
      "track type=bedGraph name=\"%s %s unique fragments per million peak counts\" visibility=full color=%s altColor=%s graphType=bar autoScale=on yLineOnOff=on",
      st, sbt, color_of[[sbt]], color_of[[sbt]]
    ), con)
    writeLines(sprintf(
      "%s\t%d\t%d\t%.6f",
      sub$chrom, sub$bin_start - 1L, sub$bin_end, sub$per_million_atac
    ), con)
    close(con)
    message("Wrote ", path)
  }
}
genome <- Sys.getenv("IGV_GENOME", unset = path.expand("~/igv/genomes/dm6/dm6.json"))
batch <- c(
  "new",
  paste("genome", genome),
  paste("load", file.path(igv, "P24_pDm8.bedGraph")),
  paste("load", file.path(igv, "P24_yDm8.bedGraph")),
  paste("load", file.path(igv, "P48_pDm8.bedGraph")),
  paste("load", file.path(igv, "P48_yDm8.bedGraph")),
  paste("load", file.path(igv, "DARs_P24.bed")),
  paste("load", file.path(igv, "DARs_P48.bed")),
  paste("load", file.path(igv, "DEG_P24_line.bed")),
  paste("load", file.path(igv, "DEG_P48_line.bed")),
  "goto chrX:5523927-5565183"
)
writeLines(batch, file.path(igv, "load_pDm8_yDm8.txt"))
message("DONE ", Sys.time())

#!/usr/bin/env Rscript
# Motif presence in differential peaks, using JASPAR 2024 insect CORE PWMs.
# motifmatchr is not installable on this R 4.3.3, so matches are Biostrings
# matchPWM calls at 80% of the maximum score on both strands.
# Background peaks are GC-matched to the differential peaks and are not significant
# in the individual-line limma-trend model.
user_lib <- path.expand("~/R/library")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(Biostrings)
  library(BSgenome.Dmelanogaster.UCSC.dm6)
  library(data.table)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
root <- Sys.getenv("PROJECT_ROOT", unset = "")
if (!nzchar(root)) {
  root <- if (length(file_arg)) normalizePath(file.path(dirname(file_arg[1]), "..")) else getwd()
}
out <- file.path(root, "results")
fig <- file.path(root, "figures")
atac <- fread(file.path(out, "de_atac_line_models.csv"))
atac <- atac[method == "trend"]
cmp <- fread(file.path(out, "focus_comparison.csv"))
vsx <- cmp[assay == "ATAC", unique(feature)]

parse_peak <- function(peak) {
  m <- regexec("^(chr[0-9A-Za-z]+)-([0-9]+)-([0-9]+)$", peak)
  r <- regmatches(peak, m)
  data.table(
    feature = peak,
    chrom = vapply(r, function(z) if (length(z) == 4) z[2] else NA_character_, ""),
    start = vapply(r, function(z) if (length(z) == 4) as.integer(z[3]) else NA_integer_, 1L),
    end = vapply(r, function(z) if (length(z) == 4) as.integer(z[4]) else NA_integer_, 1L)
  )
}

parse_jaspar <- function(path) {
  lines <- readLines(path)
  motifs <- list()
  i <- 1L
  while (i <= length(lines)) {
    if (!startsWith(lines[i], ">")) {
      i <- i + 1L
      next
    }
    header <- strsplit(sub("^>", "", lines[i]), "\t", fixed = TRUE)[[1]]
    id <- header[1]
    name <- if (length(header) > 1) header[2] else id
    i <- i + 1L
    rows <- list()
    while (i <= length(lines) && grepl("^[ACGT][[:space:]]+\\[", lines[i])) {
      letter <- substr(lines[i], 1, 1)
      nums <- as.numeric(strsplit(gsub("\\[|\\]", "", lines[i]), "[[:space:]]+")[[1]])
      rows[[letter]] <- nums[!is.na(nums)]
      i <- i + 1L
    }
    lens <- lengths(rows)
    if (!all(c("A", "C", "G", "T") %in% names(rows)) || length(unique(lens)) != 1) next
    mat <- rbind(A = rows$A, C = rows$C, G = rows$G, T = rows$T)
    if (anyNA(mat) || any(colSums(mat) == 0)) next
    target <- max(colSums(mat))
    prob <- sweep(mat, 2, colSums(mat), "/")
    scaled <- round(prob * target)
    for (j in seq_len(ncol(scaled))) {
      delta <- target - sum(scaled[, j])
      k <- which.max(prob[, j])
      scaled[k, j] <- scaled[k, j] + delta
    }
    if (any(scaled < 0) || length(unique(colSums(scaled))) != 1) next
    storage.mode(scaled) <- "integer"
    motifs[[length(motifs) + 1L]] <- list(id = id, name = name, pwm = PWM(scaled, type = "log2probratio"))
  }
  motifs
}

message("Parsing JASPAR ", Sys.time())
motifs <- parse_jaspar(file.path(root, "ref", "JASPAR2024_CORE_insects_nr_pfms.txt"))
message("Motifs ", length(motifs))

hit_fraction <- function(pwm, seqs) {
  width <- ncol(pwm)
  vapply(seq_along(seqs), function(i) {
    s <- seqs[[i]]
    if (length(s) < width) return(FALSE)
    plus <- length(matchPWM(pwm, s, min.score = "80%", with.score = FALSE)) > 0
    if (plus) return(TRUE)
    length(matchPWM(pwm, reverseComplement(s), min.score = "80%", with.score = FALSE)) > 0
  }, logical(1))
}

gc_content <- function(seqs) {
  lf <- letterFrequency(seqs, letters = c("A", "C", "G", "T"), as.prob = TRUE)
  rowSums(lf[, c("G", "C"), drop = FALSE])
}

# One GC-matched background peak for each foreground peak, without replacement.
match_background <- function(fg_id, pool_id, gc, n_per = 5L, max_bg = 2000L, tol = 0.05) {
  used <- character()
  chosen <- character()
  for (id in fg_id) {
    d <- abs(gc[pool_id] - gc[[id]])
    ord <- order(d)
    ok <- pool_id[ord][d[ord] <= tol]
    ok <- ok[!ok %in% used]
    take <- head(ok, n_per)
    if (!length(take)) next
    chosen <- c(chosen, take)
    used <- c(used, take)
    if (length(chosen) >= max_bg) break
  }
  unique(chosen)
}

set.seed(1)
message("Computing GC for tested peaks ", Sys.time())
tested_ids <- unique(atac$feature)
loc_all <- parse_peak(tested_ids)
loc_all <- loc_all[!is.na(chrom) & !is.na(start) & !is.na(end)]
gr_all <- GRanges(loc_all$chrom, IRanges(loc_all$start, loc_all$end))
seq_all <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, gr_all)
names(seq_all) <- loc_all$feature
gc <- gc_content(seq_all)
names(gc) <- names(seq_all)

tests <- list()
for (st in c("P24", "P48")) {
  sub <- atac[stage == st]
  sig <- sub[FDR < 0.05 & abs(logFC) >= 0.5]
  for (direction in c("pDm8", "yDm8")) {
    fg_id <- if (direction == "pDm8") sig[logFC > 0, feature] else sig[logFC < 0, feature]
    fg_id <- intersect(fg_id, names(gc))
    bg_pool <- setdiff(intersect(sub$feature, names(gc)), fg_id)
    bg_id <- match_background(fg_id, bg_pool, gc)
    message(st, " ", direction, " fg ", length(fg_id), " GC-matched bg ", length(bg_id))
    tests[[paste(st, direction)]] <- list(stage = st, direction = direction, fg = fg_id, bg = bg_id)
  }
}
needed <- unique(c(unlist(lapply(tests, function(z) c(z$fg, z$bg))), vsx))
loc <- parse_peak(needed)
loc <- loc[!is.na(chrom) & !is.na(start) & !is.na(end)]
message("Fetching ", nrow(loc), " peak sequences ", Sys.time())
gr <- GRanges(loc$chrom, IRanges(loc$start, loc$end))
seqs <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6, gr)
names(seqs) <- loc$feature

message("Scanning ", Sys.time())
# Scan each motif once on the union of sequences, then subset per test.
union_ids <- names(seqs)
hits <- parallel::mclapply(motifs, function(mo) {
  hit_fraction(mo$pwm, seqs)
}, mc.cores = 8)
names(hits) <- vapply(motifs, function(mo) mo$id, "")

rows <- list()
for (mo_i in seq_along(motifs)) {
  mo <- motifs[[mo_i]]
  present <- hits[[mo_i]]
  names(present) <- union_ids
  for (tz in tests) {
    fg <- present[tz$fg]
    bg <- present[tz$bg]
    fg <- fg[!is.na(fg)]
    bg <- bg[!is.na(bg)]
    a <- sum(fg)
    b <- sum(!fg)
    c <- sum(bg)
    d <- sum(!bg)
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
    rows[[length(rows) + 1L]] <- data.frame(
      stage = tz$stage,
      higher_in = tz$direction,
      motif_id = mo$id,
      motif = mo$name,
      fg_n = length(fg),
      fg_with_motif = a,
      bg_n = length(bg),
      bg_with_motif = c,
      odds_ratio = unname(ft$estimate),
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  }
}
res <- rbindlist(rows)
res[, FDR := p.adjust(p_value, method = "BH"), by = .(stage, higher_in)]
setorder(res, stage, higher_in, p_value)
fwrite(res, file.path(out, "motif_enrichment.csv"))
fwrite(res[FDR < 0.1], file.path(out, "motif_enrichment_fdr10.csv"))
message("Motifs at FDR<0.1: ", nrow(res[FDR < 0.1]))
# A presence table on the six Vsx2 peaks is not written. At an 80% score
# nearly every insect PWM hits a 500 bp peak, so the list is not a result.
message("DONE ", Sys.time())

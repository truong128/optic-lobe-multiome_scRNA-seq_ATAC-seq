#!/usr/bin/env Rscript
# Inspect GSE305940 multiome and save a Dm8-focused subset.
user_lib <- path.expand("~/R/library")
if (dir.exists(user_lib)) .libPaths(c(user_lib, .libPaths()))
suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(Matrix)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
root <- Sys.getenv("PROJECT_ROOT", unset = "")
if (!nzchar(root)) {
  root <- if (length(file_arg)) normalizePath(file.path(dirname(file_arg[1]), "..")) else getwd()
}
rds_path <- file.path(root, "data", "GSE305940_Multiome.rds")
out_rds <- file.path(root, "results", "dm8_subset.rds")
log_path <- file.path(root, "results", "01_inspect_log.txt")
meta_path <- file.path(root, "results", "01_metadata_columns.tsv")

log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
on.exit({
  sink()
  close(log_con)
}, add = TRUE)

message("START ", Sys.time())
message("RDS bytes ", file.info(rds_path)$size)
message("Seurat ", as.character(packageVersion("Seurat")),
        " Signac ", as.character(packageVersion("Signac")))

message("Reading RDS...")
obj <- readRDS(rds_path)
message("Loaded ", paste(class(obj), collapse = "/"), " at ", Sys.time())

as_seurat_list <- function(obj) {
  if (inherits(obj, "Seurat")) {
    return(list(all = obj))
  }
  if (is.list(obj)) {
    keep <- list()
    nms <- names(obj)
    if (is.null(nms)) nms <- paste0("item", seq_along(obj))
    for (i in seq_along(obj)) {
      if (inherits(obj[[i]], "Seurat")) keep[[nms[i]]] <- obj[[i]]
    }
    if (length(keep)) return(keep)
  }
  stop("Unrecognized RDS structure: ", paste(class(obj), collapse = "/"))
}

objs <- as_seurat_list(obj)
message("Seurat objects: ", paste(names(objs), collapse = ", "))

marker_genes <- c(
  "beat-IIIc", "CG15138", "CG14521", "DIP-gamma", "DIP-γ", "DIPgamma",
  "dpr11", "Hmx", "dac", "tj", "VGlut", "svp", "tup"
)

assay_has <- function(seu, genes) {
  hits <- list()
  for (a in Assays(seu)) {
    rn <- rownames(seu[[a]])
    found <- genes[genes %in% rn]
    if (length(found)) hits[[a]] <- found
  }
  hits
}

pick_cluster_col <- function(md) {
  prefs <- c(
    "celltype", "cell_type", "CellType", "annotation", "Annotation",
    "cluster_name", "ClusterName", "ident", "seurat_clusters",
    "cluster", "Cluster", "RNA_snn_res.0.8", "wsnn_res.0.8"
  )
  prefs <- intersect(prefs, colnames(md))
  # prefer a column whose values mention Dm8
  for (cn in colnames(md)) {
    x <- md[[cn]]
    if (is.character(x) || is.factor(x)) {
      if (any(grepl("Dm8", as.character(x), ignore.case = TRUE))) return(cn)
    }
  }
  if (length(prefs)) return(prefs[1])
  # fallback: factor/character with 30-400 levels
  for (cn in colnames(md)) {
    x <- md[[cn]]
    if (is.character(x) || is.factor(x)) {
      n <- length(unique(x))
      if (n >= 30 && n <= 400) return(cn)
    }
  }
  NA_character_
}

pick_stage_col <- function(md, obj_name) {
  for (cn in colnames(md)) {
    x <- md[[cn]]
    if (!(is.character(x) || is.factor(x))) next
    u <- unique(as.character(x))
    if (any(u %in% c("P0", "P24", "P48", "Adult", "adult", "L3"))) return(cn)
  }
  NA_character_
}

summaries <- list()
subsets <- list()

for (nm in names(objs)) {
  seu <- objs[[nm]]
  message("\n======== OBJECT ", nm, " ========")
  message("cells ", ncol(seu), " assays ", paste(Assays(seu), collapse = ","))
  message("default ", tryCatch(DefaultAssay(seu), error = function(e) "NA"))
  message("reductions ", paste(Reductions(seu), collapse = ","))
  md <- seu@meta.data
  message("n meta cols ", ncol(md))
  # column inventory
  inv <- data.frame(
    object = nm,
    column = colnames(md),
    class = vapply(md, function(x) paste(class(x), collapse = ","), ""),
    n_unique = vapply(md, function(x) length(unique(x)), 1L),
    example = vapply(md, function(x) paste(head(unique(as.character(x)), 8), collapse = " | "), ""),
    stringsAsFactors = FALSE
  )
  summaries[[nm]] <- inv
  write.table(inv, meta_path, sep = "\t", quote = FALSE, row.names = FALSE, append = nm != names(objs)[1], col.names = nm == names(objs)[1])

  for (i in seq_len(nrow(inv))) {
    if (inv$n_unique[i] <= 60) {
      message("COL ", inv$column[i], " (n=", inv$n_unique[i], ") ", inv$example[i])
    } else if (inv$class[i] %in% c("character", "factor")) {
      u <- unique(as.character(md[[inv$column[i]]]))
      hit <- grep("Dm8|DIP|beat|pDm|yDm", u, value = TRUE, ignore.case = TRUE)
      message("COL ", inv$column[i], " n=", inv$n_unique[i],
              if (length(hit)) paste(" MATCH", paste(head(hit, 30), collapse = " | ")) else paste(" ex", inv$example[i]))
    }
  }

  hits <- assay_has(seu, marker_genes)
  for (a in names(hits)) message("ASSAY ", a, " has ", paste(hits[[a]], collapse = ","))

  cl_col <- pick_cluster_col(md)
  st_col <- pick_stage_col(md, nm)
  message("cluster col: ", cl_col, " stage col: ", st_col)

  # RNA assay containing beat-IIIc or CG14521
  rna_assay <- NULL
  rna_gene_beat <- NULL
  rna_gene_dip <- NULL
  for (a in Assays(seu)) {
    rn <- rownames(seu[[a]])
    if (any(c("beat-IIIc", "CG15138") %in% rn) && any(c("CG14521", "DIP-gamma", "DIP-γ", "DIPgamma") %in% rn)) {
      # skip peak-like assays
      if (any(grepl("^chr[2-4X]-\\d+", head(rn, 20)))) next
      rna_assay <- a
      rna_gene_beat <- intersect(c("beat-IIIc", "CG15138"), rn)[1]
      rna_gene_dip <- intersect(c("CG14521", "DIP-gamma", "DIP-γ", "DIPgamma"), rn)[1]
      break
    }
  }
  message("RNA assay ", rna_assay, " beat ", rna_gene_beat, " dip ", rna_gene_dip)

  atac_assay <- NULL
  for (a in Assays(seu)) {
    rn <- rownames(seu[[a]])
    if (any(grepl("^chr[2-4X][-:]\\d+", head(rn, 30)))) {
      atac_assay <- a
      break
    }
  }
  message("ATAC assay ", atac_assay)

  keep_cells <- character(0)
  reason <- "none"

  if (!is.na(cl_col)) {
    labs <- as.character(md[[cl_col]])
    dm_levels <- unique(labs[grepl("Dm8", labs, ignore.case = TRUE)])
    if (length(dm_levels)) {
      keep_cells <- colnames(seu)[labs %in% dm_levels]
      reason <- paste0("annotation:", cl_col, "=", paste(dm_levels, collapse = ","))
      message("Annotated Dm8 levels: ", paste(dm_levels, collapse = " | "), " n=", length(keep_cells))
    }
  }

  if (!length(keep_cells) && !is.null(rna_assay) && !is.na(cl_col)) {
    # counts preferred
    mat <- tryCatch(
      GetAssayData(seu, assay = rna_assay, layer = "counts"),
      error = function(e) tryCatch(
        GetAssayData(seu, assay = rna_assay, layer = "data"),
        error = function(e2) GetAssayData(seu, assay = rna_assay, slot = "counts")
      )
    )
    beat <- mat[rna_gene_beat, ]
    dip <- mat[rna_gene_dip, ]
    labs <- as.character(md[[cl_col]])
    stats <- tapply(seq_along(labs), labs, function(ii) {
      c(
        n = length(ii),
        pct_beat = mean(beat[ii] > 0),
        pct_dip = mean(dip[ii] > 0),
        mean_beat = mean(beat[ii]),
        mean_dip = mean(dip[ii])
      )
    })
    st <- do.call(rbind, stats)
    st <- as.data.frame(st)
    st$cluster <- rownames(st)
    st <- st[order(-st$pct_beat * st$pct_dip, -st$n), ]
    write.csv(st, file.path(root, "results", paste0("cluster_marker_rates_", nm, ".csv")), row.names = FALSE)
    message("Top clusters by joint marker detection:")
    print(head(st, 15))
    # Dm8-like: both markers in >=25% of cells, and at least 30 cells
    cand <- st$cluster[st$pct_beat >= 0.25 & st$pct_dip >= 0.20 & st$n >= 30]
    # if none, take the single best joint cluster with n>=30
    if (!length(cand)) {
      cand <- st$cluster[which.max(st$pct_beat * st$pct_dip)]
    }
    # prefer clusters that also don't have extremely high dip in a non-joint way; keep top 3 max
    if (length(cand) > 3) cand <- cand[1:3]
    keep_cells <- colnames(seu)[labs %in% cand]
    reason <- paste0("marker_clusters:", cl_col, "=", paste(cand, collapse = ","))
    message("Keeping clusters ", paste(cand, collapse = ","), " n=", length(keep_cells))
  }

  if (!length(keep_cells)) {
    message("NO CELLS KEPT for ", nm)
    next
  }

  message("Subsetting ", length(keep_cells), " cells because ", reason)
  small <- subset(seu, cells = keep_cells)
  small$source_object <- nm
  small$dm8_keep_reason <- reason
  if (!is.na(st_col)) small$stage_used <- as.character(small@meta.data[[st_col]])
  else small$stage_used <- nm
  if (!is.na(cl_col)) small$cluster_used <- as.character(small@meta.data[[cl_col]])
  else small$cluster_used <- NA_character_

  # drop unused assays to shrink the subset
  keep_assays <- unique(c(rna_assay, atac_assay, "RNA", "ATAC", "peaks"))
  keep_assays <- intersect(keep_assays, Assays(small))
  keep_red <- intersect(Reductions(small), c("umap", "wnn.umap", "umap.rna", "umap.atac", "ref.umap"))
  if (!length(keep_red)) keep_red <- head(Reductions(small), 2)
  message("Diet assays ", paste(keep_assays, collapse = ","), " red ", paste(keep_red, collapse = ","))
  small <- DietSeurat(small, assays = keep_assays, dimreducs = keep_red)
  subsets[[nm]] <- small
  rm(small)
  gc()
}

if (!length(subsets)) stop("No Dm8 cells identified")

if (length(subsets) == 1) {
  out <- subsets[[1]]
} else {
  out <- subsets
}
saveRDS(out, out_rds)
message("Wrote ", out_rds, " bytes ", file.info(out_rds)$size)
message("DONE ", Sys.time())

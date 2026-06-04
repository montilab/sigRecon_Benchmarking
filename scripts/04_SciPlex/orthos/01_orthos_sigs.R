library(tidyverse)
library(Seurat)
library(Matrix)
library(reticulate)

Sys.setenv(RETICULATE_PYTHON_ENV = "r-orthos")
reticulate::use_condaenv("r-orthos", required = TRUE)
Sys.setenv(BASILISK_EXTERNAL_DIR = file.path(Sys.getenv("MLAB"), "personal/andrewdr/basilisk"))

library(orthos)
library(sigrecon)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
SIG_PATH <- file.path(PATH, "data/sigs/sciplex")
SAVE_PATH <- file.path(SIG_PATH, "orthos")

ORGANISM <- "Human"
CONTROL_NAME <- "Vehicle"
DRUG_COL <- "product_name"

dir.create(SAVE_PATH, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------
# Helpers
# ------------------------------------------------

clean_gene <- function(x) {
  x <- as.character(x)
  x <- stringr::str_remove(x, "\\.\\.\\..+")
  x <- stringr::str_trim(x)
  x
}

get_counts <- function(seurat_obj) {
  tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "counts"),
    error = function(e) {
      Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
    }
  )
}

collapse_count_matrix <- function(counts) {
  gene_ids <- clean_gene(rownames(counts))
  keep <- !is.na(gene_ids) & gene_ids != ""
  counts <- counts[keep, , drop = FALSE]
  gene_ids <- gene_ids[keep]

  counts <- as.matrix(counts)
  rownames(counts) <- gene_ids
  rowsum(counts, group = rownames(counts), reorder = FALSE)
}

prepare_lfc_table <- function(lfc_tbl) {
  lfc_tbl %>%
    dplyr::mutate(
      gene = clean_gene(gene),
      padj_for_score = dplyr::if_else(
        is.na(padj) | padj <= 0,
        .Machine$double.xmin,
        padj
      ),
      abs_score = abs(log2FoldChange) * -log10(padj_for_score)
    ) %>%
    dplyr::filter(!is.na(gene), gene != "", !is.na(log2FoldChange)) %>%
    dplyr::group_by(pb, gene) %>%
    dplyr::slice_max(abs_score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
}

lfc_table_to_matrix <- function(lfc_tbl, drugs) {
  wide_tbl <- lfc_tbl %>%
    dplyr::filter(pb %in% drugs) %>%
    dplyr::select(gene, pb, log2FoldChange) %>%
    tidyr::pivot_wider(
      names_from = pb,
      values_from = log2FoldChange,
      values_fill = 0
    )

  mat <- wide_tbl %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()

  mat[, drugs, drop = FALSE]
}

make_orthos_signature <- function(residual_vec, limit = 100, min_genes = 5) {
  residual_vec <- residual_vec[!is.na(residual_vec)]
  residual_vec <- residual_vec[!duplicated(names(residual_vec))]
  ranked_genes <- names(sort(residual_vec, decreasing = TRUE))
  positive_genes <- ranked_genes[residual_vec[ranked_genes] > 0]

  up <- head(positive_genes, limit)
  if (length(up) < min_genes) {
    warning(
      "Fewer than ", min_genes,
      " positive residual genes for one contrast; using the top-ranked genes."
    )
    up <- head(ranked_genes, limit)
  }

  list(
    up = up,
    up_full = ranked_genes
  )
}

run_source <- function(source, seurat_obj, source_sigs, lfc_tbl_path) {
  message("Running SciPlex orthos residualization for source: ", source)

  source_sig_up <- lapply(source_sigs, function(x) x$up)
  lfc_tbl <- prepare_lfc_table(readRDS(lfc_tbl_path))
  counts <- collapse_count_matrix(get_counts(seurat_obj))

  drugs <- intersect(names(source_sigs), unique(lfc_tbl$pb))
  if (length(drugs) == 0) {
    stop("No shared drugs for source: ", source)
  }

  ctrl_samples <- rownames(seurat_obj@meta.data)[
    seurat_obj@meta.data[[DRUG_COL]] == CONTROL_NAME
  ]
  ctrl_samples <- intersect(ctrl_samples, colnames(counts))
  if (length(ctrl_samples) == 0) {
    stop("No control samples found for source: ", source)
  }

  context_counts <- rowMeans(counts[, ctrl_samples, drop = FALSE])
  MD <- lfc_table_to_matrix(lfc_tbl, drugs)
  common_genes <- intersect(names(context_counts), rownames(MD))

  if (length(common_genes) < 1000) {
    stop(
      "Only ", length(common_genes),
      " genes overlap between source context counts and source LFCs for ", source, "."
    )
  }

  MD <- MD[common_genes, drugs, drop = FALSE]
  M <- matrix(
    round(context_counts[common_genes]),
    nrow = length(common_genes),
    ncol = length(drugs),
    dimnames = list(common_genes, drugs)
  )

  stopifnot(identical(dimnames(M), dimnames(MD)))

  dec <- orthos::decomposeVar(
    M = M,
    MD = MD,
    organism = ORGANISM,
    verbose = FALSE
  )

  residual_mat <- SummarizedExperiment::assay(dec, "RESIDUAL_CONTRASTS")

  orthos_sigs <- setNames(vector("list", ncol(residual_mat)), colnames(residual_mat))
  for (drug in colnames(residual_mat)) {
    limit <- length(source_sig_up[[drug]])
    if (is.null(limit) || is.na(limit) || limit < 1) {
      limit <- 100
    }

    orthos_sigs[[drug]] <- make_orthos_signature(
      residual_vec = residual_mat[, drug],
      limit = limit
    )
  }

  orthos_sigs <- orthos_sigs[lengths(lapply(orthos_sigs, `[[`, "up")) >= 5]

  saveRDS(orthos_sigs, file.path(SAVE_PATH, paste0(source, "_orthos_sigs.rds")))
  saveRDS(dec, file.path(SAVE_PATH, paste0(source, "_orthos_decomposeVar.rds")))

  message("Saved ", length(orthos_sigs), " ", source, " orthos residual signatures.")
}

# ------------------------------------------------
# Load and run sources
# ------------------------------------------------

a549_pb <- readRDS(file.path(PATH, "data/sci_plex/a549_filtered_pb.rds"))
k562_pb <- readRDS(file.path(PATH, "data/sci_plex/k562_filtered_pb.rds"))
mcf7_pb <- readRDS(file.path(PATH, "data/sci_plex/mcf7_filtered_pb.rds"))

run_source(
  source = "a549",
  seurat_obj = a549_pb,
  source_sigs = sigrecon::sciplex.a549,
  lfc_tbl_path = file.path(SIG_PATH, "a549_pb_deseq_tables.rds")
)

run_source(
  source = "k562",
  seurat_obj = k562_pb,
  source_sigs = sigrecon::sciplex.k562,
  lfc_tbl_path = file.path(SIG_PATH, "k562_pb_deseq_tables.rds")
)

run_source(
  source = "mcf7",
  seurat_obj = mcf7_pb,
  source_sigs = sigrecon::sciplex.mcf7,
  lfc_tbl_path = file.path(SIG_PATH, "mcf7_pb_deseq_tables.rds")
)

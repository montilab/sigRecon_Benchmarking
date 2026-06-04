library(tidyverse)
library(Matrix)
library(reticulate)

Sys.setenv(RETICULATE_PYTHON_ENV = "r-orthos")
reticulate::use_condaenv("r-orthos", required = TRUE)
Sys.setenv(BASILISK_EXTERNAL_DIR = file.path(Sys.getenv("MLAB"), "personal/andrewdr/basilisk"))

library(orthos)
library(sigrecon)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATA_PATH <- file.path(PATH, "data/perturb_seq")
SIG_PATH <- file.path(PATH, "data/sigs/perturb-seq")
SAVE_PATH <- file.path(SIG_PATH, "orthos")

ORGANISM <- "Human"
CONTROL_NAME <- "non-targeting"

dir.create(SAVE_PATH, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------
# Helpers
# ------------------------------------------------

collapse_count_matrix <- function(counts) {
  gene_ids <- as.character(rownames(counts))
  gene_ids <- stringr::str_trim(gene_ids)
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
      gene = stringr::str_trim(as.character(gene)),
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

lfc_table_to_matrix <- function(lfc_tbl, pbs, genes) {
  wide_tbl <- lfc_tbl %>%
    dplyr::filter(pb %in% pbs) %>%
    dplyr::select(gene, pb, log2FoldChange) %>%
    tidyr::pivot_wider(
      names_from = pb,
      values_from = log2FoldChange,
      values_fill = 0
    )

  observed_mat <- wide_tbl %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()

  genes <- unique(genes)
  mat <- matrix(
    0,
    nrow = length(genes),
    ncol = length(pbs),
    dimnames = list(genes, pbs)
  )

  fill_genes <- intersect(rownames(mat), rownames(observed_mat))
  mat[fill_genes, ] <- observed_mat[fill_genes, pbs, drop = FALSE]
  mat
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

load_perturbseq_counts <- function(source) {
  counts <- read.csv(
    file.path(DATA_PATH, paste0(source, "_processed_pb.csv")),
    row.names = 1
  )
  colnames(counts) <- stringr::str_replace_all(colnames(counts), "\\.", "-")

  meta <- read.csv(
    file.path(DATA_PATH, paste0(source, "_processed_pb_metadata.csv")),
    row.names = 1
  )

  if (source == "rpe1") {
    rownames(meta) <- stringr::str_replace_all(
      rownames(meta),
      "AC118549\\.",
      "AC118549-"
    )
  }

  stopifnot(all(colnames(counts) %in% rownames(meta)))
  counts <- counts[, rownames(meta), drop = FALSE]

  list(counts = collapse_count_matrix(counts), meta = meta)
}

run_source <- function(source, source_sigs, lfc_tbl_path) {
  message("Running perturb-seq orthos residualization for source: ", source)

  source_sig_up <- lapply(source_sigs, function(x) x$up)
  lfc_tbl <- prepare_lfc_table(readRDS(lfc_tbl_path))
  data <- load_perturbseq_counts(source)

  pbs <- intersect(names(source_sigs), unique(lfc_tbl$pb))
  if (length(pbs) == 0) {
    stop("No shared perturbations for source: ", source)
  }

  ctrl_samples <- rownames(data$meta)[data$meta$gene == CONTROL_NAME]
  ctrl_samples <- intersect(ctrl_samples, colnames(data$counts))
  if (length(ctrl_samples) == 0) {
    stop("No control samples found for source: ", source)
  }

  context_counts <- rowMeans(data$counts[, ctrl_samples, drop = FALSE])
  context_genes <- names(context_counts)
  MD <- lfc_table_to_matrix(lfc_tbl, pbs, genes = context_genes)
  common_genes <- intersect(context_genes, rownames(MD))
  message("Source ", source, " has ", length(context_genes), " context genes and ",
          sum(rowSums(abs(MD)) != 0), " genes with non-zero LFCs.")

  if (length(common_genes) < 1000) {
    stop(
      "Only ", length(common_genes),
      " genes overlap between source context counts and source LFCs for ", source, "."
    )
  }

  MD <- MD[common_genes, pbs, drop = FALSE]
  M <- matrix(
    round(context_counts[common_genes]),
    nrow = length(common_genes),
    ncol = length(pbs),
    dimnames = list(common_genes, pbs)
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
  for (pb in colnames(residual_mat)) {
    limit <- length(source_sig_up[[pb]])
    if (is.null(limit) || is.na(limit) || limit < 1) {
      limit <- 100
    }

    orthos_sigs[[pb]] <- make_orthos_signature(
      residual_vec = residual_mat[, pb],
      limit = limit
    )
  }

  orthos_sigs <- orthos_sigs[lengths(lapply(orthos_sigs, `[[`, "up")) >= 5]

  saveRDS(orthos_sigs, file.path(SAVE_PATH, paste0(source, "_orthos_sigs.rds")))
  saveRDS(dec, file.path(SAVE_PATH, paste0(source, "_orthos_decomposeVar.rds")))

  message("Saved ", length(orthos_sigs), " ", source, " orthos residual signatures.")
}

# ------------------------------------------------
# Run sources
# ------------------------------------------------

run_source(
  source = "k562",
  source_sigs = sigrecon::perturbseq.k562,
  lfc_tbl_path = file.path(SIG_PATH, "k562_pb_deseq_tables.rds")
)

run_source(
  source = "rpe1",
  source_sigs = sigrecon::perturbseq.rpe1,
  lfc_tbl_path = file.path(SIG_PATH, "rpe1_pb_deseq_tables.rds")
)

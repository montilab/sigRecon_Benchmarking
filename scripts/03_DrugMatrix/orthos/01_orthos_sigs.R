library(tidyverse)
library(Biobase)
library(reticulate)

Sys.setenv(RETICULATE_PYTHON_ENV = "r-orthos")
reticulate::use_condaenv("r-orthos", required = TRUE)
Sys.setenv(BASILISK_EXTERNAL_DIR = file.path(Sys.getenv("MLAB"), "personal/andrewdr/basilisk"))

library(orthos)
library(sigrecon)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
DATA_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/drugmatrix")
SIG_PATH <- file.path(PATH, "data/sigs/drugmatrix")
TABLE_PATH <- file.path(SIG_PATH, "tables")
SAVE_PATH <- file.path(SIG_PATH, "orthos")

ORGANISM <- "Mouse"
CONTROL_DOSE <- "0 mg/kg"

dir.create(SAVE_PATH, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------
# Helpers
# ------------------------------------------------

clean_gene <- function(x) {
  x <- as.character(x)
  x <- stringr::str_trim(x)
  x
}

load_lfc_tables <- function(source) {
  paths <- Sys.glob(file.path(TABLE_PATH, paste0("*_", source, "_sigs.csv")))
  if (length(paths) == 0) {
    stop("No Limma tables found for source: ", source)
  }

  tbls <- lapply(paths, function(path) {
    drug <- sub(paste0("_", source, "_sigs\\.csv$"), "", basename(path))

    read.csv(path, stringsAsFactors = FALSE) %>%
      dplyr::transmute(
        drug = drug,
        gene = clean_gene(Gene.Symbol),
        log2FoldChange = logFC,
        padj = adj.P.Val
      )
  })

  dplyr::bind_rows(tbls)
}

prepare_lfc_table <- function(lfc_tbl) {
  lfc_tbl %>%
    dplyr::mutate(
      padj_for_score = dplyr::if_else(
        is.na(padj) | padj <= 0,
        .Machine$double.xmin,
        padj
      ),
      abs_score = abs(log2FoldChange) * -log10(padj_for_score)
    ) %>%
    dplyr::filter(!is.na(gene), gene != "", !is.na(log2FoldChange)) %>%
    dplyr::group_by(drug, gene) %>%
    dplyr::slice_max(abs_score, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
}

lfc_table_to_matrix <- function(lfc_tbl, drugs) {
  wide_tbl <- lfc_tbl %>%
    dplyr::filter(drug %in% drugs) %>%
    dplyr::select(gene, drug, log2FoldChange) %>%
    tidyr::pivot_wider(
      names_from = drug,
      values_from = log2FoldChange,
      values_fill = 0
    )

  mat <- wide_tbl %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()

  mat[, drugs, drop = FALSE]
}

expression_context <- function(eset) {
  gene_symbols <- clean_gene(Biobase::fData(eset)$`Gene Symbol`)
  keep <- !is.na(gene_symbols) & gene_symbols != ""
  expr_mat <- Biobase::exprs(eset)[keep, , drop = FALSE]
  gene_symbols <- gene_symbols[keep]

  control_samples <- rownames(Biobase::pData(eset))[
    Biobase::pData(eset)$`dose:ch1` == CONTROL_DOSE
  ]
  if (length(control_samples) == 0) {
    stop("No control samples found in DrugMatrix ExpressionSet.")
  }

  # DrugMatrix expression is log-transformed microarray data, not counts.
  # Use 2^expression summed across duplicate symbols as a positive orthos context proxy.
  pseudo_counts <- 2^expr_mat[, control_samples, drop = FALSE]
  rownames(pseudo_counts) <- gene_symbols
  collapsed <- rowsum(as.matrix(pseudo_counts), group = rownames(pseudo_counts), reorder = FALSE)
  rowMeans(collapsed)
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

run_source <- function(source, eset, source_sigs) {
  message("Running DrugMatrix orthos residualization for source: ", source)

  source_sig_up <- lapply(source_sigs, function(x) x$up)
  lfc_tbl <- prepare_lfc_table(load_lfc_tables(source))

  drugs <- intersect(names(source_sigs), unique(lfc_tbl$drug))
  if (length(drugs) == 0) {
    stop("No shared drugs for source: ", source)
  }

  context_counts <- expression_context(eset)
  MD <- lfc_table_to_matrix(lfc_tbl, drugs)
  common_genes <- intersect(names(context_counts), rownames(MD))

  if (length(common_genes) < 1000) {
    stop(
      "Only ", length(common_genes),
      " genes overlap between source context and source LFCs for ", source, "."
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

liver_eset <- readRDS(file.path(DATA_PATH, "liver.rds"))
kidney_eset <- readRDS(file.path(DATA_PATH, "kidney.rds"))

run_source(
  source = "liver",
  eset = liver_eset,
  source_sigs = sigrecon::drugmatrix.liver
)

run_source(
  source = "kidney",
  eset = kidney_eset,
  source_sigs = sigrecon::drugmatrix.kidney
)

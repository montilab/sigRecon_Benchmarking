library(tidyverse)
library(Seurat)
library(Matrix)
library(orthos)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(sigrecon)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
SAVE_PATH <- file.path(PATH, "data/sigs/tahoe")

PREFERRED_SOURCE_CELL <- "NCI-H23"
CONTROL_NAME <- "DMSO_TF"
ORGANISM <- "Human"

dir.create(SAVE_PATH, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------
# Helpers
# ------------------------------------------------

map_ids <- function(genes) {
  genes <- as.character(genes)
  is_ensembl <- grepl("^ENSG", genes)

  if (any(is_ensembl)) {
    ensembl_ids <- genes[is_ensembl]
    symbols <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = ensembl_ids,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
    symbols[is.na(symbols)] <- ensembl_ids[is.na(symbols)]
    genes[is_ensembl] <- unname(symbols)
  }

  genes
}

collapse_count_matrix <- function(counts) {
  mapped_genes <- map_ids(rownames(counts))
  keep <- !is.na(mapped_genes) & mapped_genes != ""
  counts <- counts[keep, , drop = FALSE]
  mapped_genes <- mapped_genes[keep]

  counts <- as.matrix(counts)
  rownames(counts) <- mapped_genes
  rowsum(counts, group = rownames(counts), reorder = FALSE)
}

get_counts <- function(seurat_obj) {
  tryCatch(
    Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "counts"),
    error = function(e) {
      Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
    }
  )
}

resolve_source_cell <- function(tahoe_tbls, seurat_obj, celllines, preferred_source_cell) {
  available_sources <- intersect(
    unique(tahoe_tbls$cell_line),
    unique(seurat_obj$cell_name)
  )

  if (preferred_source_cell %in% available_sources) {
    return(preferred_source_cell)
  }

  if (length(celllines) > 0 && celllines[1] %in% available_sources) {
    message(
      "Preferred source cell '", preferred_source_cell,
      "' not found; using first cell_lines.csv entry '", celllines[1], "'."
    )
    return(celllines[1])
  }

  stop(
    "Could not resolve source cell line. Available shared cell lines: ",
    paste(available_sources, collapse = ", ")
  )
}

prepare_lfc_table <- function(tahoe_tbls, source_cell) {
  tahoe_tbls %>%
    dplyr::filter(cell_line == source_cell) %>%
    dplyr::filter(drug != CONTROL_NAME) %>%
    dplyr::mutate(
      gene = map_ids(gene),
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

# ------------------------------------------------
# Load inputs
# ------------------------------------------------

tahoe_tbls <- readRDS(file.path(PATH, "data/tahoe/tahoe_deseq_dfs.rds"))
seurat_obj_f <- readRDS(file.path(PATH, "data/tahoe/pseudobulk/merged_pseudobulk_filtered.rds"))

celllines <- read.csv(file.path(PATH, "data/sigs/tahoe/cell_lines.csv"))$x
targets <- celllines[2:10]
SOURCE_CELL <- resolve_source_cell(
  tahoe_tbls = tahoe_tbls,
  seurat_obj = seurat_obj_f,
  celllines = celllines,
  preferred_source_cell = PREFERRED_SOURCE_CELL
)

source_sigs <- sigrecon::tahoe.nci_h23
source_sig_up <- lapply(source_sigs, function(x) x$up)

lfc_tbl <- prepare_lfc_table(tahoe_tbls, SOURCE_CELL)
if (nrow(lfc_tbl) == 0) {
  stop("No DESeq rows found for SOURCE_CELL = ", SOURCE_CELL)
}

drugs <- Reduce(
  intersect,
  list(
    names(source_sigs),
    unique(lfc_tbl$drug)
  )
)

if (length(drugs) == 0) {
  stop("No shared drugs between NCI_H23 signatures and Tahoe DESeq LFC table.")
}

message("Running orthos decomposition for ", length(drugs), " NCI_H23 contrasts.")
message("Control benchmark targets: ", paste(targets, collapse = ", "))

# ------------------------------------------------
# Build orthos input matrices
# ------------------------------------------------

counts <- get_counts(seurat_obj_f)
counts <- collapse_count_matrix(counts)

source_ctrl_samples <- colnames(seurat_obj_f)[
  seurat_obj_f$cell_name == SOURCE_CELL &
    seurat_obj_f$drug_name == CONTROL_NAME
]

if (length(source_ctrl_samples) == 0) {
  stop("No control samples found for SOURCE_CELL = ", SOURCE_CELL)
}

source_ctrl_samples <- intersect(source_ctrl_samples, colnames(counts))
if (length(source_ctrl_samples) == 0) {
  stop("Source control samples were not found in the count matrix.")
}

context_counts <- Matrix::rowMeans(counts[, source_ctrl_samples, drop = FALSE])

MD <- lfc_table_to_matrix(lfc_tbl, drugs)
common_genes <- intersect(names(context_counts), rownames(MD))

if (length(common_genes) < 1000) {
  stop(
    "Only ", length(common_genes),
    " genes overlap between source context counts and source LFCs."
  )
}

MD <- MD[common_genes, drugs, drop = FALSE]

# orthos' MD mode expects M and MD to have identical dimensions and names.
# We use the NCI_H23 DMSO pseudobulk mean as the source context for each contrast.
M <- matrix(
  round(context_counts[common_genes]),
  nrow = length(common_genes),
  ncol = length(drugs),
  dimnames = list(common_genes, drugs)
)

stopifnot(identical(dimnames(M), dimnames(MD)))

# ------------------------------------------------
# Decompose and convert residual LFCs to signature lists
# ------------------------------------------------

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

saveRDS(
  orthos_sigs,
  file.path(SAVE_PATH, "orthos_sigs.rds")
)

saveRDS(
  dec,
  file.path(SAVE_PATH, "orthos_decomposeVar_nci_h23.rds")
)

message("Saved ", length(orthos_sigs), " orthos residual signatures to: ",
        file.path(SAVE_PATH, "orthos_sigs.rds"))

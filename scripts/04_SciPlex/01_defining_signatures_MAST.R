library(tidyverse)
library(data.table)
library(Seurat)
library(doParallel)
library("AnnotationDbi")
library("org.Hs.eg.db")
detectCores()
registerDoParallel(16)
options(Seurat.object.assay.version = 'v5')

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
# DATA_PATH <- file.path(Sys.getenv("CBM"), "sci-plex")
DATA_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/srivatsan_2019/")
# Loading Sci-Plex Data
a549 <- readRDS(file.path(DATA_PATH, "a549.rds"))
k562 <- readRDS(file.path(DATA_PATH, "k562.rds"))
mcf7 <- readRDS(file.path(DATA_PATH, "mcf7.rds"))

# Getting Groups
stopifnot(setdiff(unique(a549$product_name), unique(mcf7$product_name)) == 0)
stopifnot(setdiff(unique(a549$product_name), unique(k562$product_name)) == 0)
products <- a549$product_name %>% unique %>% as.vector
products <- products[products != "Vehicle"]
a549_ids <- list()
k562_ids <- list()
mcf7_ids <- list()

for (product in products) {
  a549_ids[[product]] <- a549[[]] %>% rownames_to_column(var= "id") %>% filter(dose == 10000, product_name == product) %>% pull(id)
  k562_ids[[product]] <- k562[[]] %>% rownames_to_column(var= "id") %>% filter(dose == 10000, product_name == product) %>% pull(id)
  mcf7_ids[[product]] <- mcf7[[]] %>% rownames_to_column(var= "id") %>% filter(dose == 10000, product_name == product) %>% pull(id)
}

# 1. Mean Signatures

## A549
control_cell_ids <- a549[[]] %>% rownames_to_column(var= "id") %>% dplyr::filter(product_name == "Vehicle") %>% pull(id)
drug_cell_ids <- a549[[]] %>% rownames_to_column(var= "id") %>% dplyr::filter(product_name != "Vehicle") %>% pull(id)
stopifnot(drug_cell_ids %in% colnames(a549))
stopifnot(control_cell_ids %in% colnames(a549))
sig_df <- Seurat::FindMarkers(a549,
                              ident.1 = drug_cell_ids,
                              ident.2 = control_cell_ids,
                              test.use = "MAST",
                              logfc.threshold = 0,
                              only.pos = FALSE)
saveRDS(sig_df, file.path(PATH, "data/sigs/sciplex/mean/a549_sig_dfs.rds"))
## K562
control_cell_ids <- k562[[]] %>% rownames_to_column(var= "id") %>% dplyr::filter(product_name == "Vehicle") %>% pull(id)
drug_cell_ids <- k562[[]] %>% rownames_to_column(var= "id") %>% dplyr::filter(product_name != "Vehicle") %>% pull(id)
stopifnot(drug_cell_ids %in% colnames(k562))
stopifnot(control_cell_ids %in% colnames(k562))
sig_df <- Seurat::FindMarkers(k562,
                              ident.1 = drug_cell_ids,
                              ident.2 = control_cell_ids,
                              test.use = "MAST",
                              logfc.threshold = 0,
                              only.pos = FALSE)
saveRDS(sig_df, file.path(PATH, "data/sigs/sciplex/mean/k562_sig_dfs.rds"))
## MCF7
control_cell_ids <- mcf7[[]] %>% rownames_to_column(var= "id") %>% dplyr::filter(product_name == "Vehicle") %>% pull(id)
drug_cell_ids <- mcf7[[]] %>% rownames_to_column(var= "id") %>% dplyr::filter(product_name != "Vehicle") %>% pull(id)
stopifnot(drug_cell_ids %in% colnames(mcf7))
stopifnot(control_cell_ids %in% colnames(mcf7))
sig_df <- Seurat::FindMarkers(mcf7,
                              ident.1 = drug_cell_ids,
                              ident.2 = control_cell_ids,
                              test.use = "MAST",
                              logfc.threshold = 0,
                              only.pos = FALSE)
saveRDS(sig_df, file.path(PATH, "data/sigs/sciplex/mean/mcf7_sig_dfs.rds"))

# 2. Drug-specific Signatures
print("Starting a549")
a549_sigs <- foreach(product = products, .combine=dplyr::bind_rows) %dopar% {
  drug_cell_ids <- a549_ids[[product]]
  if(length(drug_cell_ids) <= 30) {
    print(paste0("Skipping ", product, ". Fewer than 30 cells."))
    NULL
  } else {
    control_cell_ids <- a549[[]] %>% rownames_to_column(var= "id") %>% dplyr::filter(product_name == "Vehicle") %>% pull(id)
    stopifnot(drug_cell_ids %in% colnames(a549))
    stopifnot(control_cell_ids %in% colnames(a549))
    sig_df <- Seurat::FindMarkers(a549,
                                  ident.1 = drug_cell_ids,
                                  ident.2 = control_cell_ids,
                                  test.use = "MAST",
                                  logfc.threshold = 0,
                                  only.pos = FALSE)
    sig_df$product <- product
    sig_df
  }
}

saveRDS(a549_sigs, file.path(PATH, "data/sigs/sciplex/a549_sig_dfs.rds"))

print("Starting K562")
k562_sigs <- foreach(product = products, .combine=dplyr::bind_rows) %dopar% {
  drug_cell_ids <- k562_ids[[product]]
  if(length(drug_cell_ids) <= 30) {
    print(paste0("Skipping ", product, ". Fewer than 30 cells."))
    NULL
  } else {
    control_cell_ids <- k562[[]] %>% rownames_to_column(var= "id") %>% dplyr::filter(product_name == "Vehicle") %>% pull(id)
    stopifnot(drug_cell_ids %in% colnames(k562))
    stopifnot(control_cell_ids %in% colnames(k562))
    sig_df <- Seurat::FindMarkers(k562,
                                  ident.1 = drug_cell_ids,
                                  ident.2 = control_cell_ids,
                                  test.use = "MAST",
                                  logfc.threshold = 0,
                                  only.pos = FALSE)
    sig_df$product <- product
    sig_df
  }
}

saveRDS(k562_sigs, file.path(PATH, "data/sigs/sciplex/k562_sig_dfs.rds"))

print("Starting MCF7")
mcf7_sigs <- foreach(product = products, .combine=dplyr::bind_rows) %dopar% {
  drug_cell_ids <- mcf7_ids[[product]]
  if(length(drug_cell_ids) <= 30) {
    print(paste0("Skipping ", product, ". Fewer than 30 cells."))
    NULL
  } else {
    control_cell_ids <- mcf7[[]] %>% rownames_to_column(var= "id") %>% dplyr::filter(product_name == "Vehicle") %>% pull(id)
    stopifnot(drug_cell_ids %in% colnames(mcf7))
    stopifnot(control_cell_ids %in% colnames(mcf7))
    sig_df <- Seurat::FindMarkers(mcf7,
                                  ident.1 = drug_cell_ids,
                                  ident.2 = control_cell_ids,
                                  test.use = "MAST",
                                  logfc.threshold = 0,
                                  only.pos = FALSE)
    sig_df$product <- product
    sig_df
  }
}

saveRDS(mcf7_sigs, file.path(PATH, "data/sigs/sciplex/mcf7_sig_dfs.rds"))

# Filtering signatures
# stopifnot(all.equal(unique(mcf7_sigs$product), unique(k562_sigs$product)))
# stopifnot(all.equal(unique(a549_sigs$product), unique(k562_sigs$product)))

sig_filter_fn <- function(diff_table, products, alpha=0.05, limit=100) {
  results <- list()
  for(drug in products) {
    print(drug)
    diff_table$logFC_adjpval <- diff_table$avg_log2FC * (-log10(diff_table$p_val_adj))
    full_sig <- diff_table %>%
      dplyr::filter(product==drug) %>%
      dplyr::arrange(desc(logFC_adjpval)) %>%
      rownames
    sig_sig <- diff_table %>%
      dplyr::filter(product==drug) %>%
      dplyr::filter(avg_log2FC > 0) %>%
      dplyr::filter(p_val_adj <= alpha) %>%
      dplyr::arrange(desc(logFC_adjpval)) %>%
      dplyr::slice(1:limit) %>%
      rownames
    results[[drug]][["up"]] <- sig_sig
    results[[drug]][["up_full"]] <- full_sig
  }
  return(results)
}

# Filtering for significant signatures
a549_sigs <- readRDS(file.path(PATH, "data/sigs/sciplex/a549_sig_dfs.rds"))
k562_sigs <- readRDS(file.path(PATH, "data/sigs/sciplex/k562_sig_dfs.rds"))
mcf7_sigs <- readRDS(file.path(PATH, "data/sigs/sciplex/mcf7_sig_dfs.rds"))

mcf7_sig_list <- sig_filter_fn(mcf7_sigs, unique(mcf7_sigs$product))
k562_sig_list <- sig_filter_fn(k562_sigs, unique(k562_sigs$product))
a549_sig_list <- sig_filter_fn(a549_sigs, unique(a549_sigs$product))
# Removing extra string padding
mcf7_sig_list <- lapply(mcf7_sig_list, function(x) lapply(x, function(y) str_remove(string = y, pattern = "\\.\\.\\..+")))
k562_sig_list <- lapply(k562_sig_list, function(x) lapply(x, function(y) str_remove(string = y, pattern = "\\.\\.\\..+")))
a549_sig_list <- lapply(a549_sig_list, function(x) lapply(x, function(y) str_remove(string = y, pattern = "\\.\\.\\..+")))
# Removing NAs
mcf7_sig_list <- lapply(mcf7_sig_list, function(x) lapply(x, function(y) y[!is.na(y)]))
k562_sig_list <- lapply(k562_sig_list, function(x) lapply(x, function(y) y[!is.na(y)]))
a549_sig_list <- lapply(a549_sig_list, function(x) lapply(x, function(y) y[!is.na(y)]))

lapply(mcf7_sig_list, function(x) length(x$up)) %>% unlist %>% unname %>% fivenum
lapply(k562_sig_list, function(x) length(x$up)) %>% unlist %>% unname %>% fivenum
lapply(a549_sig_list, function(x) length(x$up)) %>% unlist %>% unname %>% fivenum
lapply(mcf7_sig_list, function(x) length(x$up_full)) %>% unlist %>% unname %>% fivenum
lapply(k562_sig_list, function(x) length(x$up_full)) %>% unlist %>% unname %>% fivenum
lapply(a549_sig_list, function(x) length(x$up_full)) %>% unlist %>% unname %>% fivenum
# Filtering to perturbations with a reasonable number of DEGS (arbitrarily set as 5)
mcf7_sig_list <- mcf7_sig_list[lapply(mcf7_sig_list, function(x) (length(x$up) >= 5)) %>% unlist]
k562_sig_list <- k562_sig_list[lapply(k562_sig_list, function(x) (length(x$up) >= 5)) %>% unlist]
a549_sig_list <- a549_sig_list[lapply(a549_sig_list, function(x) (length(x$up) >= 5)) %>% unlist]
# # Converting to hgnc symbol
# mcf7_sigs <- lapply(mcf7_sigs, function(x) lapply(x, function(y) mapIds(org.Hs.eg.db, keys=y, column="SYMBOL", keytype="ENSEMBL", multiVals="first")))
# k562_sigs <- lapply(k562_sigs, function(x) lapply(x, function(y) mapIds(org.Hs.eg.db, keys=y, column="SYMBOL", keytype="ENSEMBL", multiVals="first")))
# a549_sigs <- lapply(a549_sigs, function(x) lapply(x, function(y) mapIds(org.Hs.eg.db, keys=y, column="SYMBOL", keytype="ENSEMBL", multiVals="first")))
# # Unnaming and removing nas
# mcf7_sigs <- lapply(mcf7_sigs, function(x) lapply(x, function(y) unname(y[!is.na(y)])))
# k562_sigs <- lapply(k562_sigs, function(x) lapply(x, function(y) unname(y[!is.na(y)])))
# a549_sigs <- lapply(a549_sigs, function(x) lapply(x, function(y) unname(y[!is.na(y)])))

saveRDS(mcf7_sig_list, file.path(PATH, "data/sigs/sciplex/mcf7_sigs.rds"))
saveRDS(k562_sig_list, file.path(PATH, "data/sigs/sciplex/k562_sigs.rds"))
saveRDS(a549_sig_list, file.path(PATH, "data/sigs/sciplex/a549_sigs.rds"))

# Filtering to shared products
shared_products <- purrr::reduce(list(names(mcf7_sig_list), names(k562_sig_list), names(a549_sig_list)), intersect)
saveRDS(mcf7_sig_list[shared_products], file.path(PATH, "data/sigs/sciplex/mcf7_sigs_filtered.rds"))
saveRDS(k562_sig_list[shared_products], file.path(PATH, "data/sigs/sciplex/k562_sigs_filtered.rds"))
saveRDS(a549_sig_list[shared_products], file.path(PATH, "data/sigs/sciplex/a549_sigs_filtered.rds"))

# Mean Signatures
a549_sigs <- readRDS(file.path(PATH, "data/sigs/sciplex/mean/a549_sig_dfs.rds"))
k562_sigs <- readRDS(file.path(PATH, "data/sigs/sciplex/mean/k562_sig_dfs.rds"))
mcf7_sigs <- readRDS(file.path(PATH, "data/sigs/sciplex/mean/mcf7_sig_dfs.rds"))
a549_sigs$product <- "drug"
k562_sigs$product <- "drug"
mcf7_sigs$product <- "drug"
mcf7_sig_list <- sig_filter_fn(mcf7_sigs, unique(mcf7_sigs$product))
k562_sig_list <- sig_filter_fn(k562_sigs, unique(k562_sigs$product))
a549_sig_list <- sig_filter_fn(a549_sigs, unique(a549_sigs$product))
# Removing extra string padding
mcf7_sig_list <- lapply(mcf7_sig_list, function(x) lapply(x, function(y) str_remove(string = y, pattern = "\\.\\.\\..+")))
k562_sig_list <- lapply(k562_sig_list, function(x) lapply(x, function(y) str_remove(string = y, pattern = "\\.\\.\\..+")))
a549_sig_list <- lapply(a549_sig_list, function(x) lapply(x, function(y) str_remove(string = y, pattern = "\\.\\.\\..+")))
# Removing NAs
mcf7_sig_list <- lapply(mcf7_sig_list, function(x) lapply(x, function(y) y[!is.na(y)]))
k562_sig_list <- lapply(k562_sig_list, function(x) lapply(x, function(y) y[!is.na(y)]))
a549_sig_list <- lapply(a549_sig_list, function(x) lapply(x, function(y) y[!is.na(y)]))

# Filtering to perturbations with a reasonable number of DEGS (arbitrarily set as 5)
mcf7_sig_list <- mcf7_sig_list[lapply(mcf7_sig_list, function(x) (length(x$up) >= 5)) %>% unlist]
# k562_sig_list <- k562_sig_list[lapply(k562_sig_list, function(x) (length(x$up) >= 5)) %>% unlist]
a549_sig_list <- a549_sig_list[lapply(a549_sig_list, function(x) (length(x$up) >= 5)) %>% unlist]
saveRDS(mcf7_sig_list, file.path(PATH, "data/sigs/sciplex/mean/mcf7_sigs_filtered.rds"))
saveRDS(k562_sig_list, file.path(PATH, "data/sigs/sciplex/mean/k562_sigs_filtered.rds"))
saveRDS(a549_sig_list, file.path(PATH, "data/sigs/sciplex/mean/a549_sigs_filtered.rds"))

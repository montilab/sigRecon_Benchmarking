library(argparse)
library(Seurat)
library(tidyverse)
library(sigrecon)
library(stringdist)
library(BiocParallel)
bp <- make_bpparam(workers = 12, progress=TRUE)
register(bp)

DATA_PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/data/sigs/sciplex/stack")
SAVE_PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon/results/eval/sciplex")
SC_PATH <- file.path(Sys.getenv("AGED"), "CBMrepositoryData/perturbational_data/srivatsan_2019")

parser <- ArgumentParser()
parser$add_argument("--experiment", default="10th")

args <- parser$parse_args()
experiment <- args$experiment

SPLIT_TYPE <- experiment
print(paste0("Benchmarking ", SPLIT_TYPE))
# Signature preparation
# 1. Need to map to ensembl since ground truth was generated with ensemblIDs but predicted signatures are in gene symbols
# 2. Change name of drugs back to ground truth labels since stack slightly changes these upon generation.
split_path <- file.path(Sys.getenv("MLAB"),
                        "projects/brcameta/projects/sig_recon/data/sigs/sciplex/drug_splits.csv")
split_tbl <- read.csv(split_path)

# load seurat reference
a549 <- readRDS(file.path(SC_PATH, "a549.rds"))

ensembl_hgnc_tbl <- a549@assays$RNA@meta.features %>%
  tibble::rownames_to_column("ensembl_id")

map_hgnc_ensembl <- function(x) lapply(x, function(y) lapply(y, function(z) ensembl_hgnc_tbl$ensembl_id[match(z, ensembl_hgnc_tbl$feature_name)]))
clean_sigs <- function(x) lapply(x, function(y) lapply(y, function(z) z[!is.na(z)]))

# load predictions
a549_mcf7_sigs <- readRDS(file.path(DATA_PATH, paste0("a549_mcf7_",SPLIT_TYPE,"_sigs.rds")))
a549_k562_sigs <- readRDS(file.path(DATA_PATH, paste0("a549_k562_",SPLIT_TYPE,"_sigs.rds")))
k562_mcf7_sigs <- readRDS(file.path(DATA_PATH, paste0("k562_mcf7_",SPLIT_TYPE,"_sigs.rds")))
k562_a549_sigs <- readRDS(file.path(DATA_PATH, paste0("k562_a549_",SPLIT_TYPE,"_sigs.rds")))
mcf7_a549_sigs <- readRDS(file.path(DATA_PATH, paste0("mcf7_a549_",SPLIT_TYPE,"_sigs.rds")))
mcf7_k562_sigs <- readRDS(file.path(DATA_PATH, paste0("mcf7_k562_",SPLIT_TYPE,"_sigs.rds")))

# map to ensembl
a549_mcf7_sigs <- lapply(a549_mcf7_sigs, map_hgnc_ensembl)
a549_k562_sigs <- lapply(a549_k562_sigs, map_hgnc_ensembl)
k562_mcf7_sigs <- lapply(k562_mcf7_sigs, map_hgnc_ensembl)
k562_a549_sigs <- lapply(k562_a549_sigs, map_hgnc_ensembl)
mcf7_a549_sigs <- lapply(mcf7_a549_sigs, map_hgnc_ensembl)
mcf7_k562_sigs <- lapply(mcf7_k562_sigs, map_hgnc_ensembl) 
# renomve NAs
a549_mcf7_sigs <- lapply(a549_mcf7_sigs, clean_sigs)
a549_k562_sigs <- lapply(a549_k562_sigs, clean_sigs)
k562_mcf7_sigs <- lapply(k562_mcf7_sigs, clean_sigs)
k562_a549_sigs <- lapply(k562_a549_sigs, clean_sigs)
mcf7_a549_sigs <- lapply(mcf7_a549_sigs, clean_sigs)
mcf7_k562_sigs <- lapply(mcf7_k562_sigs, clean_sigs) 

# load true sigs
a549_true_sigs <- sigrecon::sciplex.a549
k562_true_sigs <- sigrecon::sciplex.k562
mcf7_true_sigs <- sigrecon::sciplex.mcf7
stopifnot(all.equal(names(a549_true_sigs), names(k562_true_sigs)))
stopifnot(all.equal(names(a549_true_sigs), names(mcf7_true_sigs)))

# fuzzy name matching
a549_mcf7_drug_names <- lapply(a549_mcf7_sigs, function(x) names(x)) %>% unlist %>% unique
a549_k562_drug_names <- lapply(a549_k562_sigs, function(x) names(x)) %>% unlist %>% unique
k562_mcf7_drug_names <- lapply(k562_mcf7_sigs, function(x) names(x)) %>% unlist %>% unique
k562_a549_drug_names <- lapply(k562_a549_sigs, function(x) names(x)) %>% unlist %>% unique
mcf7_a549_drug_names <- lapply(mcf7_a549_sigs, function(x) names(x)) %>% unlist %>% unique
mcf7_k562_drug_names <- lapply(mcf7_k562_sigs, function(x) names(x)) %>% unlist %>% unique
all_drugs_names <- purrr::reduce(list(a549_mcf7_drug_names, a549_k562_drug_names, k562_mcf7_drug_names, k562_a549_drug_names, mcf7_a549_drug_names, mcf7_k562_drug_names), intersect)

all_drugs_names_missing <- all_drugs_names[!(all_drugs_names %in% names(a549_true_sigs))]

matches_df <- data.frame(
  stack_name = all_drugs_names_missing,
  best_match = NA,
  distance = NA
)

for(i in seq_along(all_drugs_names_missing)){
  
  distances <- stringdist(all_drugs_names_missing[i], names(a549_true_sigs), method="lv")
  best_idx <- which.min(distances)
  
  matches_df$best_match[i] <- names(a549_true_sigs)[best_idx]
  matches_df$distance[i] <- min(distances)
}

# View(matches_df %>% arrange(distance))

matches_df <- matches_df %>% arrange(distance) %>% dplyr::slice(c(1:7,9:13,15,18,20))

rename_map <- setNames(matches_df$best_match, matches_df$stack_name)

rename_drugs <- function(split_list, rename_map){
  
  lapply(split_list, function(drug_list){
    
    old_names <- names(drug_list)
    
    idx <- match(old_names, names(rename_map))
    replace <- !is.na(idx)
    
    old_names[replace] <- rename_map[idx[replace]]
    
    names(drug_list) <- old_names
    drug_list
  })
}

a549_mcf7_sigs <- rename_drugs(a549_mcf7_sigs, rename_map)
a549_k562_sigs <- rename_drugs(a549_k562_sigs, rename_map)
k562_mcf7_sigs <- rename_drugs(k562_mcf7_sigs, rename_map)
k562_a549_sigs <- rename_drugs(k562_a549_sigs, rename_map)
mcf7_a549_sigs <- rename_drugs(mcf7_a549_sigs, rename_map)
mcf7_k562_sigs <- rename_drugs(mcf7_k562_sigs, rename_map)

# Evaluation
a549_mcf7_eval_table <- sig_eval_table(
  source_sigs = lapply(a549_true_sigs, function(x) x$up),
  pred_sigs = lapply(a549_mcf7_sigs, function(x) lapply(x, function(y) y$up)),
  true_sigs = mcf7_true_sigs,
  source = "a549_mcf7",
  splits = TRUE,
  split_file = split_path,
  split_type = SPLIT_TYPE,
  BPPARAM = bp
)

a549_k562_eval_table <- sig_eval_table(
  source_sigs = lapply(a549_true_sigs, function(x) x$up),
  pred_sigs = lapply(a549_k562_sigs, function(x) lapply(x, function(y) y$up)),
  true_sigs = k562_true_sigs,
  source = "a549_k562",
  splits = TRUE,
  split_file = split_path,
  split_type = SPLIT_TYPE,
  BPPARAM = bp
)

k562_mcf7_eval_table <- sig_eval_table(
  source_sigs = lapply(k562_true_sigs, function(x) x$up),
  pred_sigs = lapply(k562_mcf7_sigs, function(x) lapply(x, function(y) y$up)),
  true_sigs = mcf7_true_sigs,
  source = "k562_mcf7",
  splits = TRUE,
  split_file = split_path,
  split_type = SPLIT_TYPE,
  BPPARAM = bp
)

k562_a549_eval_table <- sig_eval_table(
  source_sigs = lapply(k562_true_sigs, function(x) x$up),
  pred_sigs = lapply(k562_a549_sigs, function(x) x$up),
  true_sigs = a549_true_sigs,
  source = "k562_a549",
  splits = TRUE,
  split_file = split_path,
  split_type = SPLIT_TYPE,
  BPPARAM = bp
)

mcf7_a549_eval_table <- sig_eval_table(
  source_sigs = lapply(mcf7_true_sigs, function(x) x$up),
  pred_sigs = lapply(mcf7_a549_sigs, function(x) x$up),
  true_sigs = a549_true_sigs,
  source = "mcf7_a549",
  splits = TRUE,
  split_file = split_path,
  split_type = SPLIT_TYPE,
  BPPARAM = bp
)

mcf7_k562_eval_table <- sig_eval_table(
  source_sigs = lapply(mcf7_true_sigs, function(x) x$up),
  pred_sigs = lapply(mcf7_k562_sigs, function(x) x$up),
  true_sigs = k562_true_sigs,
  source = "mcf7_k562",
  splits = TRUE,
  split_file = split_path,
  split_type = SPLIT_TYPE,
  BPPARAM = bp
)

all_eval_table <- rbind(
  a549_mcf7_eval_table,
  a549_k562_eval_table,
  k562_mcf7_eval_table,
  k562_a549_eval_table,
  mcf7_a549_eval_table,
  mcf7_k562_eval_table
)

no_change_eval_table <- readRDS("/restricted/projectnb/brcameta/projects/sig_recon/results/eval/sciplex/no_change_eval_table.rds")
combined_aggregated_df <- merge(no_change_eval_table %>% dplyr::select("source", "gene", "NES", "jacc"), 
                                all_eval_table %>% mutate(kept_alpha = kept/(kept+displaced)), 
                                by = c("source", "gene"), 
                                suffixes = c("_FALSE", "_TRUE")) %>% tibble
saveRDS(combined_aggregated_df,
        file.path(SAVE_PATH, paste0("stack_",SPLIT_TYPE,"_eval_table.rds")))

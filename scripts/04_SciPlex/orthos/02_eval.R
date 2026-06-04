library(tidyverse)
library(sigrecon)
library(BiocParallel)

bp <- make_bpparam(workers = 15, RNGseed = 123, type = "multicore")
BiocParallel::register(bp)

PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
RECON_PATH <- file.path(PATH, "data/sigs/sciplex/orthos")
SAVE_PATH <- file.path(PATH, "results/eval/sciplex")

dir.create(SAVE_PATH, recursive = TRUE, showWarnings = FALSE)

true_sigs <- list(
  a549 = sigrecon::sciplex.a549,
  k562 = sigrecon::sciplex.k562,
  mcf7 = sigrecon::sciplex.mcf7
)

source_sigs <- lapply(true_sigs, function(x) lapply(x, function(y) y$up))

orthos_sigs <- list(
  a549 = lapply(readRDS(file.path(RECON_PATH, "a549_orthos_sigs.rds")), function(x) x$up),
  k562 = lapply(readRDS(file.path(RECON_PATH, "k562_orthos_sigs.rds")), function(x) x$up),
  mcf7 = lapply(readRDS(file.path(RECON_PATH, "mcf7_orthos_sigs.rds")), function(x) x$up)
)

cell_lines <- names(true_sigs)
pairs <- expand.grid(
  source = cell_lines,
  target = cell_lines,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(source != target)

all_eval_no_change_table <- readRDS(
  file.path(SAVE_PATH, "no_change_eval_table.rds")
)

eval_tables <- list()
for (i in seq_len(nrow(pairs))) {
  src <- pairs$source[i]
  tgt <- pairs$target[i]

  eval_tables[[paste(src, tgt, sep = "_")]] <- sig_eval_table(
    source_sigs = source_sigs[[src]],
    pred_sigs = orthos_sigs[[src]],
    true_sigs = true_sigs[[tgt]],
    source = paste0(src, "_", tgt),
    target = tgt,
    BPPARAM = bp
  )
}

all_eval_table <- dplyr::bind_rows(eval_tables)

saveRDS(
  all_eval_table,
  file.path(SAVE_PATH, "ctrl_orthos_residual_raw_eval_table.rds")
)

combined_aggregated_df <- merge(
  all_eval_no_change_table %>% dplyr::select(source, gene, NES, jacc),
  all_eval_table %>% mutate(kept_alpha = kept / (kept + displaced)),
  by = c("source", "gene"),
  suffixes = c("_FALSE", "_TRUE")
) %>%
  tibble()

saveRDS(
  combined_aggregated_df,
  file.path(SAVE_PATH, "ctrl_orthos_residual_eval_table.rds")
)

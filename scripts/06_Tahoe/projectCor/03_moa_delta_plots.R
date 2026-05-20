suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(ggswarm)
  library(purrr)
  library(stringr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) default else args[[idx + 1]]
}

PROJECT_PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
EVAL_PATH <- file.path(PROJECT_PATH, "results/eval/tahoe")

metadata_path <- arg_value("--metadata", "drug_metadata.parquet")
if (!file.exists(metadata_path)) {
  metadata_path <- file.path(PROJECT_PATH, "drug_metadata.parquet")
}

output_dir <- arg_value(
  "--output-dir",
  file.path(EVAL_PATH, "projectCor_moa_delta_plots")
)
min_n <- as.integer(arg_value("--min-n", "5"))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

drug_metadata <- arrow::read_parquet(metadata_path) |>
  as_tibble() |>
  select(drug, `moa-broad`, `moa-fine`) |>
  distinct()

read_projectcor_eval <- function(regime, score) {
  readRDS(file.path(EVAL_PATH, paste0(regime, "_projectcor_", score, ".rds"))) |>
    as_tibble() |>
    mutate(
      regime = regime,
      score = score,
      delta_NES = NES_TRUE - NES_FALSE,
      delta_jacc = jacc_TRUE - jacc_FALSE
    )
}

plot_tbl <- tidyr::crossing(
  regime = c("ctrl", "10th", "90th"),
  score = c("gsva", "eigen")
) |>
  pmap_dfr(\(regime, score) read_projectcor_eval(regime, score)) |>
  left_join(drug_metadata, by = c("gene" = "drug"))

saveRDS(plot_tbl, file.path(output_dir, "projectcor_moa_delta_plot_table.rds"))

plot_moa_delta <- function(df, score_name) {
  plot_df <- bind_rows(
    df |> transmute(score, regime, gene, target, moa_type = "moa-broad", moa = `moa-broad`, delta_metric = "delta_NES", delta_value = delta_NES),
    df |> transmute(score, regime, gene, target, moa_type = "moa-broad", moa = `moa-broad`, delta_metric = "delta_jacc", delta_value = delta_jacc),
    df |> transmute(score, regime, gene, target, moa_type = "moa-fine", moa = `moa-fine`, delta_metric = "delta_NES", delta_value = delta_NES),
    df |> transmute(score, regime, gene, target, moa_type = "moa-fine", moa = `moa-fine`, delta_metric = "delta_jacc", delta_value = delta_jacc)
  ) |>
    filter(score == score_name, !is.na(moa), is.finite(delta_value)) |>
    group_by(regime, moa_type, delta_metric, moa) |>
    filter(n_distinct(gene) >= min_n) |>
    ungroup() |>
    mutate(moa = str_wrap(moa, width = 24))
  
  if (nrow(plot_df) == 0) {
    warning("No rows to plot for ", score_name)
    return(invisible(NULL))
  }
  
  plt <- ggplot(plot_df, aes(x = reorder(moa, delta_value, median), y = delta_value)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
    geom_violin(fill = "grey88", color = "grey45", linewidth = 0.35, scale = "width", trim = TRUE) +
    ggswarm::geom_quasirandom(aes(color = target), size = 0.9, alpha = 0.55, width = 0.22) +
    facet_grid(
      rows = vars(moa_type, delta_metric),
      cols = vars(regime),
      scales = "free_y",
      space = "free_y"
    ) +
    coord_flip() +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "right",
      axis.title.y = element_blank()
    ) +
    labs(
      title = paste(score_name, "projectCor delta metrics by MOA"),
      x = NULL,
      y = "delta value",
      color = "target"
    )
  
  out_file <- file.path(output_dir, paste0(score_name, "_moa_delta_metrics.png"))
  ggsave(out_file, plt, width = 16, height = 12, dpi = 300)
}

for (score_name in c("gsva", "eigen")) {
  plot_moa_delta(plot_tbl, score_name)
}

message("Saved MOA delta plots and plotting table to: ", output_dir)

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
  library(ggbeeswarm)
  library(purrr)
  library(stringr)
  library(tibble)
  library(stringdist)
})


PROJECT_PATH <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
EVAL_PATH <- file.path(PROJECT_PATH, "results/eval/tahoe")

metadata_path <- file.path(Sys.getenv("AGED"), 
                           "CBMrepositoryData/perturbational_data/tahoe/metadata/drug_metadata.parquet")

output_dir <- file.path(PROJECT_PATH, "results/plots/tahoe_projectCor_moa_delta_plots")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

drug_metadata <- arrow::read_parquet(metadata_path) |>
  as_tibble() |>
  select(drug, `moa-broad`, `moa-fine`) |>
  distinct()

test <- readRDS(file.path(EVAL_PATH, paste0("ctrl", "_projectcor_", "gsva", "_eval_table.rds")))
test90 <- readRDS(file.path(EVAL_PATH, paste0("90th", "_projectcor_", "gsva", "_eval_table.rds")))
test10 <- readRDS(file.path(EVAL_PATH, paste0("10th", "_projectcor_", "gsva", "_eval_table.rds")))
missing_names <-unique(test$gene)[!(unique(test$gene) %in% drug_metadata$drug)]
matches_df <- data.frame(
  stack_name = missing_names,
  best_match = NA,
  distance = NA
)

for (i in seq_along(missing_names)) {
  distances <- stringdist(missing_names[i], unique(drug_metadata$drug), method = "lv")
  best_idx <- which.min(distances)
  matches_df$best_match[i] <- unique(drug_metadata$drug)[best_idx]
  matches_df$distance[i] <- min(distances)
}

rename_map <- setNames(matches_df$best_match, matches_df$stack_name)
test <- test %>% dplyr::mutate(gene = if_else(gene %in% names(rename_map), rename_map[gene], gene))

read_projectcor_eval <- function(regime, score, rename_map, EVAL_PATH) {
  readRDS(file.path(EVAL_PATH, paste0(regime, "_projectcor_", score, "_eval_table.rds"))) |>
    as_tibble() |>
    mutate(
      gene = if_else(gene %in% names(rename_map), rename_map[gene], gene),
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
  pmap_dfr(\(regime, score) read_projectcor_eval(regime, score, rename_map, EVAL_PATH)) |>
  left_join(drug_metadata, by = c("gene" = "drug"))

saveRDS(plot_tbl, file.path(output_dir, "projectcor_moa_delta_plot_table.rds"))

plot_moa_delta <- function(df, score_name) {
  regime_order <- c("90th", "10th", "ctrl")
  
  plot_df <- df |>
    transmute(
      score,
      regime,
      gene,
      target,
      moa = `moa-fine`,
      delta_NES
    ) |>
    filter(score == score_name, !is.na(moa), is.finite(delta_NES)) |>
    group_by(score, regime, target, moa, gene) |>
    summarise(
      delta_NES = median(delta_NES, na.rm = TRUE),
      n_eval_rows = n(),
      .groups = "drop"
    ) |>
    group_by(regime, target) |>
    ungroup() |>
    mutate(regime = factor(regime, levels = regime_order))
  
  if (nrow(plot_df) == 0) {
    warning("No rows to plot for ", score_name)
    return(invisible(NULL))
  }
  
  walk(unique(plot_df$target), function(target_name) {
    target_df <- plot_df |> 
      filter(target == target_name)
    
    # Use 90th subset as reference for MOA ordering
    moa_order_90th <- target_df |>
      filter(regime == "90th") |>
      group_by(moa) |>
      summarise(
        median_delta_NES_90th = median(delta_NES, na.rm = TRUE),
        .groups = "drop"
      ) |>
      arrange(desc(median_delta_NES_90th)) |>
      pull(moa)
    
    # If any MOAs are absent from 90th but present in other regimes,
    # append them at the bottom alphabetically
    remaining_moas <- setdiff(unique(target_df$moa), moa_order_90th)
    moa_order <- c(moa_order_90th, sort(remaining_moas))
    
    # coord_flip() reverses the visual top-to-bottom order,
    # so use rev() to make descending median appear top-to-bottom.
    target_df <- target_df |>
      mutate(
        moa = factor(moa, levels = rev(moa_order))
      )
    
    clean_target <- target_name |>
      str_replace_all("[^[:alnum:]_-]+", "_") |>
      str_replace_all("^_+|_+$", "")
    
    plt <- ggplot(target_df, aes(x = moa, y = delta_NES)) +
      geom_hline(yintercept = 0, color = "grey70", linewidth = 0.35) +
      geom_violin(
        fill = "grey88",
        color = "grey45",
        linewidth = 0.35,
        scale = "width",
        trim = TRUE
      ) +
      ggbeeswarm::geom_quasirandom(
        color = "#2f6f9f",
        size = 0.9,
        alpha = 0.55,
        width = 0.22
      ) +
      facet_grid(. ~ regime) +
      coord_flip() +
      scale_x_discrete(
        labels = function(x) str_wrap(x, width = 24),
        drop = FALSE
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.y = element_blank()
      ) +
      labs(
        title = paste(score_name, "projectCor delta NES by moa-fine in", target_name),
        x = NULL,
        y = "delta NES"
      )
    
    out_file <- file.path(
      output_dir,
      paste0(score_name, "_", clean_target, "_moa_fine_delta_NES.png")
    )
    
    ggsave(out_file, plt, width = 12, height = 8, dpi = 300)
  })
}

for (score_name in c("gsva", "eigen")) {
  plot_moa_delta(plot_tbl, score_name)
}

message("Saved MOA delta plots and plotting table to: ", output_dir)

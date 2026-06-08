suppressPackageStartupMessages({
  library(dplyr); 
  library(ggplot2); 
  library(purrr); 
  library(stringr); 
  library(tibble)
})

parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(); i <- 1
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    if (i == length(args) || startsWith(args[[i + 1]], "--")) { out[[key]] <- TRUE; i <- i + 1
    } else { out[[key]] <- args[[i + 1]]; i <- i + 2 }
  }
  out
}

normalize_context <- function(x) str_to_lower(str_replace_all(x, "[^[:alnum:]]", ""))

clean_label <- function(path) tools::file_path_sans_ext(basename(path)) |> str_replace_all("[^[:alnum:]_\\-]+", "_")

default_arg <- function(args, key, value) if (is.null(args[[key]])) value else args[[key]]

trade_one <- function(tbl, n_sample = NULL) {
  
  context <- unique(tbl$cell_line)[1]; context_key <- unique(tbl$context_key)[1]; drug <- unique(tbl$drug)[1]
  
  trade_input <- tbl |>
    filter(is.finite(log2FoldChange), is.finite(lfcSE), is.finite(pvalue), abs(log2FoldChange) <= 10) |>
    distinct(gene, .keep_all = TRUE) |> as.data.frame()
  
  rownames(trade_input) <- make.unique(as.character(trade_input$gene))
  
  failed <- function(msg) tibble(context = context, context_key = context_key, drug = drug, n_genes_trade = nrow(trade_input), TI = NA_real_, pi_DEG = NA_real_, trade_error = msg)
  
  if (nrow(trade_input) < 10) return(failed("Fewer than 10 finite genes retained for TRADE"))
  out <- tryCatch(
    TRADEtools::TRADE(mode = "univariate", results1 = trade_input, log2FoldChange = "log2FoldChange", lfcSE = "lfcSE", pvalue = "pvalue", n_sample = n_sample, verbose = FALSE),
    error = function(e) e
  )
  
  if (inherits(out, "error")) return(failed(conditionMessage(out)))
  
  tibble(context = context, context_key = context_key, drug = drug, n_genes_trade = nrow(trade_input), TI = out$distribution_summary$transcriptome_wide_impact, pi_DEG = out$distribution_summary$Me, trade_error = NA_character_)
}

compute_trade_metrics <- function(de_rds, n_sample = NULL) {
  
  if (!requireNamespace("TRADEtools", quietly = TRUE)) stop("TRADEtools is required. Install with remotes::install_github('ajaynadig/TRADEtools').")
  
  message("Loading DE table: ", de_rds)
  
  readRDS(de_rds) |>
    as_tibble() |>
    mutate(cell_line = as.character(cell_line), drug = as.character(drug), context_key = normalize_context(cell_line)) |>
    group_split(context_key, drug) |>
    map_dfr(trade_one, n_sample = n_sample)
}

add_trade_deltas <- function(eval_tbl, trade_metrics, source_cell_line) {
  
  source_metrics <- trade_metrics |>
    filter(context_key == normalize_context(source_cell_line)) |>
    select(drug, source_trade_context = context, source_n_genes_trade = n_genes_trade, source_TI = TI, source_pi_DEG = pi_DEG, source_trade_error = trade_error)
  
  target_metrics <- trade_metrics |>
    select(target_key = context_key, drug, target_trade_context = context, target_n_genes_trade = n_genes_trade, target_TI = TI, target_pi_DEG = pi_DEG, target_trade_error = trade_error)
  
  eval_tbl |>
    mutate(delta_NES = NES_TRUE - NES_FALSE, delta_jacc = jacc_TRUE - jacc_FALSE, target_key = normalize_context(target)) |>
    left_join(source_metrics, by = c("gene" = "drug")) |>
    left_join(target_metrics, by = c("gene" = "drug", "target_key")) |>
    mutate(delta_TI = target_TI - source_TI, delta_pi_DEG = target_pi_DEG - source_pi_DEG)
}

spearman_summary <- function(df, x_col, y_col, group_col = NULL) {
  
  cor_group <- function(dat) {
    dat <- dat |> filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]]))
    if (nrow(dat) < 3 || n_distinct(dat[[x_col]]) < 2 || n_distinct(dat[[y_col]]) < 2) return(tibble(n = nrow(dat), rho = NA_real_, p_value = NA_real_))
    test <- suppressWarnings(cor.test(dat[[x_col]], dat[[y_col]], method = "spearman"))
    tibble(n = nrow(dat), rho = unname(test$estimate), p_value = test$p.value)
  }
  
  if (is.null(group_col)) cor_group(df) |> mutate(group = "pooled", .before = 1) else df |> group_by(.data[[group_col]]) |> group_modify(~ cor_group(.x)) |> ungroup() |> rename(group = all_of(group_col))
}

plot_scatter <- function(df, x_col, y_col, target_col, title, subtitle, out_file, facet = FALSE) {
  
  plot_df <- df |> filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]]))
  
  if (nrow(plot_df) == 0) return(invisible(NULL))
  
  plt <- ggplot(plot_df, aes(.data[[x_col]], .data[[y_col]])) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_point(aes(color = .data[[target_col]]), alpha = 0.72, size = 1.8) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.45) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank(), legend.position = if (facet) "none" else "right") +
    labs(x = x_col, y = y_col, title = title, subtitle = subtitle, color = target_col)
  
  if (facet) plt <- plt + facet_wrap(stats::as.formula(paste("~", target_col)), scales = "free")
  
  ggsave(out_file, plt, width = if (facet) 12 else 7.5, height = if (facet) 8 else 5.5, dpi = 300)
}

run_trade_eval_correlations <- function(eval_paths, de_rds, output_dir, source_cell_line = "NCI-H23", n_sample = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  # trade_metrics <- compute_trade_metrics(de_rds, n_sample = n_sample)
  # saveRDS(trade_metrics, file.path(output_dir, "context_trade_metrics.rds"))
  trade_metrics <- readRDS(file.path(output_dir, "context_trade_metrics.rds"))
  all_cor <- list()
  for (eval_path in eval_paths) {
    label <- clean_label(eval_path); message("Processing eval table: ", eval_path)
    joined <- readRDS(eval_path) |> as_tibble() |> add_trade_deltas(trade_metrics, source_cell_line)
    saveRDS(joined, file.path(output_dir, paste0(label, "_trade_joined.rds")))
    for (x_col in c("delta_TI", "source_TI", "target_TI")) for (y_col in c("delta_NES")) {
      pooled <- spearman_summary(joined, x_col, y_col) |> mutate(eval_table = label, x_metric = x_col, y_metric = y_col, scope = "pooled", .before = 1)
      by_target <- spearman_summary(joined, x_col, y_col, "target") |> mutate(eval_table = label, x_metric = x_col, y_metric = y_col, scope = "by_target", .before = 1)
      all_cor[[length(all_cor) + 1]] <- bind_rows(pooled, by_target)
      subtitle <- sprintf("Spearman rho = %.3f, p = %.3g, n = %s", pooled$rho, pooled$p_value, pooled$n)
      base <- file.path(output_dir, paste0(label, "_", x_col, "_vs_", y_col))
      # plot_scatter(joined, x_col, y_col, "target", paste(label, x_col, "vs", y_col), subtitle, paste0(base, "_pooled.png"))
      plot_scatter(joined, x_col, y_col, "target", paste(label, x_col, "vs", y_col), "Faceted by target cell line", paste0(base, "_by_target.png"), TRUE)
    }
  }
  cor_tbl <- bind_rows(all_cor)
  saveRDS(cor_tbl, file.path(output_dir, "trade_eval_spearman_correlations.rds"))
  invisible(list(trade_metrics = trade_metrics, correlation_summary = cor_tbl))
}

if (sys.nframe() == 0) {
  project_path <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
  args <- parse_cli_args()
  eval_dir     <- default_arg(args, "eval-dir", file.path(project_path, "results/eval/tahoe"))
  method       <- default_arg(args, "method", "projectcor_gsva")
  eval_pattern <- default_arg(args, "eval-pattern", paste0("*", method, "*_eval_table.rds"))
  eval_paths   <- Sys.glob(file.path(eval_dir, eval_pattern))
  if (length(eval_paths) == 0) eval_paths <- list.files(eval_dir, pattern = gsub("\\*", ".*", eval_pattern), full.names = TRUE)
  run_trade_eval_correlations(
    eval_paths = eval_paths,
    de_rds = default_arg(args, "source-de-rds", file.path(project_path, "data/tahoe/tahoe_deseq_dfs.rds")),
    output_dir = default_arg(args, "output-dir", file.path(project_path, "results/plots/trade_metric_correlations/tahoe_projectCor")),
    source_cell_line = default_arg(args, "source-cell-line", "NCI-H23"),
    n_sample = if (is.null(args[["n-sample"]])) NULL else as.integer(args[["n-sample"]])
  )
}

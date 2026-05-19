suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(stringr)
  library(tibble)
})

parse_cli_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    key <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1
    } else {
      out[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

arg_or <- function(args, name, default) {
  if (!is.null(args[[name]])) args[[name]] else default
}

split_arg <- function(x) {
  if (is.null(x) || identical(x, "")) {
    character(0)
  } else {
    str_split(x, ",", simplify = FALSE)[[1]] |> str_trim() |> discard(~ .x == "")
  }
}

normalize_context <- function(x) {
  str_to_lower(str_replace_all(x, "[^[:alnum:]]", ""))
}

clean_label <- function(path) {
  tools::file_path_sans_ext(basename(path)) |>
    str_replace_all("[^[:alnum:]_\\-]+", "_")
}

compute_source_trade_metrics <- function(
  source_de_rds,
  source_cell_line = "NCI-H23",
  source_cell_col = "cell_line",
  source_pert_col = "drug",
  source_gene_col = "gene",
  n_sample = NULL
) {
  if (!requireNamespace("TRADEtools", quietly = TRUE)) {
    stop(
      "TRADEtools is required. Install it with remotes::install_github('ajaynadig/TRADEtools') ",
      "or run this script in an environment where TRADEtools is already installed."
    )
  }

  message("Loading source DE table: ", source_de_rds)
  source_de <- readRDS(source_de_rds) |> as_tibble()
  required_source_cols <- c(
    source_cell_col,
    source_pert_col,
    source_gene_col,
    "log2FoldChange",
    "lfcSE",
    "pvalue"
  )
  missing_source_cols <- setdiff(required_source_cols, names(source_de))
  if (length(missing_source_cols) > 0) {
    stop("Source DE table is missing columns: ", paste(missing_source_cols, collapse = ", "))
  }

  source_de <- source_de |>
    mutate(.source_context_norm = normalize_context(.data[[source_cell_col]])) |>
    filter(.source_context_norm == normalize_context(source_cell_line)) |>
    select(-.source_context_norm)

  if (nrow(source_de) == 0) {
    stop("No source DE rows found for source cell line: ", source_cell_line)
  }

  message("Computing TRADE metrics for ", n_distinct(source_de[[source_pert_col]]), " perturbations")
  source_de |>
    group_split(.data[[source_pert_col]]) |>
    map_dfr(function(tbl) {
      perturbation <- unique(tbl[[source_pert_col]])
      if (length(perturbation) != 1 || is.na(perturbation)) {
        stop("Could not determine unique perturbation name for a source DE group.")
      }

      trade_input <- tbl |>
        filter(
          is.finite(log2FoldChange),
          is.finite(lfcSE),
          is.finite(pvalue),
          abs(log2FoldChange) <= 10
        ) |>
        distinct(.data[[source_gene_col]], .keep_all = TRUE) |>
        as.data.frame()

      rownames(trade_input) <- make.unique(as.character(trade_input[[source_gene_col]]))

      if (nrow(trade_input) < 10) {
        return(tibble(
          !!source_pert_col := perturbation,
          n_genes_trade = nrow(trade_input),
          TI = NA_real_,
          pi_DEG = NA_real_,
          trade_error = "Fewer than 10 finite genes retained for TRADE"
        ))
      }

      trade_out <- tryCatch(
        TRADEtools::TRADE(
          mode = "univariate",
          results1 = trade_input,
          log2FoldChange = "log2FoldChange",
          lfcSE = "lfcSE",
          pvalue = "pvalue",
          n_sample = n_sample,
          verbose = FALSE
        ),
        error = function(e) e
      )

      if (inherits(trade_out, "error")) {
        return(tibble(
          !!source_pert_col := perturbation,
          n_genes_trade = nrow(trade_input),
          TI = NA_real_,
          pi_DEG = NA_real_,
          trade_error = conditionMessage(trade_out)
        ))
      }

      tibble(
        !!source_pert_col := perturbation,
        n_genes_trade = nrow(trade_input),
        TI = trade_out$distribution_summary$transcriptome_wide_impact,
        pi_DEG = trade_out$distribution_summary$Me,
        trade_error = NA_character_
      )
    })
}

spearman_summary <- function(df, x_col, y_col, group_col = NULL) {
  summarize_group <- function(dat) {
    dat <- dat |>
      filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]]))
    if (nrow(dat) < 3 || n_distinct(dat[[x_col]]) < 2 || n_distinct(dat[[y_col]]) < 2) {
      return(tibble(n = nrow(dat), rho = NA_real_, p_value = NA_real_))
    }
    test <- suppressWarnings(cor.test(dat[[x_col]], dat[[y_col]], method = "spearman"))
    tibble(n = nrow(dat), rho = unname(test$estimate), p_value = test$p.value)
  }

  if (is.null(group_col) || !group_col %in% names(df)) {
    summarize_group(df) |> mutate(group = "pooled", .before = 1)
  } else {
    df |>
      group_by(.data[[group_col]]) |>
      group_modify(~ summarize_group(.x)) |>
      ungroup() |>
      rename(group = all_of(group_col))
  }
}

plot_trade_scatter <- function(
  df,
  x_col,
  y_col,
  group_col,
  eval_pert_col,
  facet,
  title,
  subtitle,
  out_file,
  label_top = 0
) {
  plot_df <- df |>
    filter(is.finite(.data[[x_col]]), is.finite(.data[[y_col]]))

  if (nrow(plot_df) == 0) {
    warning("No finite rows for plot: ", out_file)
    return(invisible(NULL))
  }

  plt <- ggplot(plot_df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_hline(yintercept = 0, color = "grey75", linewidth = 0.35) +
    geom_point(aes(color = .data[[group_col]]), alpha = 0.72, size = 1.8) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.45) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = if (facet) "none" else "right"
    ) +
    labs(
      x = x_col,
      y = y_col,
      title = title,
      subtitle = subtitle,
      color = group_col
    )

  if (label_top > 0) {
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      label_df <- plot_df |>
        slice_max(order_by = abs(.data[[y_col]]), n = label_top, with_ties = FALSE)
      plt <- plt +
        ggrepel::geom_text_repel(
          data = label_df,
          aes(label = .data[[eval_pert_col]]),
          size = 2.7,
          max.overlaps = Inf,
          show.legend = FALSE
        )
    } else {
      warning("ggrepel is not installed; skipping labels for ", out_file)
    }
  }

  if (facet) {
    plt <- plt + facet_wrap(stats::as.formula(paste("~", group_col)), scales = "free")
  }

  ggsave(out_file, plt, width = if (facet) 12 else 7.5, height = if (facet) 8 else 5.5, dpi = 300)
}

run_trade_eval_correlations <- function(
  eval_paths,
  source_de_rds,
  output_dir,
  source_cell_line = "NCI-H23",
  source_cell_col = "cell_line",
  source_pert_col = "drug",
  source_gene_col = "gene",
  eval_pert_col = "gene",
  target_col = "target",
  n_sample = NULL,
  label_top = 0
) {
  if (length(eval_paths) == 0) {
    stop("No eval RDS files supplied.")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  trade_metrics <- compute_source_trade_metrics(
    source_de_rds = source_de_rds,
    source_cell_line = source_cell_line,
    source_cell_col = source_cell_col,
    source_pert_col = source_pert_col,
    source_gene_col = source_gene_col,
    n_sample = n_sample
  )

  write.csv(trade_metrics, file.path(output_dir, "source_trade_metrics.csv"), row.names = FALSE)
  saveRDS(trade_metrics, file.path(output_dir, "source_trade_metrics.rds"))

  all_correlations <- list()
  for (eval_path in eval_paths) {
    label <- clean_label(eval_path)
    message("Processing eval table: ", eval_path)
    eval_tbl <- readRDS(eval_path) |> as_tibble()
    required_eval_cols <- c(eval_pert_col, "NES_FALSE", "NES_TRUE", "jacc_FALSE", "jacc_TRUE")
    missing_eval_cols <- setdiff(required_eval_cols, names(eval_tbl))
    if (length(missing_eval_cols) > 0) {
      stop("Eval table ", eval_path, " is missing columns: ", paste(missing_eval_cols, collapse = ", "))
    }

    if (!target_col %in% names(eval_tbl)) {
      eval_tbl[[target_col]] <- "all_targets"
    }

    joined <- eval_tbl |>
      mutate(
        delta_NES = NES_TRUE - NES_FALSE,
        delta_jacc = jacc_TRUE - jacc_FALSE
      ) |>
      left_join(
        trade_metrics,
        by = setNames(source_pert_col, eval_pert_col)
      )

    joined_out <- file.path(output_dir, paste0(label, "_trade_joined.rds"))
    csv_out <- file.path(output_dir, paste0(label, "_trade_joined.csv"))
    saveRDS(joined, joined_out)
    write.csv(joined, csv_out, row.names = FALSE)

    missing_trade <- joined |>
      filter(is.na(TI) | is.na(pi_DEG)) |>
      distinct(.data[[eval_pert_col]]) |>
      nrow()
    if (missing_trade > 0) {
      warning(label, ": ", missing_trade, " perturbations did not join to finite TRADE metrics.")
    }

    for (x_col in c("TI", "pi_DEG")) {
      for (y_col in c("delta_NES", "delta_jacc")) {
        pooled_cor <- spearman_summary(joined, x_col, y_col) |>
          mutate(eval_table = label, x_metric = x_col, y_metric = y_col, scope = "pooled", .before = 1)
        target_cor <- spearman_summary(joined, x_col, y_col, target_col) |>
          mutate(eval_table = label, x_metric = x_col, y_metric = y_col, scope = "by_target", .before = 1)
        all_correlations[[length(all_correlations) + 1]] <- bind_rows(pooled_cor, target_cor)

        subtitle <- sprintf(
          "Spearman rho = %.3f, p = %.3g, n = %s",
          pooled_cor$rho,
          pooled_cor$p_value,
          pooled_cor$n
        )

        plot_trade_scatter(
          joined,
          x_col = x_col,
          y_col = y_col,
          group_col = target_col,
          eval_pert_col = eval_pert_col,
          facet = FALSE,
          title = paste(label, x_col, "vs", y_col),
          subtitle = subtitle,
          out_file = file.path(output_dir, paste0(label, "_", x_col, "_vs_", y_col, "_pooled.png")),
          label_top = label_top
        )
        plot_trade_scatter(
          joined,
          x_col = x_col,
          y_col = y_col,
          group_col = target_col,
          eval_pert_col = eval_pert_col,
          facet = TRUE,
          title = paste(label, x_col, "vs", y_col),
          subtitle = "Faceted by target cell line; correlations are reported in the CSV summary.",
          out_file = file.path(output_dir, paste0(label, "_", x_col, "_vs_", y_col, "_by_target.png")),
          label_top = label_top
        )
      }
    }
  }

  correlation_summary <- bind_rows(all_correlations)
  write.csv(correlation_summary, file.path(output_dir, "trade_eval_spearman_correlations.csv"), row.names = FALSE)
  saveRDS(correlation_summary, file.path(output_dir, "trade_eval_spearman_correlations.rds"))

  message("Done. Outputs written to: ", output_dir)
  invisible(list(trade_metrics = trade_metrics, correlation_summary = correlation_summary))
}

if (sys.nframe() == 0) {
  project_path <- file.path(Sys.getenv("MLAB"), "projects/brcameta/projects/sig_recon")
  args <- parse_cli_args()

  source_de_rds <- arg_or(
    args,
    "source-de-rds",
    file.path(project_path, "data/tahoe/tahoe_deseq_dfs.rds")
  )
  source_cell_line <- arg_or(args, "source-cell-line", "NCI-H23")
  source_cell_col <- arg_or(args, "source-cell-col", "cell_line")
  source_pert_col <- arg_or(args, "source-pert-col", "drug")
  source_gene_col <- arg_or(args, "source-gene-col", "gene")
  eval_pert_col <- arg_or(args, "eval-pert-col", "gene")
  target_col <- arg_or(args, "target-col", "target")
  output_dir <- arg_or(
    args,
    "output-dir",
    file.path(project_path, "results/eval/tahoe/trade_metric_correlations/projectCor")
  )
  eval_dir <- arg_or(args, "eval-dir", file.path(project_path, "results/eval/tahoe"))
  eval_pattern <- arg_or(args, "eval-pattern", "^(ctrl|10th|90th)_projectcor_(gsva|eigen)\\.rds$")
  n_sample_arg <- arg_or(args, "n-sample", "")
  n_sample <- if (identical(n_sample_arg, "")) NULL else as.integer(n_sample_arg)
  label_top <- as.integer(arg_or(args, "label-top", "0"))

  eval_paths <- split_arg(args[["eval-rds"]])
  if (length(eval_paths) == 0) {
    eval_paths <- list.files(eval_dir, pattern = eval_pattern, full.names = TRUE)
  }
  if (length(eval_paths) == 0) {
    stop("No eval RDS files found. Supply --eval-rds or adjust --eval-dir/--eval-pattern.")
  }

  run_trade_eval_correlations(
    eval_paths = eval_paths,
    source_de_rds = source_de_rds,
    output_dir = output_dir,
    source_cell_line = source_cell_line,
    source_cell_col = source_cell_col,
    source_pert_col = source_pert_col,
    source_gene_col = source_gene_col,
    eval_pert_col = eval_pert_col,
    target_col = target_col,
    n_sample = n_sample,
    label_top = label_top
  )
}

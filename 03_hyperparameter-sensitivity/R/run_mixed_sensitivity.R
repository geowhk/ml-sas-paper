suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(yaml)
  library(future)
})

source("03_hyperparameter-sensitivity/R/sensitivity_helpers.R")

expand_values <- function(values) {
  unlist(values, use.names = FALSE)
}

make_grid_from_config <- function(config) {
  tidyr::expand_grid(
    proximity_type = expand_values(config$proximity_types),
    ntree = as.integer(expand_values(config$ntree)),
    mtry = as.integer(expand_values(config$mtry)),
    nodesize = as.integer(expand_values(config$nodesize)),
    rep = seq_len(as.integer(config$replications))
  ) |>
    mutate(seed = as.integer(config$seed) + row_number())
}

cleanup_legacy_reference_outputs <- function(output_dir) {
  legacy_files <- file.path(
    output_dir,
    c(
      "objects/reference_object.rds",
      "figures/local_stat_rank_stability_heatmap.png",
      "figures/top10_spatial_overlap_by_ntree.png"
    )
  )
  unlink(legacy_files[file.exists(legacy_files)])
}

run_one_condition <- function(row, inputs, config) {
  set.seed(row$seed)
  S <- compute_urf_proximity_hp(
    data = inputs$data,
    n_tree = row$ntree,
    mtry = row$mtry,
    nodesize = row$nodesize,
    proximity_type = row$proximity_type,
    workers = as.integer(config$tree_workers)
  )
  local_stat <- compute_local_stat(inputs$nb, S, ids = inputs$dgp$sf$id)
  metrics <- summarise_true_cluster_recovery(local_stat, inputs$true_cluster)

  list(
    metrics = bind_cols(
      tibble(
        dgp_id = config$dgp_id,
        variables = paste(config$variables, collapse = ","),
        proximity_type = row$proximity_type,
        ntree = row$ntree,
        mtry = row$mtry,
        nodesize = row$nodesize,
        rep = row$rep,
        seed = row$seed
      ),
      metrics
    ),
    local_stat = local_stat
  )
}

summarise_metrics <- function(results) {
  metric_cols <- c(
    "top10_true_jaccard", "top15_true_jaccard",
    "true_local_mean", "false_local_mean", "true_false_local_gap",
    "true_rank_q25", "true_rank_median", "true_rank_q75"
  )

  results |>
    group_by(dgp_id, variables, proximity_type, ntree, mtry, nodesize) |>
    summarise(
      across(
        all_of(metric_cols),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ stats::sd(.x, na.rm = TRUE),
          q25 = ~ stats::quantile(.x, 0.25, na.rm = TRUE),
          q50 = ~ stats::quantile(.x, 0.50, na.rm = TRUE),
          q75 = ~ stats::quantile(.x, 0.75, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
}

summarise_pairwise_metrics <- function(pairwise_df) {
  metric_cols <- c(
    "local_pearson", "local_spearman", "local_rmse",
    "top10_pairwise_jaccard", "top15_pairwise_jaccard"
  )

  pairwise_df |>
    group_by(comparison_group, same_proximity_type) |>
    summarise(
      n_pairs = dplyr::n(),
      across(
        all_of(metric_cols),
        list(
          mean = ~ mean(.x, na.rm = TRUE),
          sd = ~ stats::sd(.x, na.rm = TRUE),
          q25 = ~ stats::quantile(.x, 0.25, na.rm = TRUE),
          q50 = ~ stats::quantile(.x, 0.50, na.rm = TRUE),
          q75 = ~ stats::quantile(.x, 0.75, na.rm = TRUE)
        ),
        .names = "{.col}_{.fn}"
      ),
      .groups = "drop"
    )
}

save_sensitivity_plots <- function(summary_df, pairwise_df, output_dir) {
  suppressPackageStartupMessages(library(ggplot2))
  dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

  common_theme <- theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, margin = margin(b = 8)),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      strip.text = element_text(size = 9)
    )

  p1 <- summary_df |>
    ggplot(aes(factor(mtry), factor(nodesize), fill = top10_true_jaccard_mean)) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(proximity_type ~ ntree) +
    scale_fill_viridis_c(option = "C", na.value = "grey85") +
    labs(x = "mtry", y = "nodesize", fill = "Top 10% Jaccard", title = "True-cluster recovery across hyperparameters") +
    common_theme +
    guides(fill = guide_colorbar(title.position = "top", barwidth = 10, barheight = 0.5))

  p2 <- pairwise_df |>
    ggplot(aes(comparison_group, local_spearman, fill = comparison_group)) +
    geom_boxplot(outlier.alpha = 0.15, linewidth = 0.25) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = NULL, y = "Pairwise Spearman correlation", fill = NULL, title = "Pairwise local-statistic stability") +
    common_theme +
    theme(axis.text.x = element_text(angle = 15, hjust = 1)) +
    guides(fill = guide_legend(nrow = 1))

  p3 <- summary_df |>
    ggplot(aes(factor(nodesize), top10_true_jaccard_mean, color = proximity_type, group = proximity_type)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.7) +
    geom_errorbar(aes(ymin = top10_true_jaccard_q25, ymax = top10_true_jaccard_q75), width = 0.15) +
    facet_grid(mtry ~ ntree) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "nodesize", y = "Top 10% true-cluster Jaccard", color = "Proximity", title = "Path vs terminal true-cluster recovery") +
    common_theme +
    guides(color = guide_legend(nrow = 1))

  ggsave(file.path(output_dir, "figures", "true_cluster_recovery_heatmap.png"), p1, width = 12, height = 6.2, dpi = 300)
  ggsave(file.path(output_dir, "figures", "pairwise_local_stat_stability.png"), p2, width = 12, height = 6.8, dpi = 300)
  ggsave(file.path(output_dir, "figures", "path_vs_terminal_nodesize_stability.png"), p3, width = 12, height = 7.4, dpi = 300)
}

run_sensitivity <- function(config_path) {
  config <- yaml::read_yaml(config_path)
  output_dir <- if (length(commandArgs(trailingOnly=TRUE)) == 2L) commandArgs(trailingOnly=TRUE)[2] else config$output_dir
  resume <- identical(Sys.getenv("MLSAS_RESUME", "0"), "1")
  if (file.exists(output_dir) && !resume) stop("Output must not exist: ", output_dir)
  signature <- list(config=config,code=tools::md5sum(c(
    "03_hyperparameter-sensitivity/R/run_mixed_sensitivity.R",
    "03_hyperparameter-sensitivity/R/sensitivity_helpers.R",
    "02_revision-simulation/R/dgp.R","02_revision-simulation/R/mlsas_helpers.R")))
  signature_file <- file.path(output_dir,"checkpoint_config.rds")
  if (resume && (!file.exists(signature_file) || !identical(readRDS(signature_file),signature))) stop("Checkpoint code/config mismatch")
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  if (!resume) saveRDS(signature,signature_file)
  checkpoint_dir <- file.path(output_dir,"conditions")
  dir.create(checkpoint_dir,showWarnings=FALSE)
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "objects"), recursive = TRUE, showWarnings = FALSE)
  cleanup_legacy_reference_outputs(output_dir)

  condition_workers <- as.integer(Sys.getenv("MLSAS_WORKERS", as.character(config$condition_workers)))
  tree_workers <- as.integer(config$tree_workers)
  if (condition_workers > 1 && tree_workers > 1) {
    stop("condition_workers와 tree_workers는 동시에 2 이상으로 설정하지 않는다. 조건 단위 병렬화는 tree_workers: 1을 사용한다.", call. = FALSE)
  }

  worker_count <- max(condition_workers, tree_workers)
  if (worker_count > 1) {
    future::plan(future::multisession, workers = worker_count)
  } else {
    future::plan(future::sequential)
  }
  on.exit(future::plan(future::sequential), add = TRUE)

  inputs <- load_mixed_dgp43_inputs(config$dgp_id, config$variables)
  grid <- make_grid_from_config(config)
  readr::write_csv(grid, file.path(output_dir, "tables", "sensitivity_design_grid.csv"))

  message("조건 수: ", nrow(grid))
  rows <- split(grid, seq_len(nrow(grid)))
  run_row <- function(row_df) {
    row <- as.list(row_df[1, ])
    message(
      "실행: proximity=", row$proximity_type,
      ", ntree=", row$ntree,
      ", mtry=", row$mtry,
      ", nodesize=", row$nodesize,
      ", rep=", row$rep
    )
    checkpoint <- file.path(checkpoint_dir,paste0("seed_",row$seed,".rds"))
    if (file.exists(checkpoint)) {
      result<-readRDS(checkpoint)
      stopifnot(identical(as.integer(result$metrics$seed),as.integer(row$seed)))
      return(result)
    }
    result<-run_one_condition(row,inputs,config)
    temporary<-paste0(checkpoint,".partial")
    saveRDS(result,temporary)
    if(!file.rename(temporary,checkpoint))stop("Checkpoint finalization failed")
    result
  }
  # Retain the original future RNG kind even with a sequential execution plan.
  run_outputs <- future.apply::future_lapply(
      rows,
      run_row,
      future.seed = TRUE,
      future.scheduling = Inf,
      future.packages = c("dplyr", "tibble", "randomForest", "spdep", "sf")
    )

  results <- purrr::map_dfr(run_outputs, "metrics")
  pairwise_df <- pairwise_local_stability(run_outputs)
  summary_df <- summarise_metrics(results)
  pairwise_summary <- summarise_pairwise_metrics(pairwise_df)

  readr::write_csv(results, file.path(output_dir, "tables", "sensitivity_results_raw.csv"))
  readr::write_csv(summary_df, file.path(output_dir, "tables", "sensitivity_results_summary.csv"))
  readr::write_csv(pairwise_df, file.path(output_dir, "tables", "pairwise_stability_raw.csv"))
  readr::write_csv(pairwise_summary, file.path(output_dir, "tables", "pairwise_stability_summary.csv"))
  saveRDS(
    list(config = config, results = results, summary = summary_df, pairwise = pairwise_df, pairwise_summary = pairwise_summary),
    file.path(output_dir, "objects", "sensitivity_results.rds")
  )
  save_sensitivity_plots(summary_df, pairwise_df, output_dir)

  message("Sensitivity results saved under: ", output_dir)
  invisible(results)
}

if (sys.nframe() == 0L) {
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L,2L)) {
  stop("사용법: Rscript 03_hyperparameter-sensitivity/R/run_mixed_sensitivity.R <config.yml>", call. = FALSE)
}
run_sensitivity(args[[1]])

}

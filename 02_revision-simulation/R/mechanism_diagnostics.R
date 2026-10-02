suppressPackageStartupMessages({
  library(sf)
  library(tidyverse)
  library(future)
  library(future.apply)
  library(patchwork)
  library(randomForest)
  library(spdep)
})

source("02_revision-simulation/R/mlsas_helpers.R")
source("02_revision-simulation/R/dgp.R")

default_diagnostic_specs <- function() {
  tibble::tribble(
    ~dgp_id,   ~run_dir,
    "dgp3",   "02_revision-simulation/runs/production/perm99999/dgp3_ntree1000_perm99999_20260706_081857",
    "dgp3.1", "02_revision-simulation/runs/production/perm99999/dgp3.1_ntree1000_perm99999_20260707_020849",
    "dgp4",   "02_revision-simulation/runs/production/perm99999/dgp4_ntree1000_perm99999_20260707_024628",
    "dgp4.1", "02_revision-simulation/runs/production/perm99999/dgp4.1_ntree1000_perm99999_20260707_051107",
    "dgp4.2", "02_revision-simulation/runs/production/perm99999/dgp4.2_ntree1000_perm99999_20260707_070416"
  )
}

auc_rank <- function(score, truth) {
  truth <- as.logical(truth)
  ok <- is.finite(score) & !is.na(truth)
  score <- score[ok]
  truth <- truth[ok]
  n_pos <- sum(truth)
  n_neg <- sum(!truth)
  if (n_pos == 0 || n_neg == 0) {
    return(NA_real_)
  }
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[truth]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

standardized_difference <- function(x, g) {
  g <- as.logical(g)
  x1 <- x[g]
  x0 <- x[!g]
  pooled_sd <- sqrt((stats::var(x1, na.rm = TRUE) + stats::var(x0, na.rm = TRUE)) / 2)
  if (!is.finite(pooled_sd) || pooled_sd == 0) {
    return(NA_real_)
  }
  (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / pooled_sd
}

make_pair_table <- function(S, true_cluster, nb_obj) {
  n <- nrow(S)
  upper <- which(upper.tri(S), arr.ind = TRUE)
  neighbor_key <- purrr::imap_dfr(
    nb_obj,
    function(neighbors, i) {
      neighbors <- neighbors[neighbors > i]
      if (!length(neighbors)) {
        return(tibble::tibble(i = integer(), j = integer()))
      }
      tibble::tibble(i = as.integer(i), j = as.integer(neighbors))
    }
  )
  neighbor_id <- paste(neighbor_key$i, neighbor_key$j, sep = "_")

  pair_df <- tibble::tibble(
    i = upper[, "row"],
    j = upper[, "col"],
    similarity = S[upper],
    i_cluster = true_cluster[i],
    j_cluster = true_cluster[j],
    neighbor_pair = paste(i, j, sep = "_") %in% neighbor_id,
    pair_type = dplyr::case_when(
      i_cluster & j_cluster ~ "cluster-cluster",
      xor(i_cluster, j_cluster) ~ "cluster-background",
      TRUE ~ "background-background"
    )
  )

  pair_df |>
    dplyr::mutate(
      pair_type = factor(
        pair_type,
        levels = c("cluster-cluster", "cluster-background", "background-background")
      )
    )
}

summarise_pair_similarity <- function(pair_df, dgp_id, set_name) {
  pair_df |>
    group_by(pair_type, neighbor_pair) |>
    summarise(
      n_pairs = n(),
      mean_similarity = mean(similarity, na.rm = TRUE),
      median_similarity = median(similarity, na.rm = TRUE),
      q25_similarity = quantile(similarity, 0.25, na.rm = TRUE),
      q75_similarity = quantile(similarity, 0.75, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(dgp_id = dgp_id, set_name = set_name) |>
    relocate(dgp_id, set_name)
}

summarise_local_results <- function(result_sf, dgp_id, set_name) {
  df <- sf::st_drop_geometry(result_sf) |>
    mutate(true_cluster = as.logical(true_cluster))

  tibble::tibble(
    dgp_id = dgp_id,
    set_name = set_name,
    n_true = sum(df$true_cluster, na.rm = TRUE),
    n_background = sum(!df$true_cluster, na.rm = TRUE),
    local_stat_auc = auc_rank(df$local.Moran, df$true_cluster),
    local_stat_std_diff = standardized_difference(df$local.Moran, df$true_cluster),
    raw_true = sum(df$p_perm <= 0.05 & df$true_cluster, na.rm = TRUE),
    raw_false = sum(df$p_perm <= 0.05 & !df$true_cluster, na.rm = TRUE),
    fdr_true = sum(df$p_perm_adj <= 0.05 & df$true_cluster, na.rm = TRUE),
    fdr_false = sum(df$p_perm_adj <= 0.05 & !df$true_cluster, na.rm = TRUE),
    min_p_true = min(df$p_perm[df$true_cluster], na.rm = TRUE),
    min_p_background = min(df$p_perm[!df$true_cluster], na.rm = TRUE),
    min_p_adj_true = min(df$p_perm_adj[df$true_cluster], na.rm = TRUE),
    min_p_adj_background = min(df$p_perm_adj[!df$true_cluster], na.rm = TRUE)
  )
}

make_fdr_ladder <- function(result_sf, dgp_id, set_name, alpha = 0.05) {
  df <- sf::st_drop_geometry(result_sf) |>
    mutate(true_cluster = as.logical(true_cluster)) |>
    arrange(p_perm) |>
    mutate(
      rank = row_number(),
      n_tests = n(),
      bh_threshold = rank / n_tests * alpha,
      bh_pass = p_perm <= bh_threshold,
      dgp_id = dgp_id,
      set_name = set_name
    )

  df |>
    select(dgp_id, set_name, id, true_cluster, rank, n_tests, p_perm, p_perm_adj, bh_threshold, bh_pass)
}

summarise_mechanism_contrast <- function(pair_summary, local_summary, fdr_ladder) {
  proximity_gap <- pair_summary |>
    filter(
      !neighbor_pair,
      pair_type %in% c("cluster-cluster", "cluster-background", "background-background")
    ) |>
    select(dgp_id, set_name, pair_type, mean_similarity) |>
    tidyr::pivot_wider(names_from = pair_type, values_from = mean_similarity) |>
    mutate(
      proximity_gap_cluster_vs_boundary = `cluster-cluster` - `cluster-background`,
      proximity_gap_cluster_vs_background = `cluster-cluster` - `background-background`
    ) |>
    select(dgp_id, set_name, starts_with("proximity_gap"))

  fdr_rank_summary <- fdr_ladder |>
    group_by(dgp_id, set_name) |>
    mutate(max_bh_pass_rank = ifelse(any(bh_pass), max(rank[bh_pass], na.rm = TRUE), NA_integer_)) |>
    summarise(
      true_in_top20_p = sum(true_cluster[rank <= 20], na.rm = TRUE),
      true_in_top40_p = sum(true_cluster[rank <= 40], na.rm = TRUE),
      bh_pass_true = sum(rank <= max_bh_pass_rank & true_cluster, na.rm = TRUE),
      bh_pass_false = sum(rank <= max_bh_pass_rank & !true_cluster, na.rm = TRUE),
      max_bh_pass_rank = dplyr::first(max_bh_pass_rank),
      .groups = "drop"
    )

  local_summary |>
    select(
      dgp_id, set_name, local_stat_auc, local_stat_std_diff,
      raw_true, raw_false, fdr_true, fdr_false,
      min_p_true, min_p_background, min_p_adj_true, min_p_adj_background
    ) |>
    left_join(proximity_gap, by = c("dgp_id", "set_name")) |>
    left_join(fdr_rank_summary, by = c("dgp_id", "set_name")) |>
    arrange(dgp_id, set_name)
}

diagnose_one_set <- function(analysis_result, result_obj, diagnostic_ntree, workers) {
  dgp_id <- analysis_result$dgp$id
  set_name <- result_obj$set_name
  variables <- result_obj$variables
  sf_obj <- analysis_result$dgp$sf
  true_cluster <- as.logical(sf_obj$true_cluster)

  message("진단 proximity 계산: ", dgp_id, " / ", set_name, " / ", paste(variables, collapse = ", "))
  data_df <- prepare_attribute_data(sf_obj, variables)
  S <- compute_path_proximity(
    data = data_df,
    n_tree = diagnostic_ntree,
    workers = workers,
    validate = FALSE
  )

  nb_obj <- make_nb(sf_obj, queen = TRUE)
  pair_df <- make_pair_table(S, true_cluster, nb_obj) |>
    mutate(dgp_id = dgp_id, set_name = set_name)

  list(
    pair_summary = summarise_pair_similarity(pair_df, dgp_id, set_name),
    pair_sample = pair_df |>
      group_by(pair_type, neighbor_pair) |>
      group_modify(~ dplyr::slice_sample(.x, n = min(2000, nrow(.x)), replace = FALSE)) |>
      ungroup(),
    local_summary = summarise_local_results(result_obj$sf, dgp_id, set_name),
    fdr_ladder = make_fdr_ladder(result_obj$sf, dgp_id, set_name),
    mds = {
      dist_obj <- stats::as.dist(1 - S)
      coords <- stats::cmdscale(dist_obj, k = 2)
      tibble::tibble(
        dgp_id = dgp_id,
        set_name = set_name,
        id = sf_obj$id,
        true_cluster = true_cluster,
        dim1 = coords[, 1],
        dim2 = coords[, 2]
      )
    }
  )
}

save_diagnostic_plots <- function(pair_summary, pair_sample, local_long, fdr_ladder, mds_df, output_dir) {
  dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

  p1 <- pair_summary |>
    filter(!neighbor_pair) |>
    ggplot(aes(pair_type, mean_similarity, fill = pair_type)) +
    geom_col(width = 0.75, show.legend = FALSE) +
    facet_grid(dgp_id ~ set_name, scales = "free_y") +
    labs(x = NULL, y = "Mean proximity", title = "전체 pair proximity coherence") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  p2 <- pair_summary |>
    filter(neighbor_pair) |>
    ggplot(aes(pair_type, mean_similarity, fill = pair_type)) +
    geom_col(width = 0.75, show.legend = FALSE) +
    facet_grid(dgp_id ~ set_name, scales = "free_y") +
    labs(x = NULL, y = "Mean proximity", title = "공간 이웃 pair proximity coherence") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  p3 <- local_long |>
    ggplot(aes(true_cluster, value, fill = true_cluster)) +
    geom_boxplot(outlier.alpha = 0.25, width = 0.65, show.legend = FALSE) +
    facet_grid(metric ~ dgp_id + set_name, scales = "free_y") +
    labs(x = NULL, y = NULL, title = "Local ML-SAS statistic/p-value 분포") +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))

  p4 <- fdr_ladder |>
    ggplot(aes(rank, p_perm)) +
    geom_line(aes(y = bh_threshold), color = "grey35", linewidth = 0.4) +
    geom_point(aes(color = true_cluster), alpha = 0.75, size = 0.9) +
    facet_grid(dgp_id ~ set_name, scales = "free_y") +
    labs(x = "p-value rank", y = "Permutation p-value", title = "FDR ladder: p_perm vs BH threshold") +
    theme_bw(base_size = 10)

  p5 <- mds_df |>
    ggplot(aes(dim1, dim2, color = true_cluster)) +
    geom_point(alpha = 0.8, size = 1.3) +
    facet_grid(dgp_id ~ set_name, scales = "free") +
    labs(x = "MDS 1", y = "MDS 2", title = "Proximity matrix MDS") +
    theme_bw(base_size = 10)

  ggsave(file.path(output_dir, "figures", "similarity_coherence_all_pairs.png"), p1, width = 11, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "figures", "similarity_coherence_neighbor_pairs.png"), p2, width = 11, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "figures", "local_stat_distribution.png"), p3, width = 14, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "figures", "fdr_ladder.png"), p4, width = 11, height = 8, dpi = 300)
  ggsave(file.path(output_dir, "figures", "proximity_mds.png"), p5, width = 11, height = 8, dpi = 300)

  invisible(NULL)
}

run_mechanism_diagnostics <- function(
    specs = default_diagnostic_specs(),
    output_dir = "02_revision-simulation/runs/diagnostics/mechanism_perm99999",
    diagnostic_ntree = 300,
    workers = 4
) {
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "objects"), recursive = TRUE, showWarnings = FALSE)

  if (workers > 1) {
    future::plan(future::multisession, workers = workers)
  } else {
    future::plan(future::sequential)
  }
  on.exit(future::plan(future::sequential), add = TRUE)

  results <- list()
  for (i in seq_len(nrow(specs))) {
    analysis_path <- file.path(specs$run_dir[[i]], "objects", paste0(specs$dgp_id[[i]], "_analysis_result.rds"))
    if (!file.exists(analysis_path)) {
      warning("분석 결과 객체가 없어 건너뜀: ", analysis_path)
      next
    }
    analysis_result <- readRDS(analysis_path)
    for (result_obj in analysis_result$mlsas_results) {
      key <- paste(analysis_result$dgp$id, result_obj$set_name, sep = "__")
      results[[key]] <- diagnose_one_set(
        analysis_result = analysis_result,
        result_obj = result_obj,
        diagnostic_ntree = diagnostic_ntree,
        workers = workers
      )
    }
  }

  pair_summary <- purrr::map_dfr(results, "pair_summary")
  pair_sample <- purrr::map_dfr(results, "pair_sample")
  local_summary <- purrr::map_dfr(results, "local_summary")
  fdr_ladder <- purrr::map_dfr(results, "fdr_ladder")
  mds_df <- purrr::map_dfr(results, "mds")
  mechanism_contrast <- summarise_mechanism_contrast(pair_summary, local_summary, fdr_ladder)

  local_long <- purrr::imap_dfr(
    results,
    function(res, nm) {
      parts <- strsplit(nm, "__", fixed = TRUE)[[1]]
      res_obj <- readRDS(file.path(
        specs$run_dir[match(parts[[1]], specs$dgp_id)],
        "objects",
        paste0(parts[[1]], "_analysis_result.rds")
      ))
      set_obj <- res_obj$mlsas_results[[parts[[2]]]]
      sf::st_drop_geometry(set_obj$sf) |>
        transmute(
          dgp_id = parts[[1]],
          set_name = parts[[2]],
          true_cluster = as.logical(true_cluster),
          local_stat = local.Moran,
          p_perm = p_perm,
          p_perm_adj = p_perm_adj
        ) |>
        pivot_longer(c(local_stat, p_perm, p_perm_adj), names_to = "metric", values_to = "value")
    }
  )

  readr::write_csv(pair_summary, file.path(output_dir, "tables", "similarity_pair_summary.csv"))
  readr::write_csv(pair_sample, file.path(output_dir, "tables", "similarity_pair_sample.csv"))
  readr::write_csv(local_summary, file.path(output_dir, "tables", "local_stat_summary.csv"))
  readr::write_csv(fdr_ladder, file.path(output_dir, "tables", "fdr_ladder.csv"))
  readr::write_csv(mds_df, file.path(output_dir, "tables", "proximity_mds.csv"))
  readr::write_csv(mechanism_contrast, file.path(output_dir, "tables", "mechanism_contrast_summary.csv"))
  saveRDS(results, file.path(output_dir, "objects", "mechanism_diagnostic_objects.rds"))

  save_diagnostic_plots(
    pair_summary = pair_summary,
    pair_sample = pair_sample,
    local_long = local_long,
    fdr_ladder = fdr_ladder,
    mds_df = mds_df,
    output_dir = output_dir
  )

  message("Mechanism diagnostics saved under: ", output_dir)
  invisible(list(
    pair_summary = pair_summary,
    local_summary = local_summary,
    mechanism_contrast = mechanism_contrast,
    fdr_ladder = fdr_ladder,
    mds = mds_df,
    output_dir = output_dir
  ))
}

if (identical(environment(), globalenv())) {
  args <- commandArgs(trailingOnly = TRUE)
  diagnostic_ntree <- 300
  workers <- 4
  output_dir <- "02_revision-simulation/runs/diagnostics/mechanism_perm99999"

  for (arg in args) {
    if (startsWith(arg, "--ntree=")) {
      diagnostic_ntree <- as.integer(sub("^--ntree=", "", arg))
    }
    if (startsWith(arg, "--workers=")) {
      workers <- as.integer(sub("^--workers=", "", arg))
    }
    if (startsWith(arg, "--output_dir=")) {
      output_dir <- sub("^--output_dir=", "", arg)
    }
  }

  run_mechanism_diagnostics(
    diagnostic_ntree = diagnostic_ntree,
    workers = workers,
    output_dir = output_dir
  )
}

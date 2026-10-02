# Post-run corrected production helpers for the professor-source repeated analysis.
# Rank AUC is oriented so that larger scores always indicate stronger target evidence.

suppressPackageStartupMessages({
  library(sf)
  library(spdep)
  library(randomForest)
  library(Matrix)
  library(future)
  library(future.apply)
  library(dplyr)
  library(tibble)
  library(forcats)
  library(yaml)
  library(rgeoda)
})

resolve_production_paths <- function(script_file) {
  protocol_dir <- normalizePath(file.path(dirname(script_file), ".."))
  repo_dir <- normalizePath(file.path(protocol_dir, "../.."))
  list(
    production_dir = protocol_dir,
    revision_dir = repo_dir,
    repo_dir = repo_dir,
    freeze_dir = protocol_dir,
    protocol_dir = protocol_dir
  )
}

load_production_sources <- function(paths) {
  source(file.path(paths$protocol_dir, "R", "generate_professor_dgps.R"), local = .GlobalEnv)
  source(file.path(paths$repo_dir, "02_revision-simulation", "R", "mlsas_helpers.R"), local = .GlobalEnv)
  source(file.path(paths$repo_dir, "02_revision-simulation", "R", "mlsas_optimized.R"), local = .GlobalEnv)
  source(file.path(paths$repo_dir, "02_revision-simulation", "R", "multigeary_helpers.R"), local = .GlobalEnv)
  spec <- yaml::read_yaml(file.path(paths$protocol_dir, "repeated_simulation_spec.yml"))
  if (as.integer(spec$inference$permutations) != 99999L) {
    stop("Production specification must use 99,999 permutations", call. = FALSE)
  }
  spec
}

expected_methods <- function(dgp) {
  switch(
    as.character(as.integer(dgp)),
    `2` = c("path_continuous", "path_mixed"),
    `3` = c("local_geary", "path_continuous"),
    `4` = c("local_geary", "path_continuous", "path_mixed", "terminal_mixed"),
    stop("Unsupported DGP", call. = FALSE)
  )
}

prepare_method_data <- function(sf_obj, variables) {
  out <- as.data.frame(sf::st_drop_geometry(sf_obj)[, variables, drop = FALSE])
  rownames(out) <- as.character(sf_obj$id)
  out
}

compute_terminal_proximity <- function(data, n_tree, mtry, nodesize) {
  n <- nrow(data)
  synthetic <- as.data.frame(lapply(data, function(column) sample(column, replace = TRUE)))
  data_all <- rbind(data, synthetic)
  y_all <- factor(c(rep("real", n), rep("synthetic", n)))
  rf <- randomForest::randomForest(
    x = data_all,
    y = y_all,
    ntree = as.integer(n_tree),
    mtry = as.integer(mtry),
    nodesize = as.integer(nodesize)
  )
  nodes <- attr(stats::predict(rf, data_all, nodes = TRUE), "nodes")
  nodes <- nodes[seq_len(n), , drop = FALSE]
  proximity <- matrix(0, nrow = n, ncol = n)
  for (tree_num in seq_len(n_tree)) {
    proximity <- proximity + outer(nodes[, tree_num], nodes[, tree_num], `==`)
  }
  proximity <- proximity / n_tree
  rownames(proximity) <- rownames(data)
  colnames(proximity) <- rownames(data)
  proximity
}

safe_rank_auc <- function(score, truth) {
  truth <- as.logical(truth)
  if (!any(truth) || !any(!truth)) return(NA_real_)
  positive <- score[truth]
  negative <- score[!truth]
  if (!all(is.finite(c(positive, negative)))) return(NA_real_)
  ranks <- rank(c(positive, negative), ties.method = "average")
  n_positive <- length(positive)
  n_negative <- length(negative)
  (sum(ranks[seq_len(n_positive)]) - n_positive * (n_positive + 1) / 2) /
    (n_positive * n_negative)
}

orient_rank_auc_score <- function(statistic, stronger_target_evidence = c("higher", "lower")) {
  stronger_target_evidence <- match.arg(stronger_target_evidence)
  if (identical(stronger_target_evidence, "higher")) statistic else -statistic
}

make_method_metrics <- function(local, truth, dgp, replication_id, method, score) {
  truth <- as.logical(truth)
  detected <- as.logical(local$detected)
  discoveries <- sum(detected)
  true_positive <- sum(detected & truth)
  false_positive <- sum(detected & !truth)
  false_negative <- sum(!detected & truth)
  true_count <- sum(truth)
  background_count <- sum(!truth)
  precision <- if (discoveries > 0L) true_positive / discoveries else NA_real_
  recall <- if (true_count > 0L) true_positive / true_count else NA_real_
  # v6: no discoveries contribute F1=0 for non-null targets.
  denominator <- 2 * true_positive + false_positive + false_negative
  f1 <- if (true_count > 0L) 2 * true_positive / denominator else NA_real_
  union_count <- true_positive + false_positive + false_negative

  data.frame(
    dgp = as.integer(dgp),
    replication_id = as.integer(replication_id),
    method = method,
    discoveries = discoveries,
    any_rejection = discoveries > 0L,
    true_positive = true_positive,
    false_positive = false_positive,
    false_negative = false_negative,
    recall = recall,
    background_fpr = if (background_count > 0L) false_positive / background_count else NA_real_,
    precision = precision,
    no_discovery = discoveries == 0L,
    f1 = f1,
    jaccard = if (union_count > 0L) true_positive / union_count else NA_real_,
    rank_auc = safe_rank_auc(score, truth),
    stringsAsFactors = FALSE
  )
}

run_mlsas_method <- function(sf_obj, dgp, replication_id, method, variables,
                             method_code, proximity_type, spec) {
  data <- prepare_method_data(sf_obj, variables)
  rf_seed <- as.integer(3000000L + 10000L * dgp + 10L * replication_id + method_code)
  permutation_seed <- as.integer(4000000L + 10000L * dgp + 10L * replication_id + method_code)
  set.seed(rf_seed)
  similarity <- if (identical(proximity_type, "path")) {
    compute_path_proximity_optimized(
      data,
      n_tree = as.integer(spec$inference$rf_ntree),
      workers = 1L,
      validate = FALSE,
      mtry = as.integer(spec$inference$rf_mtry),
      nodesize = as.integer(spec$inference$rf_nodesize)
    )
  } else {
    compute_terminal_proximity(
      data,
      n_tree = as.integer(spec$inference$rf_ntree),
      mtry = as.integer(spec$inference$rf_mtry),
      nodesize = as.integer(spec$inference$rf_nodesize)
    )
  }
  nb <- spdep::poly2nb(sf_obj, queen = TRUE, snap = sqrt(.Machine$double.eps))
  set.seed(permutation_seed)
  raw <- test.local.mlsas(
    id = sf_obj$id,
    nb = nb,
    S = similarity,
    style = spec$spatial_design$mlsas_weight_style,
    method = spec$inference$mlsas_method,
    alternative = spec$inference$mlsas_alternative,
    diag.zero = isTRUE(spec$spatial_design$diag_zero),
    p.adj = spec$inference$p_adjust,
    nperm = as.integer(spec$inference$permutations)
  )
  local <- data.frame(
    id = raw$id,
    statistic = raw$local.Moran,
    permutation_mean = raw$E,
    p_upper = raw$p_perm,
    p_adjusted = raw$p_perm_adj,
    detected = raw$p_perm_adj <= spec$inference$q & raw$local.Moran > raw$E,
    stringsAsFactors = FALSE
  )
  metrics <- make_method_metrics(
    # Larger local ML-SAS values indicate stronger positive local similarity.
    local, sf_obj$true_cluster, dgp, replication_id, method,
    score = orient_rank_auc_score(local$statistic, "higher")
  )
  list(local = local, metrics = metrics)
}

run_local_geary_method <- function(sf_obj, dgp, replication_id, method, variables,
                                   method_code, spec) {
  permutation_seed <- as.integer(4000000L + 10000L * dgp + 10L * replication_id + method_code)
  weights <- rgeoda::queen_weights(sf_obj)
  raw <- rgeoda::local_multigeary(
    weights,
    sf::st_drop_geometry(sf_obj[, variables, drop = FALSE]),
    permutations = as.integer(spec$inference$permutations),
    significance_cutoff = spec$inference$q,
    cpu_threads = 1L,
    seed = permutation_seed
  )
  p_value <- as.numeric(raw[["p_vals"]])
  statistic <- as.numeric(raw[["lisa_vals"]])
  labels <- as.character(raw[["labels"]][as.integer(raw[["c_vals"]]) + 1L])
  p_adjusted <- stats::p.adjust(p_value, method = spec$inference$p_adjust)
  local <- data.frame(
    id = sf_obj$id,
    statistic = statistic,
    p_value = p_value,
    p_adjusted = p_adjusted,
    class = labels,
    detected = labels == "Positive" & p_adjusted <= spec$inference$q,
    stringsAsFactors = FALSE
  )
  metrics <- make_method_metrics(
    # Smaller Local Geary values indicate stronger local similarity. Reverse
    # the score so all rank-AUC inputs use a higher-is-more-target-like scale.
    local, sf_obj$true_cluster, dgp, replication_id, method,
    score = orient_rank_auc_score(local$statistic, "lower")
  )
  list(local = local, metrics = metrics)
}

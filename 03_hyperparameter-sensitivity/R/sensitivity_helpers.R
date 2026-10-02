suppressPackageStartupMessages({
  library(dplyr)
  library(future)
  library(future.apply)
  library(randomForest)
  library(spdep)
  library(sf)
  library(tibble)
})

source("02_revision-simulation/R/dgp.R")
source("02_revision-simulation/R/mlsas_helpers.R")

make_addcl1_data <- function(data) {
  synthetic <- as.data.frame(lapply(data, function(col) sample(col, replace = TRUE)))
  data_all <- rbind(data, synthetic)
  y_all <- factor(c(rep("real", nrow(data)), rep("synthetic", nrow(synthetic))))
  list(data_all = data_all, y_all = y_all, n = nrow(data))
}

extract_tree_paths <- function(rf, data_all, tree_num) {
  tree <- randomForest::getTree(rf, k = tree_num, labelVar = FALSE)
  paths <- vector("list", nrow(data_all))

  for (i in seq_len(nrow(data_all))) {
    path <- integer()
    node <- 1L
    repeat {
      path <- c(path, node)
      split_var <- tree[node, "split var"]
      if (split_var == 0) {
        break
      }

      split_val <- tree[node, "split point"]
      split_var_name <- colnames(data_all)[split_var]
      value <- data_all[i, split_var_name]

      split_column <- data_all[[split_var_name]]

      if (is.factor(split_column) && !is.ordered(split_column)) {
        split_bits <- intToBits(split_val)
        value_level_index <- as.integer(value)
        if (split_bits[value_level_index] == as.raw(1)) {
          node <- tree[node, "left daughter"]
        } else {
          node <- tree[node, "right daughter"]
        }
      } else if ((if (is.ordered(split_column)) as.integer(value) else value) <= split_val) {
        node <- tree[node, "left daughter"]
      } else {
        node <- tree[node, "right daughter"]
      }
    }
    paths[[i]] <- path
  }

  paths
}

compute_path_proximity_hp <- function(data, n_tree, mtry, nodesize, workers = 1) {
  addcl <- make_addcl1_data(data)
  rf <- randomForest::randomForest(
    x = addcl$data_all,
    y = addcl$y_all,
    ntree = n_tree,
    mtry = mtry,
    nodesize = nodesize
  )

  tree_apply <- if (workers > 1) {
    function(x, fun, ...) future.apply::future_lapply(x, fun, ..., future.packages = "randomForest")
  } else {
    function(x, fun, ...) lapply(x, fun, ...)
  }

  all_paths <- tree_apply(
    seq_len(n_tree),
    function(tree_num) extract_tree_paths(rf, addcl$data_all, tree_num)
  )

  tree_jaccard <- function(tree_num) {
    paths <- all_paths[[tree_num]]
    mtx <- matrix(0, nrow = addcl$n, ncol = addcl$n)
    for (i in seq_len(addcl$n)) {
      Ti <- paths[[i]]
      for (j in i:addcl$n) {
        Tj <- paths[[j]]
        val <- length(intersect(Ti, Tj)) / length(union(Ti, Tj))
        mtx[i, j] <- val
        if (i != j) {
          mtx[j, i] <- val
        }
      }
    }
    mtx
  }

  mtx_list <- tree_apply(seq_len(n_tree), tree_jaccard)
  proximity <- Reduce(`+`, mtx_list) / n_tree
  rownames(proximity) <- rownames(data)
  colnames(proximity) <- rownames(data)
  proximity
}

compute_terminal_proximity_hp <- function(data, n_tree, mtry, nodesize) {
  addcl <- make_addcl1_data(data)
  rf <- randomForest::randomForest(
    x = addcl$data_all,
    y = addcl$y_all,
    ntree = n_tree,
    mtry = mtry,
    nodesize = nodesize
  )

  nodes <- attr(stats::predict(rf, addcl$data_all, nodes = TRUE), "nodes")
  nodes <- nodes[seq_len(addcl$n), , drop = FALSE]
  proximity <- matrix(0, nrow = addcl$n, ncol = addcl$n)

  for (tree_num in seq_len(n_tree)) {
    same_terminal <- outer(nodes[, tree_num], nodes[, tree_num], `==`)
    proximity <- proximity + same_terminal
  }

  proximity <- proximity / n_tree
  rownames(proximity) <- rownames(data)
  colnames(proximity) <- rownames(data)
  proximity
}

compute_urf_proximity_hp <- function(data, n_tree, mtry, nodesize, proximity_type, workers = 1) {
  proximity_type <- match.arg(proximity_type, c("path", "terminal"))
  if (proximity_type == "path") {
    return(compute_path_proximity_hp(data, n_tree, mtry, nodesize, workers = workers))
  }
  compute_terminal_proximity_hp(data, n_tree, mtry, nodesize)
}

upper_values <- function(mtx) {
  mtx[upper.tri(mtx)]
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3 || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
    return(NA_real_)
  }
  stats::cor(x[ok], y[ok], method = method)
}

rmse <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  sqrt(mean((x[ok] - y[ok])^2))
}

top_ids <- function(x, proportion) {
  n_top <- max(1L, ceiling(length(x) * proportion))
  order(x, decreasing = TRUE)[seq_len(n_top)]
}

jaccard_index <- function(a, b) {
  a <- unique(a)
  b <- unique(b)
  denom <- length(union(a, b))
  if (denom == 0) {
    return(NA_real_)
  }
  length(intersect(a, b)) / denom
}

summarise_true_cluster_recovery <- function(local_stat, true_cluster) {
  true_ids <- which(as.logical(true_cluster))
  false_ids <- which(!as.logical(true_cluster))
  ranks <- rank(-local_stat, ties.method = "average")

  tibble::tibble(
    top10_true_jaccard = jaccard_index(top_ids(local_stat, 0.10), true_ids),
    top15_true_jaccard = jaccard_index(top_ids(local_stat, 0.15), true_ids),
    true_local_mean = mean(local_stat[true_ids], na.rm = TRUE),
    false_local_mean = mean(local_stat[false_ids], na.rm = TRUE),
    true_false_local_gap = true_local_mean - false_local_mean,
    true_rank_q25 = stats::quantile(ranks[true_ids], 0.25, na.rm = TRUE),
    true_rank_median = stats::quantile(ranks[true_ids], 0.50, na.rm = TRUE),
    true_rank_q75 = stats::quantile(ranks[true_ids], 0.75, na.rm = TRUE)
  )
}

make_condition_label <- function(proximity_type, ntree, mtry, nodesize) {
  paste0(proximity_type, "_ntree", ntree, "_mtry", mtry, "_nodesize", nodesize)
}

comparison_group <- function(proximity_type_a, proximity_type_b, ntree_a, ntree_b, mtry_a, mtry_b, nodesize_a, nodesize_b) {
  changed <- c(
    proximity_type = proximity_type_a != proximity_type_b,
    ntree = ntree_a != ntree_b,
    mtry = mtry_a != mtry_b,
    nodesize = nodesize_a != nodesize_b
  )
  if (sum(changed) == 1 && changed[["proximity_type"]]) {
    return("proximity만 다름")
  }
  if (sum(changed) == 1 && changed[["ntree"]]) {
    return("ntree만 다름")
  }
  if (sum(changed) == 1 && changed[["mtry"]]) {
    return("mtry만 다름")
  }
  if (sum(changed) == 1 && changed[["nodesize"]]) {
    return("nodesize만 다름")
  }
  "여러 조건 다름"
}

pairwise_local_stability <- function(run_outputs) {
  meta <- dplyr::bind_rows(lapply(seq_along(run_outputs), function(i) {
    dplyr::mutate(run_outputs[[i]]$metrics, run_id = i, .before = 1)
  }))
  local_stats <- lapply(run_outputs, `[[`, "local_stat")
  local_stats <- do.call(cbind, local_stats)
  colnames(local_stats) <- meta$run_id

  dplyr::bind_rows(lapply(split(meta, meta$rep), function(rep_meta) {
    if (nrow(rep_meta) < 2) {
      return(tibble::tibble())
    }
    pairs <- utils::combn(seq_len(nrow(rep_meta)), 2)
    dplyr::bind_rows(lapply(seq_len(ncol(pairs)), function(pair_id) {
      a <- rep_meta[pairs[1, pair_id], ]
      b <- rep_meta[pairs[2, pair_id], ]
      stat_a <- local_stats[, as.character(a$run_id)]
      stat_b <- local_stats[, as.character(b$run_id)]

      tibble::tibble(
        rep = a$rep,
        condition_a = make_condition_label(a$proximity_type, a$ntree, a$mtry, a$nodesize),
        condition_b = make_condition_label(b$proximity_type, b$ntree, b$mtry, b$nodesize),
        proximity_type_a = a$proximity_type,
        proximity_type_b = b$proximity_type,
        ntree_a = a$ntree,
        ntree_b = b$ntree,
        mtry_a = a$mtry,
        mtry_b = b$mtry,
        nodesize_a = a$nodesize,
        nodesize_b = b$nodesize,
        comparison_group = comparison_group(
          a$proximity_type, b$proximity_type,
          a$ntree, b$ntree,
          a$mtry, b$mtry,
          a$nodesize, b$nodesize
        ),
        same_proximity_type = a$proximity_type == b$proximity_type,
        ntree_ratio = max(a$ntree, b$ntree) / min(a$ntree, b$ntree),
        mtry_distance = abs(a$mtry - b$mtry),
        nodesize_ratio = max(a$nodesize, b$nodesize) / min(a$nodesize, b$nodesize),
        local_pearson = safe_cor(stat_a, stat_b, "pearson"),
        local_spearman = safe_cor(stat_a, stat_b, "spearman"),
        local_rmse = rmse(stat_a, stat_b),
        top10_pairwise_jaccard = jaccard_index(top_ids(stat_a, 0.10), top_ids(stat_b, 0.10)),
        top15_pairwise_jaccard = jaccard_index(top_ids(stat_a, 0.15), top_ids(stat_b, 0.15))
      )
    }))
  }))
}

load_mixed_dgp43_inputs <- function(dgp_id = "dgp4.3", variables = c("x1", "x2", "cat1")) {
  dgp_obj <- load_dgp(dgp_id)
  list(
    dgp = dgp_obj,
    data = prepare_attribute_data(dgp_obj$sf, variables),
    nb = make_nb(dgp_obj$sf),
    true_cluster = sf::st_drop_geometry(dgp_obj$sf)$true_cluster
  )
}

compute_local_stat <- function(nb, S, ids) {
  local.mlsas(id = ids, nb = nb, S = S, style = "W", diag.zero = TRUE)$local.mlsas
}

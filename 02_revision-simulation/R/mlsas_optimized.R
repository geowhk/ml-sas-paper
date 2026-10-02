# Mechanism-preserving path-proximity engine for ML-SAS.
#
# This file intentionally lives beside, rather than replacing,
# `mlsas_helpers.R`.  The reference implementation remains the authority for
# regression tests.  The optimized engine keeps the same Addcl1 synthetic
# sample, randomForest fit, full root-to-terminal paths, per-tree Jaccard
# similarity, and arithmetic mean across trees.

.extract_real_tree_paths_optimized <- function(tree, data_all, n_real) {
  paths <- vector("list", n_real)

  for (i in seq_len(n_real)) {
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

.path_jaccard_sparse_optimized <- function(paths, n_nodes) {
  n <- length(paths)
  path_lengths <- lengths(paths)

  incidence <- Matrix::sparseMatrix(
    i = rep.int(seq_len(n), path_lengths),
    j = unlist(paths, use.names = FALSE),
    x = 1,
    dims = c(n, n_nodes)
  )

  intersection_sizes <- as.matrix(Matrix::tcrossprod(incidence))
  union_sizes <- outer(path_lengths, path_lengths, "+") - intersection_sizes
  intersection_sizes / union_sizes
}

.process_tree_chunk_optimized <- function(
  tree_ids,
  rf,
  data_all,
  n_real,
  predicted_nodes = NULL
) {
  chunk_sum <- matrix(0, nrow = n_real, ncol = n_real)
  terminal_mismatches <- 0L

  for (tree_num in tree_ids) {
    tree <- randomForest::getTree(rf, k = tree_num, labelVar = FALSE)
    paths <- .extract_real_tree_paths_optimized(tree, data_all, n_real)

    if (!is.null(predicted_nodes)) {
      extracted_terminal_nodes <- vapply(paths, function(path) path[length(path)], numeric(1))
      terminal_mismatches <- terminal_mismatches + sum(
        extracted_terminal_nodes != predicted_nodes[seq_len(n_real), tree_num]
      )
    }

    chunk_sum <- chunk_sum + .path_jaccard_sparse_optimized(paths, nrow(tree))
  }

  list(sum = chunk_sum, terminal_mismatches = terminal_mismatches)
}

.chunk_tree_ids_optimized <- function(n_tree, chunk_size) {
  tree_ids <- seq_len(n_tree)
  split(tree_ids, ceiling(tree_ids / chunk_size))
}

.apply_chunks_optimized <- function(chunks, fun, workers) {
  if (workers <= 1L || length(chunks) == 1L) {
    return(lapply(chunks, fun))
  }

  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  strategy <- if (future::supportsMulticore()) {
    future::multicore
  } else {
    future::multisession
  }
  future::plan(strategy, workers = workers)

  future.apply::future_lapply(
    chunks,
    fun,
    future.packages = c("randomForest", "Matrix"),
    future.seed = TRUE
  )
}

compute_path_proximity_optimized <- function(
  data,
  n_tree = 500,
  workers = 1,
  validate = FALSE,
  chunk_size = NULL,
  mtry = NULL,
  nodesize = NULL,
  synthetic_data = NULL
) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  n <- nrow(data)
  n_tree <- as.integer(n_tree)
  workers <- as.integer(workers)
  if (is.null(chunk_size)) {
    chunk_size <- max(1L, ceiling(n_tree / max(1L, workers * 2L)))
  } else {
    chunk_size <- as.integer(chunk_size)
  }

  if (n < 2L) {
    stop("`data` must contain at least two observations.", call. = FALSE)
  }
  if (n_tree < 1L || workers < 1L || chunk_size < 1L) {
    stop("`n_tree`, `workers`, and `chunk_size` must be positive integers.", call. = FALSE)
  }

  # The default keeps the Addcl1 construction byte-for-byte equivalent to the
  # reference. An externally generated Addcl1 sample can be supplied to make
  # paired specification/ablation comparisons use identical shared columns.
  if (is.null(synthetic_data)) {
    synthetic <- as.data.frame(lapply(data, function(col) sample(col, replace = TRUE)))
  } else {
    if (!is.data.frame(synthetic_data) ||
        nrow(synthetic_data) != n ||
        !identical(names(synthetic_data), names(data))) {
      stop("`synthetic_data` must match `data` in rows and column names.", call. = FALSE)
    }
    synthetic <- synthetic_data
  }
  data_all <- rbind(data, synthetic)
  y_all <- factor(c(rep("real", n), rep("synthetic", n)))

  rf_args <- list(x = data_all, y = y_all, ntree = n_tree)
  if (!is.null(mtry)) {
    rf_args$mtry <- mtry
  }
  if (!is.null(nodesize)) {
    rf_args$nodesize <- nodesize
  }
  rf <- do.call(randomForest::randomForest, rf_args)

  predicted_nodes <- NULL
  if (isTRUE(validate)) {
    predicted_nodes <- attr(
      stats::predict(rf, data_all, nodes = TRUE),
      "nodes"
    )
  }

  chunks <- .chunk_tree_ids_optimized(n_tree, chunk_size)
  process_chunk <- function(tree_ids) {
    .process_tree_chunk_optimized(
      tree_ids = tree_ids,
      rf = rf,
      data_all = data_all,
      n_real = n,
      predicted_nodes = predicted_nodes
    )
  }

  proximity_sum <- matrix(0, nrow = n, ncol = n)
  terminal_mismatches <- 0L

  if (workers <= 1L) {
    # Stream chunk sums into a single accumulator instead of retaining them.
    for (tree_ids in chunks) {
      result <- process_chunk(tree_ids)
      proximity_sum <- proximity_sum + result$sum
      terminal_mismatches <- terminal_mismatches + result$terminal_mismatches
    }
  } else {
    chunk_results <- .apply_chunks_optimized(chunks, process_chunk, workers = workers)
    terminal_mismatches <- sum(vapply(
      chunk_results,
      `[[`,
      integer(1),
      "terminal_mismatches"
    ))
    for (result in chunk_results) {
      proximity_sum <- proximity_sum + result$sum
    }
  }

  if (isTRUE(validate)) {
    if (terminal_mismatches > 0L) {
      stop(
        sprintf("Path validation failed: %d terminal-node mismatches.", terminal_mismatches),
        call. = FALSE
      )
    }
    message("Validation Success: optimized paths match randomForest terminal nodes.")
  }

  proximity <- proximity_sum / n_tree
  rownames(proximity) <- rownames(data)
  colnames(proximity) <- rownames(data)
  proximity
}

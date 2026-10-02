mlsas_components_optimized <- function(nb, similarity, style = "W") {
  stopifnot(is.matrix(similarity), nrow(similarity) == ncol(similarity))
  listw <- spdep::nb2listw(nb, style = style, zero.policy = FALSE)
  weights <- as.matrix(as_dgRMatrix_listw(listw))
  diag(weights) <- 0
  similarity <- as.matrix(similarity)
  diag(similarity) <- 0

  if (any(!is.finite(similarity)) || any(similarity < 0)) {
    stop("The similarity matrix must contain finite, non-negative values.")
  }
  if (sum(weights) <= 0 || sum(similarity) <= 0) {
    stop("Weights and off-diagonal similarities must have positive sums.")
  }

  n <- nrow(similarity)
  P <- n * (n - 1) * weights / sum(weights)
  Q <- similarity / sum(similarity)
  list(n = n, P = P, Q = Q)
}

mlsas_observed_optimized <- function(nb, similarity, style = "W") {
  parts <- mlsas_components_optimized(nb, similarity, style)
  products <- parts$P * parts$Q
  list(
    global = sum(products),
    local = parts$n * rowSums(products),
    P = parts$P,
    Q = parts$Q
  )
}

make_permutation_indices <- function(n, nperm, seed) {
  set.seed(seed)
  vapply(seq_len(nperm), function(i) sample.int(n), integer(n))
}

mlsas_permutation_inference_optimized <- function(
  nb,
  similarity,
  permutations,
  workers = 1L,
  style = "W",
  alpha = 0.05,
  p_adjust_method = "BH"
) {
  observed <- mlsas_observed_optimized(nb, similarity, style)
  n <- length(observed$local)
  stopifnot(nrow(permutations) == n)
  nperm <- ncol(permutations)

  chunk_ids <- split(
    seq_len(nperm),
    cut(seq_len(nperm), breaks = min(as.integer(workers), nperm), labels = FALSE)
  )
  compute_chunk <- function(ids) {
    local <- matrix(NA_real_, nrow = n, ncol = length(ids))
    global <- numeric(length(ids))
    for (j in seq_along(ids)) {
      pi <- permutations[, ids[j]]
      qp <- observed$Q[pi, pi, drop = FALSE]
      products <- observed$P * qp
      global[j] <- sum(products)
      local[, j] <- n * rowSums(products)
    }
    list(global = global, local = local)
  }

  if (workers > 1L && length(chunk_ids) > 1L) {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    strategy <- if (future::supportsMulticore()) future::multicore else future::multisession
    future::plan(strategy, workers = min(as.integer(workers), length(chunk_ids)))
    chunks <- future.apply::future_lapply(chunk_ids, compute_chunk, future.seed = TRUE)
  } else {
    chunks <- lapply(chunk_ids, compute_chunk)
  }

  global_perm <- unlist(lapply(chunks, `[[`, "global"), use.names = FALSE)
  local_perm <- do.call(cbind, lapply(chunks, `[[`, "local"))
  global_high <- (sum(global_perm >= observed$global) + 1) / (nperm + 1)
  global_low <- (sum(global_perm <= observed$global) + 1) / (nperm + 1)
  global_two <- min(1, 2 * min(global_high, global_low))
  local_high <- (rowSums(local_perm >= observed$local) + 1) / (nperm + 1)
  local_low <- (rowSums(local_perm <= observed$local) + 1) / (nperm + 1)
  local_two <- pmin(1, 2 * pmin(local_high, local_low))
  local_adjusted <- stats::p.adjust(local_two, method = p_adjust_method)
  local_mean <- rowMeans(local_perm)
  local_sd <- apply(local_perm, 1, stats::sd)
  direction <- ifelse(
    observed$local >= local_mean,
    "Similarity cluster",
    "Dissimilarity/outlier"
  )
  classification <- ifelse(local_adjusted <= alpha, direction, "Not significant")

  list(
    global = data.frame(
      observed = observed$global,
      permutation_mean = mean(global_perm),
      permutation_sd = stats::sd(global_perm),
      p_high = global_high,
      p_low = global_low,
      p_two_sided = global_two,
      nperm = nperm
    ),
    local = data.frame(
      local_mlsas = observed$local,
      permutation_mean = local_mean,
      permutation_sd = local_sd,
      p_high = local_high,
      p_low = local_low,
      p_two_sided = local_two,
      p_adjusted = local_adjusted,
      classification = classification
    )
  )
}

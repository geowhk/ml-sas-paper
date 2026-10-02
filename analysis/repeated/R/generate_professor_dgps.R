# Section 4 repeated-DGP generator based on the professor's original script.
# Source supplied by the user: source_snapshot/DGP_generation.R
# This file keeps the original data-generating defaults and RNG draw order.

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(spdep)
})

make_professor_grid_data <- function(n_row = 20,
                                     n_col = 20,
                                     cell_type = c("square", "hex"),
                                     queen = TRUE,
                                     cell_width = 1) {
  cell_type <- match.arg(cell_type)

  if (cell_type == "square") {
    bbox <- sf::st_bbox(c(xmin = 0, ymin = 0, xmax = n_col, ymax = n_row))
    grd <- sf::st_make_grid(
      sf::st_as_sfc(bbox),
      n = c(n_col, n_row),
      square = TRUE
    )
    grd <- sf::st_as_sf(grd) |>
      dplyr::mutate(id = dplyr::row_number())
    cent <- sf::st_coordinates(sf::st_centroid(grd))
    grd <- grd |>
      dplyr::mutate(cx = cent[, 1], cy = cent[, 2])
  }

  if (cell_type == "hex") {
    w <- cell_width
    h <- 2 * w / sqrt(3)
    y_step <- 3 * h / 4
    centers <- expand.grid(row = seq_len(n_row), col = seq_len(n_col)) |>
      dplyr::arrange(row, col) |>
      dplyr::mutate(
        cx = (col - 1) * w + w / 2 + dplyr::if_else(row %% 2 == 0, w / 2, 0),
        cy = (n_row - row) * y_step + h / 2
      )

    make_hex <- function(cx, cy, w, h) {
      coords <- matrix(
        c(
          cx,         cy + h / 2,
          cx + w / 2, cy + h / 4,
          cx + w / 2, cy - h / 4,
          cx,         cy - h / 2,
          cx - w / 2, cy - h / 4,
          cx - w / 2, cy + h / 4,
          cx,         cy + h / 2
        ),
        ncol = 2,
        byrow = TRUE
      )
      sf::st_polygon(list(coords))
    }

    geom <- sf::st_sfc(Map(function(x, y) make_hex(x, y, w, h), centers$cx, centers$cy))
    grd <- sf::st_sf(
      id = seq_len(nrow(centers)),
      row = centers$row,
      col = centers$col,
      cx = centers$cx,
      cy = centers$cy,
      geometry = geom
    )
  }

  nb <- spdep::poly2nb(grd, queen = queen, snap = sqrt(.Machine$double.eps))
  lw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  list(sf = grd, nb = nb, lw = lw)
}

make_professor_true_cluster <- function(sf_data,
                                        center_col = 10.5,
                                        center_row = 10.5,
                                        radius = 4.5) {
  sf_data |>
    dplyr::mutate(
      dist_to_center = sqrt((col - center_col)^2 + (row - center_row)^2),
      true_cluster = dist_to_center <= radius,
      cluster_type = dplyr::if_else(true_cluster, "cluster", "background")
    )
}

simulate_professor_dgp1 <- function(grid_obj,
                                    seed = 123,
                                    radius = 3.2,
                                    high_mean = 2.5,
                                    low_mean = -2.5,
                                    mid_mean = 0,
                                    cluster_sd = 0.15,
                                    background_mean = 0,
                                    background_sd = 1.25) {
  set.seed(seed)
  sf_data <- grid_obj$sf
  n <- nrow(sf_data)
  hh_center <- c(col = 5, row = 5)
  ll_center <- c(col = 16, row = 16)
  mm_center <- c(col = 10.5, row = 10.5)
  sf_data <- sf_data |>
    dplyr::mutate(
      d_hh = sqrt((col - hh_center["col"])^2 + (row - hh_center["row"])^2),
      d_ll = sqrt((col - ll_center["col"])^2 + (row - ll_center["row"])^2),
      d_mm = sqrt((col - mm_center["col"])^2 + (row - mm_center["row"])^2),
      cluster_type = dplyr::case_when(
        d_hh <= radius ~ "HH",
        d_ll <= radius ~ "LL",
        d_mm <= radius ~ "MM",
        TRUE ~ "background"
      ),
      true_cluster = cluster_type != "background"
    )
  x <- stats::rnorm(n, mean = background_mean, sd = background_sd)
  idx_hh <- sf_data$cluster_type == "HH"
  idx_ll <- sf_data$cluster_type == "LL"
  idx_mm <- sf_data$cluster_type == "MM"
  x[idx_hh] <- stats::rnorm(sum(idx_hh), mean = high_mean, sd = cluster_sd)
  x[idx_ll] <- stats::rnorm(sum(idx_ll), mean = low_mean, sd = cluster_sd)
  x[idx_mm] <- stats::rnorm(sum(idx_mm), mean = mid_mean, sd = cluster_sd)
  sf_data |>
    dplyr::mutate(
      x = x,
      cluster_type = factor(cluster_type, levels = c("background", "HH", "LL", "MM")),
      true_hh = cluster_type == "HH",
      true_ll = cluster_type == "LL",
      true_mm = cluster_type == "MM",
      dgp = "DGP1_univariate_geary_type_clusters"
    ) |>
    dplyr::select(
      id, row, col, cx, cy, x, cluster_type, true_cluster,
      true_hh, true_ll, true_mm, dgp, geometry
    )
}

simulate_professor_dgp2 <- function(grid_obj, seed = 123, n_cat = 3) {
  set.seed(seed)
  sf_data <- grid_obj$sf
  n <- nrow(sf_data)
  sf_data |>
    dplyr::mutate(
      x1 = stats::rnorm(n, 0, 1),
      x2 = stats::rnorm(n, 0, 1),
      x3 = stats::rnorm(n, 0, 1),
      cat1 = factor(sample(paste0("C", seq_len(n_cat)), size = n, replace = TRUE)),
      cluster_type = factor("background"),
      true_cluster = FALSE,
      dgp = "DGP2_null_mixed_CSR"
    )
}

simulate_professor_dgp3 <- function(grid_obj,
                                    seed = 123,
                                    radius = 4.5,
                                    effect_size = 1.5,
                                    cluster_sd = 0.35,
                                    background_sd = 1.0) {
  set.seed(seed)
  sf_data <- make_professor_true_cluster(grid_obj$sf, radius = radius)
  n <- nrow(sf_data)
  x1 <- stats::rnorm(n, 0, background_sd)
  x2 <- stats::rnorm(n, 0, background_sd)
  x3 <- stats::rnorm(n, 0, background_sd)
  idx_cluster <- sf_data$true_cluster
  x1[idx_cluster] <- stats::rnorm(sum(idx_cluster), mean = effect_size, sd = cluster_sd)
  x2[idx_cluster] <- stats::rnorm(sum(idx_cluster), mean = effect_size, sd = cluster_sd)
  sf_data |>
    dplyr::mutate(
      x1 = x1,
      x2 = x2,
      x3 = x3,
      cluster_type = factor(cluster_type, levels = c("background", "cluster")),
      dgp = "DGP3_linear_multivariate_continuous"
    ) |>
    dplyr::select(-dist_to_center)
}

simulate_professor_dgp4 <- function(grid_obj,
                                    seed = 123,
                                    radius = 4.5,
                                    center_col = NULL,
                                    center_row = NULL,
                                    cluster_center = 0.70,
                                    cluster_latent_sd = 0.25,
                                    cluster_noise_sd = 0.12,
                                    background_sd = 1.20,
                                    nonlinear_strength = 0.35,
                                    cat_prob_cluster = c(A = 0.90, B = 0.05, C = 0.05),
                                    cat_prob_background = c(A = 0.20, B = 0.40, C = 0.40),
                                    include_x3 = FALSE) {
  set.seed(seed)
  sf_data <- grid_obj$sf
  n <- nrow(sf_data)
  if (is.null(center_col)) center_col <- mean(range(sf_data$col, na.rm = TRUE))
  if (is.null(center_row)) center_row <- mean(range(sf_data$row, na.rm = TRUE))
  sf_data <- sf_data |>
    dplyr::mutate(
      dist_to_center = sqrt((col - center_col)^2 + (row - center_row)^2),
      true_cluster = dist_to_center <= radius,
      cluster_type = dplyr::if_else(true_cluster, "cluster", "background")
    )
  idx_cluster <- sf_data$true_cluster
  idx_background <- !sf_data$true_cluster
  x1 <- stats::rnorm(n, mean = 0, sd = background_sd)
  x2 <- stats::rnorm(n, mean = 0, sd = background_sd)
  u <- stats::rnorm(sum(idx_cluster), mean = 0, sd = cluster_latent_sd)
  x1[idx_cluster] <- cluster_center + u +
    stats::rnorm(sum(idx_cluster), mean = 0, sd = cluster_noise_sd)
  x2[idx_cluster] <- cluster_center +
    nonlinear_strength * (u^2 - cluster_latent_sd^2) +
    stats::rnorm(sum(idx_cluster), mean = 0, sd = cluster_noise_sd)
  x1 <- as.numeric(scale(x1))
  x2 <- as.numeric(scale(x2))
  x1x2 <- x1 * x2
  interaction_sign <- ifelse(x1x2 > 0, "same_sign", "opposite_sign")
  if (include_x3) x3 <- stats::rnorm(n, mean = 0, sd = 1)

  # Preserve the professor's RNG order: cluster categories first, background second.
  cat1 <- character(n)
  cat1[idx_cluster] <- sample(
    names(cat_prob_cluster),
    size = sum(idx_cluster),
    replace = TRUE,
    prob = cat_prob_cluster
  )
  cat1[idx_background] <- sample(
    names(cat_prob_background),
    size = sum(idx_background),
    replace = TRUE,
    prob = cat_prob_background
  )

  out <- sf_data |>
    dplyr::mutate(
      x1 = x1,
      x2 = x2,
      x1x2 = x1x2,
      interaction_sign = factor(interaction_sign, levels = c("opposite_sign", "same_sign")),
      cat1 = factor(cat1, levels = c("A", "B", "C")),
      cluster_type = factor(cluster_type, levels = c("background", "cluster")),
      dgp = "DGP4_nonlinear_mixed_interaction"
    )
  if (include_x3) {
    out |>
      dplyr::mutate(x3 = x3) |>
      dplyr::select(
        id, row, col, cx, cy, x1, x2, x3, x1x2, interaction_sign,
        cat1, true_cluster, cluster_type, dgp, geometry
      )
  } else {
    out |>
      dplyr::select(
        id, row, col, cx, cy, x1, x2, x1x2, interaction_sign,
        cat1, true_cluster, cluster_type, dgp, geometry
      )
  }
}

professor_dataset_seed <- function(dgp, replication_id) {
  dgp <- as.integer(dgp)
  replication_id <- as.integer(replication_id)
  stopifnot(dgp %in% 1:4, replication_id >= 1L)
  as.integer(2000000L + dgp * 100000L + replication_id)
}

generate_professor_dgp <- function(dgp, replication_id) {
  dgp <- as.integer(dgp)
  seed <- professor_dataset_seed(dgp, replication_id)
  grid_obj <- make_professor_grid_data(
    n_row = 20,
    n_col = 20,
    cell_type = "hex",
    queen = TRUE,
    cell_width = 1
  )
  switch(
    as.character(dgp),
    `1` = simulate_professor_dgp1(grid_obj, seed = seed),
    `2` = simulate_professor_dgp2(grid_obj, seed = seed),
    `3` = simulate_professor_dgp3(grid_obj, seed = seed),
    `4` = simulate_professor_dgp4(grid_obj, seed = seed),
    stop("Unsupported DGP: ", dgp, call. = FALSE)
  )
}

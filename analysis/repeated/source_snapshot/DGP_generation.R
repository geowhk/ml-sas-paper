# Packages ----------------------------------------------------------------

library(tidyverse)
library(sf)
library(tmap)
library(spdep)
library(sfdep)
library(spatialreg)
library(Guerry)
library(moments)
options(scipen = 999)


make_grid_data <- function(n_row = 20,
                           n_col = 20,
                           cell_type = c("square", "hex"),
                           queen = TRUE,
                           cell_width = 1) {
  
  cell_type <- match.arg(cell_type)
  
  if (cell_type == "square") {
    
    bbox <- st_bbox(c(
      xmin = 0, ymin = 0,
      xmax = n_col, ymax = n_row
    ))
    
    grd <- st_make_grid(
      st_as_sfc(bbox),
      n = c(n_col, n_row),
      square = TRUE
    )
    
    grd <- st_as_sf(grd) |>
      mutate(id = row_number())
    
    cent <- st_coordinates(st_centroid(grd))
    
    grd <- grd |>
      mutate(
        cx = cent[, 1],
        cy = cent[, 2]
      )
  }
  
  if (cell_type == "hex") {
    
    # pointy-top hexagon
    # cell_width는 hexagon의 좌우 폭
    w <- cell_width
    h <- 2 * w / sqrt(3)
    
    # pointy-top hexagon에서 행 간 중심 간격
    y_step <- 3 * h / 4
    
    centers <- expand.grid(
      row = seq_len(n_row),
      col = seq_len(n_col)
    ) |>
      arrange(row, col) |>
      mutate(
        cx = (col - 1) * w + w / 2 +
          if_else(row %% 2 == 0, w / 2, 0),
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
      
      st_polygon(list(coords))
    }
    
    geom <- st_sfc(
      Map(
        function(x, y) make_hex(x, y, w, h),
        centers$cx,
        centers$cy
      )
    )
    
    grd <- st_sf(
      id = seq_len(nrow(centers)),
      row = centers$row,
      col = centers$col,
      cx = centers$cx,
      cy = centers$cy,
      geometry = geom
    )
  }
  
  nb <- poly2nb(grd, queen = queen, snap = sqrt(.Machine$double.eps))
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  list(
    sf = grd,
    nb = nb,
    lw = lw
  )
}

# DGP ---------------------------------------------------------------------

# ------------------------------------------------------------
# Helper: central true cluster
# ------------------------------------------------------------

make_true_cluster <- function(sf_data,
                              center_col = 10.5,
                              center_row = 10.5,
                              radius = 4.5) {
  
  sf_data |>
    mutate(
      dist_to_center = sqrt((col - center_col)^2 + 
                              (row - center_row)^2),
      true_cluster = dist_to_center <= radius,
      cluster_type = if_else(true_cluster, "cluster", "background")
    )
}

# DGP 1: Univariate continuous Geary-type clusters

# ------------------------------------------------------------
# DGP 1
# Univariate continuous Geary-type clusters
#
# Purpose:
# - HH: high-high cluster
# - LL: low-low cluster
# - MM: moderate-moderate cluster near the global mean
# - Moran and Getis-Ord are expected to detect mainly HH and LL.
# - Local Geary and local ML-SAS are expected to detect HH, LL, and MM.
# ------------------------------------------------------------

simulate_dgp1_univariate_geary_clusters <- function(grid_obj,
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
  
  # Cluster centers on a 20 x 20 hexagonal grid
  hh_center <- c(col = 5,    row = 5)
  ll_center <- c(col = 16,   row = 16)
  mm_center <- c(col = 10.5, row = 10.5)
  
  sf_data <- sf_data |>
    mutate(
      d_hh = sqrt((col - hh_center["col"])^2 +
                    (row - hh_center["row"])^2),
      d_ll = sqrt((col - ll_center["col"])^2 +
                    (row - ll_center["row"])^2),
      d_mm = sqrt((col - mm_center["col"])^2 +
                    (row - mm_center["row"])^2),
      
      cluster_type = case_when(
        d_hh <= radius ~ "HH",
        d_ll <= radius ~ "LL",
        d_mm <= radius ~ "MM",
        TRUE ~ "background"
      ),
      
      true_cluster = cluster_type != "background"
    )
  
  # Background: spatially random continuous values
  x <- rnorm(n, mean = background_mean, sd = background_sd)
  
  # Cluster indices
  idx_hh <- sf_data$cluster_type == "HH"
  idx_ll <- sf_data$cluster_type == "LL"
  idx_mm <- sf_data$cluster_type == "MM"
  
  # HH and LL: extreme homogeneous clusters
  x[idx_hh] <- rnorm(sum(idx_hh), mean = high_mean, sd = cluster_sd)
  x[idx_ll] <- rnorm(sum(idx_ll), mean = low_mean,  sd = cluster_sd)
  
  # MM: locally homogeneous but not globally extreme
  x[idx_mm] <- rnorm(sum(idx_mm), mean = mid_mean,  sd = cluster_sd)
  
  sf_data |>
    mutate(
      x = x,
      cluster_type = factor(
        cluster_type,
        levels = c("background", "HH", "LL", "MM")
      ),
      true_hh = cluster_type == "HH",
      true_ll = cluster_type == "LL",
      true_mm = cluster_type == "MM",
      dgp = "DGP1_univariate_geary_type_clusters"
    ) |>
    select(
      id, row, col, cx, cy,
      x,
      cluster_type,
      true_cluster,
      true_hh,
      true_ll,
      true_mm,
      dgp,
      geometry
    )
}

# DGP 2: Null mixed-type CSR

# ------------------------------------------------------------
# DGP 2
# Null mixed-type CSR
#
# Purpose:
# - No true spatial cluster.
# - Continuous and categorical variables are spatially random.
# - Used to evaluate false positive behavior.
# ------------------------------------------------------------

simulate_dgp2_null_mixed <- function(grid_obj,
                                     seed = 123,
                                     n_cat = 3) {
  
  set.seed(seed)
  
  sf_data <- grid_obj$sf
  n <- nrow(sf_data)
  
  sf_data |>
    mutate(
      x1 = rnorm(n, 0, 1),
      x2 = rnorm(n, 0, 1),
      x3 = rnorm(n, 0, 1),
      cat1 = factor(
        sample(
          paste0("C", seq_len(n_cat)),
          size = n,
          replace = TRUE
        )
      ),
      cluster_type = factor("background"),
      true_cluster = FALSE,
      dgp = "DGP2_null_mixed_CSR"
    )
}

# DGP 3: Linear multivariate continuous cluster

# ------------------------------------------------------------
# DGP 3
# Linear multivariate continuous cluster
#
# Purpose:
# - A simple multivariate continuous cluster.
# - x1 and x2 have mean shifts inside the true cluster.
# - x3 is a noise variable.
# - Used for comparison with Multivariate Local Geary.
# ------------------------------------------------------------

simulate_dgp3_linear_multivariate_continuous <- function(grid_obj,
                                                         seed = 123,
                                                         radius = 4.5,
                                                         effect_size = 1.5,
                                                         cluster_sd = 0.35,
                                                         background_sd = 1.0) {
  
  set.seed(seed)
  
  sf_data <- grid_obj$sf |>
    make_true_cluster(radius = radius)
  
  n <- nrow(sf_data)
  
  x1 <- rnorm(n, 0, background_sd)
  x2 <- rnorm(n, 0, background_sd)
  x3 <- rnorm(n, 0, background_sd)
  
  idx_cluster <- sf_data$true_cluster
  
  x1[idx_cluster] <- rnorm(
    sum(idx_cluster),
    mean = effect_size,
    sd = cluster_sd
  )
  
  x2[idx_cluster] <- rnorm(
    sum(idx_cluster),
    mean = effect_size,
    sd = cluster_sd
  )
  
  sf_data |>
    mutate(
      x1 = x1,
      x2 = x2,
      x3 = x3,
      cluster_type = factor(
        cluster_type,
        levels = c("background", "cluster")
      ),
      dgp = "DGP3_linear_multivariate_continuous"
    ) |>
    select(-dist_to_center)
}

# DGP 4: Nonlinear mixed-type interaction cluster

simulate_dgp4_nonlinear_mixed_interaction <- function(grid_obj,
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
  
  # ------------------------------------------------------------
  # Define the true cluster
  # ------------------------------------------------------------
  
  if (is.null(center_col)) {
    center_col <- mean(range(sf_data$col, na.rm = TRUE))
  }
  
  if (is.null(center_row)) {
    center_row <- mean(range(sf_data$row, na.rm = TRUE))
  }
  
  sf_data <- sf_data |>
    mutate(
      dist_to_center = sqrt((col - center_col)^2 +
                              (row - center_row)^2),
      true_cluster = dist_to_center <= radius,
      cluster_type = if_else(true_cluster, "cluster", "background")
    )
  
  idx_cluster <- sf_data$true_cluster
  idx_background <- !sf_data$true_cluster
  
  # ------------------------------------------------------------
  # Background values
  #
  # Outside the true cluster, x1 and x2 are spatially independent
  # random variables with relatively large dispersion.
  # ------------------------------------------------------------
  
  x1 <- rnorm(n, mean = 0, sd = background_sd)
  x2 <- rnorm(n, mean = 0, sd = background_sd)
  
  # ------------------------------------------------------------
  # Cluster values
  #
  # Inside the true cluster, x1 and x2 form a moderately homogeneous
  # joint structure. The signal is not intended to be a simple
  # high-high univariate cluster, but a joint nonlinear structure
  # that is easier to detect when the two variables are considered
  # together.
  #
  # u is a local latent variable with small variance.
  # x1 varies around a moderate value.
  # x2 varies around the same moderate level but also contains
  # a weak nonlinear component of u.
  # ------------------------------------------------------------
  
  u <- rnorm(sum(idx_cluster), mean = 0, sd = cluster_latent_sd)
  
  x1[idx_cluster] <- cluster_center +
    u +
    rnorm(sum(idx_cluster), mean = 0, sd = cluster_noise_sd)
  
  x2[idx_cluster] <- cluster_center +
    nonlinear_strength * (u^2 - cluster_latent_sd^2) +
    rnorm(sum(idx_cluster), mean = 0, sd = cluster_noise_sd)
  
  # ------------------------------------------------------------
  # Standardize marginal scales
  #
  # This prevents the result from being driven simply by different
  # raw measurement scales of x1 and x2.
  # ------------------------------------------------------------
  
  x1 <- as.numeric(scale(x1))
  x2 <- as.numeric(scale(x2))
  
  # ------------------------------------------------------------
  # Diagnostic nonlinear interaction
  #
  # x1x2 is generated only for checking whether the nonlinear
  # joint structure is present.
  #
  # Important:
  # x1x2 should NOT be used as an input variable in ML-SAS.
  # ------------------------------------------------------------
  
  x1x2 <- x1 * x2
  
  interaction_sign <- ifelse(
    x1x2 > 0,
    "same_sign",
    "opposite_sign"
  )
  
  # ------------------------------------------------------------
  # Optional CSR-like noise variable
  #
  # For the main DGP 4 analysis, I recommend include_x3 = FALSE.
  # x3 can be used later as an additional robustness check.
  # ------------------------------------------------------------
  
  if (include_x3) {
    x3 <- rnorm(n, mean = 0, sd = 1)
  }
  
  # ------------------------------------------------------------
  # Categorical variable
  #
  # Inside the true cluster, category A is much more likely but
  # not deterministic. Outside the cluster, categories are more
  # balanced and spatially random.
  #
  # This is the key mixed-type component:
  # it should strengthen full ML-SAS relative to continuous-only ML-SAS.
  # ------------------------------------------------------------
  
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
  
  # ------------------------------------------------------------
  # Return sf data
  # ------------------------------------------------------------
  
  out <- sf_data |>
    mutate(
      x1 = x1,
      x2 = x2,
      x1x2 = x1x2,
      interaction_sign = factor(
        interaction_sign,
        levels = c("opposite_sign", "same_sign")
      ),
      cat1 = factor(cat1, levels = c("A", "B", "C")),
      cluster_type = factor(
        cluster_type,
        levels = c("background", "cluster")
      ),
      dgp = "DGP4_nonlinear_mixed_interaction"
    )
  
  if (include_x3) {
    out <- out |>
      mutate(x3 = x3) |>
      select(
        id, row, col, cx, cy,
        x1, x2, x3,
        x1x2, interaction_sign,
        cat1,
        true_cluster,
        cluster_type,
        dgp,
        geometry
      )
  } else {
    out <- out |>
      select(
        id, row, col, cx, cy,
        x1, x2,
        x1x2, interaction_sign,
        cat1,
        true_cluster,
        cluster_type,
        dgp,
        geometry
      )
  }
  
  out
}

# ------------------------------------------------------------
# Generate all DGPs using one common grid
# ------------------------------------------------------------

simulate_all_dgps <- function(n_row = 20,
                              n_col = 20,
                              cell_type = c("hex", "square"),
                              queen = TRUE,
                              cell_width = 1,
                              seed = 2026) {
  
  cell_type <- match.arg(cell_type)
  
  # Create one common grid and spatial-weights object
  grid_obj <- make_grid_data(
    n_row = n_row,
    n_col = n_col,
    cell_type = cell_type,
    queen = queen,
    cell_width = cell_width
  )
  
  # Generate the four DGPs on the same spatial support
  dgp1 <- simulate_dgp1_univariate_geary_clusters(
    grid_obj = grid_obj,
    seed = seed
  )
  
  dgp2 <- simulate_dgp2_null_mixed(
    grid_obj = grid_obj,
    seed = seed + 1L
  )
  
  dgp3 <- simulate_dgp3_linear_multivariate_continuous(
    grid_obj = grid_obj,
    seed = seed + 2L
  )
  
  dgp4 <- simulate_dgp4_nonlinear_mixed_interaction(
    grid_obj = grid_obj,
    seed = seed + 3L
  )
  
  list(
    grid = grid_obj,
    dgp1 = dgp1,
    dgp2 = dgp2,
    dgp3 = dgp3,
    dgp4 = dgp4
  )
}

# ------------------------------------------------------------
# Generate simulation data
# ------------------------------------------------------------

sim_data <- simulate_all_dgps(
  n_row = 20,
  n_col = 20,
  cell_type = "hex",
  queen = TRUE,
  cell_width = 1,
  seed = 2026
)

sim_dgp1 <- sim_data$dgp1
sim_dgp2 <- sim_data$dgp2
sim_dgp3 <- sim_data$dgp3
sim_dgp4 <- sim_data$dgp4

# 시각화 ---------------------------------------------------------------------

# DGP 1 확인 

# True cluster map: HH, LL, MM
ggplot(sim_dgp1) +
  geom_sf(aes(fill = cluster_type),
          color = "grey70",
          linewidth = 0.1) +
  coord_sf(expand = FALSE) +
  labs(fill = "Cluster Types") + 
  # ggtitle("DGP 1: True H-H, L-L, and M-M clusters") +
  theme_void()

# Continuous variable x
ggplot(sim_dgp1) +
  geom_sf(aes(fill = x),
          color = "grey70",
          linewidth = 0.1) +
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("DGP 1: Univariate continuous variable")


# DGP 2 확인 

# x1: CSR-like continuous variable
ggplot(sim_dgp2) +
  geom_sf(aes(fill = x1),
          color = "grey70",
          linewidth = 0.1) +
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("DGP 2: Null mixed-type CSR - x1")

# x2: CSR-like continuous variable
ggplot(sim_dgp2) +
  geom_sf(aes(fill = x2),
          color = "grey70",
          linewidth = 0.1) +
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("DGP 2: Null mixed-type CSR - x2")

# x3: CSR-like continuous variable
ggplot(sim_dgp2) +
  geom_sf(aes(fill = x3),
          color = "grey70",
          linewidth = 0.1) +
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("DGP 2: Null mixed-type CSR - x3")

# cat1: random categorical variable
ggplot(sim_dgp2) +
  geom_sf(aes(fill = cat1),
          color = "grey70",
          linewidth = 0.1) +
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("DGP 2: Null mixed-type CSR - categorical variable")


# DGP 3 확인 

# True cluster map
ggplot(sim_dgp3) +
  geom_sf(aes(fill = cluster_type),
          color = "grey70",
          linewidth = 0.1) +
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("DGP 3: True linear multivariate cluster")

# x1: signal variable
ggplot(sim_dgp3) +
  geom_sf(aes(fill = x1),
          color = "grey70",
          linewidth = 0.1) +
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("DGP 3: x1 with linear cluster")

# x2: signal variable
ggplot(sim_dgp3) +
  geom_sf(aes(fill = x2),
          color = "grey70",
          linewidth = 0.1) +
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("DGP 3: x2 with linear cluster")

# x3: CSR-like noise variable
ggplot(sim_dgp3) +
  geom_sf(aes(fill = x3),
          color = "grey70",
          linewidth = 0.1) +
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("DGP 3: x3 as CSR-like noise")


# DGP 4 확인 

# True cluster map
ggplot(sim_dgp4) +
  geom_sf(aes(fill = cluster_type),
          color = "grey70",
          linewidth = 0.1) +
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("DGP 4: True nonlinear mixed-type interaction cluster")

# x1
ggplot(sim_dgp4) +
  geom_sf(aes(fill = x1),
          color = "grey70",
          linewidth = 0.1) +
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("DGP 4: x1")

# x2
ggplot(sim_dgp4) +
  geom_sf(aes(fill = x2),
          color = "grey70",
          linewidth = 0.1) +
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("DGP 4: x2")

# x1 × x2 interaction
ggplot(sim_dgp4) +
  geom_sf(aes(fill = x1x2),
          color = "grey70",
          linewidth = 0.1) +
  scale_fill_gradient2(midpoint = 0) +
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("DGP 4: Diagnostic interaction x1 × x2")

# Categorical variable
ggplot(sim_dgp4) +
  geom_sf(aes(fill = cat1),
          color = "grey70",
          linewidth = 0.1) +
  coord_sf(expand = FALSE) +
  theme_void() +
  ggtitle("DGP 4: Categorical variable")

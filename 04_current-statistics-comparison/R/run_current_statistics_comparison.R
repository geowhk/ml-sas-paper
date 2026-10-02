suppressPackageStartupMessages({
  library(dplyr)
  library(future)
  library(ggplot2)
  library(readr)
  library(sf)
  library(spdep)
  library(tibble)
  library(tidyr)
  library(yaml)
})

source("03_hyperparameter-sensitivity/R/sensitivity_helpers.R")

표준화_유사도 <- function(S) {
  diag(S) <- 0
  S / sum(S)
}

global_mlsas_계산 <- function(P, Q, 순열) {
  n <- nrow(P)
  n * (n - 1) * sum(P * Q[순열, 순열])
}

moran_i_계산 <- function(y, W) {
  z <- y - mean(y)
  n <- length(y)
  as.numeric((n / sum(W)) * crossprod(z, W %*% z) / crossprod(z))
}

geary_c_계산 <- function(y, W) {
  z <- y - mean(y)
  n <- length(y)
  차이제곱 <- outer(y, y, "-")^2
  as.numeric(((n - 1) / (2 * sum(W))) * sum(W * 차이제곱) / crossprod(z))
}

bb_join_count_계산 <- function(y, A) {
  b <- as.numeric(y == "B")
  as.numeric(0.5 * crossprod(b, A %*% b))
}

상관요약 <- function(x, y, 비교명) {
  검정 <- stats::cor.test(x, y, method = "pearson")
  tibble(
    comparison = 비교명,
    n = sum(stats::complete.cases(x, y)),
    pearson_r = unname(검정$estimate),
    ci_low = unname(검정$conf.int[1]),
    ci_high = unname(검정$conf.int[2]),
    p_value = 검정$p.value
  )
}

그림저장 <- function(결과, 상관, 출력폴더) {
  연속형 <- 결과 |>
    select(permutation_id, continuous_mlsas, moran_i, geary_c) |>
    pivot_longer(
      cols = c(moran_i, geary_c),
      names_to = "comparison",
      values_to = "statistic"
    ) |>
    mutate(
      comparison = recode(
        comparison,
        moran_i = "Moran's I",
        geary_c = "Geary's c"
      ),
      comparison = factor(comparison, levels = c("Moran's I", "Geary's c"))
    )

  연속형_라벨 <- 상관 |>
    filter(comparison %in% c("Moran's I", "Geary's c")) |>
    mutate(
      comparison = factor(comparison, levels = c("Moran's I", "Geary's c")),
      label = sprintf("Pearson r = %.3f", pearson_r)
    )

  p1 <- ggplot(연속형, aes(continuous_mlsas, statistic)) +
    geom_point(size = 1.0, alpha = 0.55) +
    geom_smooth(method = "lm", formula = y ~ x, color = "#2f6df6", linewidth = 0.8) +
    geom_text(
      data = 연속형_라벨,
      aes(x = Inf, y = Inf, label = label),
      hjust = 1.08,
      vjust = 1.4,
      inherit.aes = FALSE,
      size = 3.5
    ) +
    facet_wrap(~ comparison, scales = "free_y", nrow = 1) +
    labs(
      x = "Global ML-SAS (continuous)",
      y = NULL,
      title = "Global ML-SAS and Classical Continuous Statistics",
      subtitle = "20 x 20 hexagonal grid; 1,000 common spatial permutations"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.border = element_rect(color = "grey45", fill = NA, linewidth = 0.45),
      panel.grid.minor = element_blank()
    )

  ggsave(
    file.path(출력폴더, "figures", "continuous_mlsas_moran_geary.png"),
    p1,
    width = 10,
    height = 4.8,
    dpi = 300,
    bg = "white"
  )

  범주형_상관 <- 상관 |>
    filter(comparison == "BB join-count") |>
    mutate(label = sprintf("Pearson r = %.3f", pearson_r))

  p2 <- ggplot(결과, aes(categorical_mlsas, bb_join_count)) +
    geom_point(size = 1.0, alpha = 0.55) +
    geom_smooth(method = "lm", formula = y ~ x, color = "#2f6df6", linewidth = 0.8) +
    geom_text(
      data = 범주형_상관,
      aes(x = Inf, y = Inf, label = label),
      hjust = 1.08,
      vjust = 1.4,
      inherit.aes = FALSE,
      size = 3.5
    ) +
    labs(
      x = "Global ML-SAS (categorical)",
      y = "BB join-count",
      title = "Global ML-SAS and BB Join-Count",
      subtitle = "20 x 20 hexagonal grid; 1,000 common spatial permutations"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.border = element_rect(color = "grey45", fill = NA, linewidth = 0.45),
      panel.grid.minor = element_blank()
    )

  ggsave(
    file.path(출력폴더, "figures", "categorical_mlsas_bb_join_count.png"),
    p2,
    width = 6.2,
    height = 5.2,
    dpi = 300,
    bg = "white"
  )
}

분석실행 <- function(설정파일) {
  설정 <- yaml::read_yaml(설정파일)
  설정$tree_workers <- as.integer(Sys.getenv("MLSAS_WORKERS", as.character(설정$tree_workers)))
  출력폴더 <- if (length(commandArgs(trailingOnly=TRUE)) == 2L) commandArgs(trailingOnly=TRUE)[2] else 설정$output_dir
  if (file.exists(출력폴더)) stop("Output must not exist: ", 출력폴더)
  dir.create(file.path(출력폴더, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(출력폴더, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(출력폴더, "objects"), recursive = TRUE, showWarnings = FALSE)

  격자 <- readRDS(설정$dgp_file)
  n <- nrow(격자)
  if (n != 400L || n %% 2L != 0L) {
    stop("20 x 20 격자와 균형 이진 표본을 위해 관측치 수가 400이어야 한다.", call. = FALSE)
  }

  nb <- make_nb(격자, queen = TRUE)
  listw <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)
  W <- spdep::listw2mat(listw)
  diag(W) <- 0
  A <- spdep::nb2mat(nb, style = "B", zero.policy = TRUE)
  diag(A) <- 0
  if (!isTRUE(all.equal(A, t(A)))) {
    stop("BB join-count requires a symmetric binary adjacency matrix.", call. = FALSE)
  }
  P <- W / sum(W)

  set.seed(as.integer(설정$master_seed))
  연속형_master <- stats::rnorm(n)
  범주형_master <- factor(sample(rep(c("B", "W"), each = n / 2L)), levels = c("B", "W"))

  workers <- as.integer(설정$tree_workers)
  if (workers > 1L) {
    future::plan(future::multisession, workers = workers)
  } else {
    future::plan(future::sequential)
  }
  on.exit(future::plan(future::sequential), add = TRUE)

  message("연속형 master sample의 path proximity를 계산한다.")
  set.seed(as.integer(설정$continuous_urf_seed))
  S_연속형 <- compute_path_proximity_hp(
    data.frame(x = 연속형_master),
    n_tree = as.integer(설정$ntree),
    mtry = as.integer(설정$mtry),
    nodesize = as.integer(설정$nodesize),
    workers = workers
  )

  message("범주형 master sample의 path proximity를 계산한다.")
  set.seed(as.integer(설정$categorical_urf_seed))
  S_범주형 <- compute_path_proximity_hp(
    data.frame(x = 범주형_master),
    n_tree = as.integer(설정$ntree),
    mtry = as.integer(설정$mtry),
    nodesize = as.integer(설정$nodesize),
    workers = workers
  )

  Q_연속형 <- 표준화_유사도(S_연속형)
  Q_범주형 <- 표준화_유사도(S_범주형)

  nperm <- as.integer(설정$permutations)
  set.seed(as.integer(설정$permutation_seed))
  순열행렬 <- t(replicate(nperm, sample.int(n), simplify = "matrix"))

  message("공통 공간 순열 ", nperm, "개에 대한 통계량을 계산한다.")
  결과목록 <- vector("list", nperm)
  for (b in seq_len(nperm)) {
    순열 <- 순열행렬[b, ]
    연속형_y <- 연속형_master[순열]
    범주형_y <- 범주형_master[순열]

    결과목록[[b]] <- tibble(
      permutation_id = b,
      continuous_mlsas = global_mlsas_계산(P, Q_연속형, 순열),
      moran_i = moran_i_계산(연속형_y, W),
      geary_c = geary_c_계산(연속형_y, W),
      categorical_mlsas = global_mlsas_계산(P, Q_범주형, 순열),
      bb_join_count = bb_join_count_계산(범주형_y, A)
    )
  }
  결과 <- bind_rows(결과목록)

  상관 <- bind_rows(
    상관요약(결과$continuous_mlsas, 결과$moran_i, "Moran's I"),
    상관요약(결과$continuous_mlsas, 결과$geary_c, "Geary's c"),
    상관요약(결과$categorical_mlsas, 결과$bb_join_count, "BB join-count")
  )

  메타데이터 <- tibble(
    item = c(
      "n", "grid", "mlsas_weights_style", "join_count_weights_style",
      "total_directed_links", "total_undirected_links", "ntree", "mtry",
      "nodesize", "tree_workers", "permutations", "master_seed",
      "continuous_urf_seed", "categorical_urf_seed", "permutation_seed"
    ),
    value = as.character(c(
      n, "20 x 20 hexagonal", "W", "B (symmetric binary adjacency)",
      sum(card(nb)), sum(A) / 2, 설정$ntree, 설정$mtry,
      설정$nodesize, workers, nperm, 설정$master_seed,
      설정$continuous_urf_seed, 설정$categorical_urf_seed, 설정$permutation_seed
    ))
  )

  readr::write_csv(결과, file.path(출력폴더, "tables", "permutation_statistics.csv"))
  readr::write_csv(상관, file.path(출력폴더, "tables", "correlation_summary.csv"))
  readr::write_csv(메타데이터, file.path(출력폴더, "tables", "run_metadata.csv"))
  saveRDS(
    list(
      config = 설정,
      master_samples = list(continuous = 연속형_master, categorical = 범주형_master),
      spatial = list(nb = nb, W = W, A = A),
      proximity = list(continuous = S_연속형, categorical = S_범주형),
      permutations = 순열행렬,
      results = 결과,
      correlations = 상관
    ),
    file.path(출력폴더, "objects", "current_statistics_comparison.rds")
  )

  그림저장(결과, 상관, 출력폴더)
  print(상관)
  message("Current-statistics comparison completed: ", 출력폴더)
  invisible(결과)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) %in% c(1L,2L)) {
    stop("사용법: Rscript 04_current-statistics-comparison/R/run_current_statistics_comparison.R <config.yml>", call. = FALSE)
  }

  분석실행(args[[1]])
}

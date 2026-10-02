# DGP 등록 안내
#
# 새 시뮬레이션 시나리오는 dgp_catalog()에 추가한다.
# 1. 시나리오 .rds 파일을 02_revision-simulation/sim_dgp/에 넣는다.
# 2. dgp5 = list(...) 같은 새 항목을 추가한다. 필요한 필드는 다음과 같다.
#    - file: .rds 파일명.
#    - label: 사람이 읽기 쉬운 짧은 시나리오 설명.
#    - continuous_sets: Multivariate Local Geary와 ML-SAS에 사용할 연속형 변수 세트.
#    - mixed_sets: ML-SAS에만 사용할 mixed-type 변수 세트.
#    - compare_univariate_lisa: 단일 연속형 변수 시나리오에서 Moran, Geary,
#      Getis-Ord 비교도 함께 돌릴 때만 TRUE로 둔다.
# 3. 새 시나리오는 아래 명령으로 가벼운 실행 검증을 한다.
#    Rscript 02_revision-simulation/R/run_compare_analysis.R \
#      --config 02_revision-simulation/config/smoke.yml --dgp_id dgp5
# 4. AWS 일괄 실행에 포함하려면 dgp id를 아래 파일에 추가한다.
#    02_revision-simulation/config/aws_all.yml.
#
# .rds 객체 조건:
# - sf 객체여야 한다.
# - 가능하면 id, true_cluster, cluster_type 컬럼을 포함한다.
# - id가 없으면 1:n으로 자동 생성한다.
# - true_cluster 또는 cluster_type이 없으면 보수적인 기본값을 추가한다.

dgp_catalog <- function() {
  list(
    dgp1 = list(
      file = "sim_dgp1.rds",
      label = "DGP 1: univariate continuous Geary-type clusters",
      continuous_sets = list(univariate = "x"),
      mixed_sets = list(),
      compare_univariate_lisa = TRUE
    ),
    dgp2 = list(
      file = "sim_dgp2.rds",
      label = "DGP 2: null mixed-type CSR",
      continuous_sets = list(continuous = c("x1", "x2", "x3")),
      mixed_sets = list(mixed = c("x1", "x2", "x3", "cat1")),
      compare_univariate_lisa = FALSE
    ),
    dgp3 = list(
      file = "sim_dgp3.rds",
      label = "DGP 3: linear multivariate continuous cluster",
      continuous_sets = list(continuous = c("x1", "x2", "x3")),
      mixed_sets = list(),
      compare_univariate_lisa = FALSE
    ),
    dgp3.1 = list(
      file = "sim_dgp3.1.rds",
      label = "DGP 3.1: linear multivariate continuous cluster",
      continuous_sets = list(continuous = c("x1", "x2", "x3")),
      mixed_sets = list(),
      compare_univariate_lisa = FALSE
    ),
    dgp4 = list(
      file = "sim_dgp4.rds",
      label = "DGP 4: nonlinear mixed-type interaction cluster",
      continuous_sets = list(continuous = c("x1", "x2")),
      mixed_sets = list(mixed = c("x1", "x2", "cat1")),
      compare_univariate_lisa = FALSE
    ),
    dgp4.1 = list(
      file = "sim_dgp4.1.rds",
      label = "DGP 4.1: nonlinear mixed-type interaction cluster",
      continuous_sets = list(continuous = c("x1", "x2", "x3")),
      mixed_sets = list(mixed = c("x1", "x2", "x3", "cat1")),
      compare_univariate_lisa = FALSE
    ),
    dgp4.2 = list(
      file = "sim_dgp4.2.rds",
      label = "DGP 4.2: nonlinear mixed-type interaction cluster",
      continuous_sets = list(continuous = c("x1", "x2")),
      mixed_sets = list(mixed = c("x1", "x2", "cat1")),
      compare_univariate_lisa = FALSE
    ),
    dgp4.3 = list(
      file = "sim_dgp4.3.rds",
      label = "DGP 4.3: nonlinear mixed-type interaction cluster",
      continuous_sets = list(continuous = c("x1", "x2")),
      mixed_sets = list(mixed = c("x1", "x2", "cat1")),
      compare_univariate_lisa = FALSE
    ),
    dgp4.4 = list(
      file = "sim_dgp4.4.rds",
      label = "DGP 4.4: nonlinear mixed-type interaction cluster",
      continuous_sets = list(continuous = c("x1", "x2")),
      mixed_sets = list(mixed = c("x1", "x2", "cat1")),
      compare_univariate_lisa = FALSE
    )
  )
}

load_dgp <- function(dgp_id, input_dir = "02_revision-simulation/sim_dgp") {
  catalog <- dgp_catalog()
  dgp_id <- match.arg(dgp_id, names(catalog))
  spec <- catalog[[dgp_id]]
  path <- file.path(input_dir, spec$file)

  if (!file.exists(path)) {
    stop("DGP file not found: ", path, call. = FALSE)
  }

  sf_obj <- readRDS(path)
  if (!inherits(sf_obj, "sf")) {
    stop("Expected an sf object in: ", path, call. = FALSE)
  }
  if (!"id" %in% names(sf_obj)) {
    sf_obj$id <- seq_len(nrow(sf_obj))
  }
  if (!"true_cluster" %in% names(sf_obj)) {
    sf_obj$true_cluster <- FALSE
  }
  if (!"cluster_type" %in% names(sf_obj)) {
    sf_obj$cluster_type <- ifelse(sf_obj$true_cluster, "cluster", "background")
  }

  list(id = dgp_id, spec = spec, sf = sf_obj)
}

load_professor_dgp <- load_dgp

prepare_attribute_data <- function(sf_obj, variables) {
  data_df <- sf::st_drop_geometry(sf_obj)[, variables, drop = FALSE]
  data_df <- as.data.frame(data_df)

  for (nm in names(data_df)) {
    if (is.character(data_df[[nm]]) || is.logical(data_df[[nm]])) {
      data_df[[nm]] <- factor(data_df[[nm]])
    }
  }

  rownames(data_df) <- as.character(sf_obj$id)
  data_df
}

make_nb <- function(sf_obj, queen = TRUE) {
  spdep::poly2nb(sf_obj, queen = queen, snap = sqrt(.Machine$double.eps))
}

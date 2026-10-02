suppressPackageStartupMessages({
  library(sf)
  library(tidyverse)
  library(spdep)
  library(future)
  library(future.apply)
  library(randomForest)
  library(rgeoda)
  library(forcats)
  library(patchwork)
})

source("02_revision-simulation/R/mlsas_helpers.R")
source("02_revision-simulation/R/multigeary_helpers.R")
source("02_revision-simulation/R/dgp.R")

default_analysis_config <- function(
    dgp_id = "dgp4",
    input_dir = "02_revision-simulation/sim_dgp",
    output_root = "02_revision-simulation/runs",
    seed = 123,
    rf_ntree = 100,
    rf_mtry = 1,
    rf_nodesize = 1,
    perm_nsim = 999,
    workers = 1,
    alpha = 0.05,
    mlsas_method = "permute_cond",
    mlsas_alternative = "two.sided",
    multigeary_alternative = "two.sided",
    local_geary_alternative = "two.sided",
    moran_alternative = "two.sided",
    require_moran_cluster_direction = FALSE,
    getis_ord_inference = "analytic",
    getis_ord_alternative = "two.sided",
    style = "W",
    diag_zero = TRUE,
    queen = TRUE,
    run_univariate_lisa = TRUE,
    save_similarity_matrix = FALSE,
    continuous_variables = NULL,
    continuous_set_name = "x1_x2",
    include_mixed_sets = TRUE,
    save_diagnostic_figures = TRUE,
    run_tag = NULL
) {
  list(
    dgp_id = dgp_id,
    input_dir = input_dir,
    output_root = output_root,
    seed = seed,
    rf_ntree = as.integer(rf_ntree),
    rf_mtry = as.integer(rf_mtry),
    rf_nodesize = as.integer(rf_nodesize),
    perm_nsim = as.integer(perm_nsim),
    workers = as.integer(workers),
    alpha = alpha,
    mlsas_method = mlsas_method,
    mlsas_alternative = match.arg(mlsas_alternative, c("greater", "less", "two.sided")),
    multigeary_alternative = match.arg(multigeary_alternative, c("less", "two.sided")),
    # spdep::localC_perm() names the positive-spatial-autocorrelation test
    # "greater", even though positive similarity corresponds to a small raw
    # local C value. Keep the package's API value in metadata for traceability.
    local_geary_alternative = match.arg(local_geary_alternative, c("greater", "less", "two.sided")),
    moran_alternative = match.arg(moran_alternative, c("greater", "less", "two.sided")),
    require_moran_cluster_direction = isTRUE(require_moran_cluster_direction),
    getis_ord_inference = match.arg(getis_ord_inference, c("analytic", "permutation")),
    getis_ord_alternative = match.arg(getis_ord_alternative, c("greater", "less", "two.sided")),
    style = style,
    diag_zero = diag_zero,
    queen = queen,
    run_univariate_lisa = run_univariate_lisa,
    save_similarity_matrix = save_similarity_matrix,
    continuous_variables = if (is.null(continuous_variables)) {
      NULL
    } else {
      as.character(unlist(continuous_variables, use.names = FALSE))
    },
    continuous_set_name = continuous_set_name,
    include_mixed_sets = include_mixed_sets,
    save_diagnostic_figures = save_diagnostic_figures,
    run_tag = run_tag
  )
}

create_run_dirs <- function(config) {
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  tagged_dgp_id <- config$dgp_id
  if (!is.null(config$run_tag) && nzchar(config$run_tag)) {
    safe_tag <- gsub("[^[:alnum:]_.-]+", "_", config$run_tag)
    tagged_dgp_id <- paste(config$dgp_id, safe_tag, sep = "_")
  }
  run_id <- sprintf(
    "%s_ntree%s_perm%s_%s",
    tagged_dgp_id,
    config$rf_ntree,
    config$perm_nsim,
    stamp
  )
  run_dir <- file.path(config$output_root, run_id)
  dirs <- list(
    run = run_dir,
    tables = file.path(run_dir, "tables"),
    figures = file.path(run_dir, "figures"),
    objects = file.path(run_dir, "objects"),
    logs = file.path(run_dir, "logs")
  )
  purrr::walk(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
  dirs
}

write_run_metadata <- function(config, dgp_obj, dirs) {
  metadata <- list(
    config = config,
    dgp_label = dgp_obj$spec$label,
    dgp_columns = names(dgp_obj$sf),
    dgp_n = nrow(dgp_obj$sf),
    cluster_table = as.list(table(dgp_obj$sf$cluster_type, useNA = "ifany")),
    created_at = as.character(Sys.time())
  )
  saveRDS(metadata, file.path(dirs$objects, "run_metadata.rds"))
  writeLines(
    capture.output(str(metadata, max.level = 3)),
    con = file.path(dirs$logs, "run_metadata.txt")
  )
}

classify_mlsas_clusters <- function(local_result, alpha = 0.05, alternative = "two.sided") {
  alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
  exploratory_p_col <- if ("p_perm" %in% names(local_result)) "p_perm" else "P"
  strict_p_col <- if ("p_perm_adj" %in% names(local_result)) {
    "p_perm_adj"
  } else if ("P.adj" %in% names(local_result)) {
    "P.adj"
  } else {
    "P"
  }

  has_direction_columns <- all(c("local.Moran", "E") %in% names(local_result))
  if (alternative == "greater" && !has_direction_columns) {
    stop("Upper-tail ML-SAS classification requires `local.Moran` and `E` columns.", call. = FALSE)
  }

  local_result |>
    mutate(
      exploratory_significance_source = exploratory_p_col,
      strict_significance_source = strict_p_col,
      similarity_direction = if (has_direction_columns) local.Moran > E else TRUE,
      direction_rule = if (alternative == "greater") "observed > permutation expectation" else "not applied",
      cluster_flag_exploratory = if_else(
        .data[[exploratory_p_col]] <= alpha &
          (alternative != "greater" | similarity_direction),
        "Significant spatial cluster",
        "Not significant"
      ),
      cluster_flag_strict = if_else(
        .data[[strict_p_col]] <= alpha &
          (alternative != "greater" | similarity_direction),
        "Significant spatial cluster",
        "Not significant"
      ),
      cluster_flag_exploratory = factor(
        cluster_flag_exploratory,
        levels = c("Significant spatial cluster", "Not significant")
      ),
      cluster_flag_strict = factor(
        cluster_flag_strict,
        levels = c("Significant spatial cluster", "Not significant")
      )
    )
}

summarise_detection <- function(result_sf, method_name, significance_col) {
  sf::st_drop_geometry(result_sf) |>
    mutate(true_cluster = as.logical(true_cluster)) |>
    group_by(true_cluster) |>
    summarise(
      n_cells = n(),
      detected_cells = sum(.data[[significance_col]] == "Significant spatial cluster", na.rm = TRUE),
      detection_rate = detected_cells / n_cells,
      .groups = "drop"
    ) |>
    mutate(method = method_name, significance_rule = significance_col) |>
    relocate(method, significance_rule, true_cluster)
}

make_binary_map <- function(sf_data, fill_col, title_text) {
  palette <- c(
    "Significant spatial cluster" = "#d7301f",
    "Not significant" = "#eeeeee"
  )

  ggplot(sf_data) +
    geom_sf(aes(fill = .data[[fill_col]]), color = "grey75", linewidth = 0.1) +
    scale_fill_manual(values = palette, drop = FALSE) +
    labs(title = title_text, fill = NULL) +
    theme_void(base_size = 11) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.key.size = unit(0.35, "cm"),
      plot.title = element_text(
        face = "bold",
        size = 9,
        lineheight = 0.95,
        margin = ggplot2::margin(b = 4)
      )
    )
}

format_set_label <- function(set_name) {
  dplyr::case_when(
    set_name == "continuous_interaction" ~ "cont. interaction",
    set_name == "mixed_interaction" ~ "mixed interaction",
    TRUE ~ gsub("_", " ", set_name)
  )
}

format_map_title <- function(method_name, set_name, rule_name = NULL) {
  parts <- c(method_name, rule_name, format_set_label(set_name))
  paste(parts[!is.na(parts) & nzchar(parts)], collapse = ": ")
}

run_mlsas_for_variable_set <- function(dgp_obj, nb_obj, variables, set_name, config, dirs) {
  message("Running ML-SAS for ", dgp_obj$id, " / ", set_name, ": ", paste(variables, collapse = ", "))
  data_df <- prepare_attribute_data(dgp_obj$sf, variables)

  S <- compute_path_proximity(
    data = data_df,
    n_tree = config$rf_ntree,
    workers = config$workers,
    validate = FALSE,
    mtry = config$rf_mtry,
    nodesize = config$rf_nodesize
  )

  local_res <- test.local.mlsas(
    id = dgp_obj$sf$id,
    nb = nb_obj,
    S = S,
    style = config$style,
    method = config$mlsas_method,
    alternative = config$mlsas_alternative,
    diag.zero = config$diag_zero,
    p.adj = "fdr",
    nperm = config$perm_nsim
  ) |>
    classify_mlsas_clusters(alpha = config$alpha, alternative = config$mlsas_alternative)

  result_sf <- bind_cols(dgp_obj$sf, local_res |> select(-id))
  prefix <- paste(dgp_obj$id, set_name, "mlsas", sep = "_")

  readr::write_csv(
    sf::st_drop_geometry(result_sf),
    file.path(dirs$tables, paste0(prefix, "_local_results.csv"))
  )
  saveRDS(
    list(local = local_res, sf = result_sf, variables = variables),
    file.path(dirs$objects, paste0(prefix, "_object.rds"))
  )
  if (isTRUE(config$save_similarity_matrix)) {
    saveRDS(S, file.path(dirs$objects, paste0(prefix, "_similarity_matrix.rds")))
  }

  list(
    set_name = set_name,
    variables = variables,
    local = local_res,
    sf = result_sf,
    exploratory_summary = summarise_detection(
      result_sf,
      paste0("Local ML-SAS (", set_name, ")"),
      "cluster_flag_exploratory"
    ),
    strict_summary = summarise_detection(
      result_sf,
      paste0("Local ML-SAS (", set_name, ")"),
      "cluster_flag_strict"
    ),
    exploratory_map = make_binary_map(
      result_sf,
      "cluster_flag_exploratory",
      format_map_title("ML-SAS", set_name, "raw")
    ),
    strict_map = make_binary_map(
      result_sf,
      "cluster_flag_strict",
      format_map_title("ML-SAS", set_name, "FDR")
    )
  )
}

run_multigeary_for_variable_set <- function(dgp_obj, variables, set_name, config, dirs) {
  message("Running Multivariate Local Geary for ", dgp_obj$id, " / ", set_name)
  w <- rgeoda::queen_weights(dgp_obj$sf)
  raw <- rgeoda::local_multigeary(
    w,
    dgp_obj$sf[, variables],
    permutations = config$perm_nsim,
    significance_cutoff = config$alpha
  )
  processed <- snu_multigeary(
    raw,
    alternative = config$multigeary_alternative,
    p.adj.type = "fdr",
    q.value = config$alpha
  ) |>
    classify_multigeary_clusters()

  result_sf <- bind_cols(dgp_obj$sf, processed)
  prefix <- paste(dgp_obj$id, set_name, "multigeary", sep = "_")
  readr::write_csv(
    sf::st_drop_geometry(result_sf),
    file.path(dirs$tables, paste0(prefix, "_local_results.csv"))
  )
  saveRDS(
    list(raw = raw, processed = processed, sf = result_sf, variables = variables),
    file.path(dirs$objects, paste0(prefix, "_object.rds"))
  )

  list(
    set_name = set_name,
    variables = variables,
    processed = processed,
    sf = result_sf,
    summary = summarise_detection(
      result_sf,
      paste0("Multivariate Local Geary (", set_name, ")"),
      "cluster_flag"
    ),
    map = make_binary_map(
      result_sf,
      "cluster_flag",
      format_map_title("Multigeary", set_name, "FDR positive")
    )
  )
}

extract_simulation_p <- function(x, source_name) {
  x <- as.matrix(x)
  p_col <- grep("Pr\\(.*\\) Sim$", colnames(x), value = TRUE)
  p_col <- setdiff(p_col, grep("folded", p_col, value = TRUE, ignore.case = TRUE))
  if (length(p_col) != 1L) {
    stop("Could not identify one simulation p-value column for ", source_name, call. = FALSE)
  }
  as.numeric(x[, p_col])
}

run_univariate_lisa_comparison <- function(dgp_obj, nb_obj, config, dirs) {
  if (!isTRUE(dgp_obj$spec$compare_univariate_lisa) || !isTRUE(config$run_univariate_lisa)) {
    return(NULL)
  }
  if (!"x" %in% names(dgp_obj$sf)) {
    return(NULL)
  }

  message("Running univariate LISA comparison for ", dgp_obj$id)
  listw <- spdep::nb2listw(nb_obj, style = "W", zero.policy = TRUE)
  listw_self <- spdep::nb2listw(spdep::include.self(nb_obj), style = "W", zero.policy = TRUE)
  x <- dgp_obj$sf$x

  moran <- spdep::localmoran_perm(
    x,
    listw,
    nsim = config$perm_nsim,
    alternative = config$moran_alternative,
    iseed = config$seed
  )
  geary <- spdep::localC_perm(
    x,
    listw,
    nsim = config$perm_nsim,
    alternative = config$local_geary_alternative,
    iseed = config$seed
  )
  gi <- if (config$getis_ord_inference == "permutation") {
    spdep::localG_perm(
      x,
      listw_self,
      nsim = config$perm_nsim,
      alternative = config$getis_ord_alternative,
      iseed = config$seed
    )
  } else {
    spdep::localG(x, listw_self, zero.policy = TRUE)
  }
  gi_internals <- attr(gi, "internals")

  moran_p <- extract_simulation_p(moran, "Local Moran")
  geary_p <- extract_simulation_p(attr(geary, "pseudo-p"), "Local Geary")
  gi_p <- if (config$getis_ord_inference == "permutation") {
    extract_simulation_p(gi_internals, "Getis-Ord Gi*")
  } else {
    as.numeric(gi_internals[, "Pr(z != E(Gi))"])
  }
  z <- x - mean(x)
  spatial_lag_z <- spdep::lag.listw(listw, z, zero.policy = TRUE)
  moran_similarity_direction <- z * spatial_lag_z > 0

  out <- tibble(
    id = dgp_obj$sf$id,
    moran_p = moran_p,
    geary_p = geary_p,
    gi_z = gi_internals[, 4],
    gi_p = gi_p,
    moran_similarity_direction = moran_similarity_direction,
    getis_ord_inference = config$getis_ord_inference
  ) |>
    mutate(
      moran_p_adj = p.adjust(moran_p, method = "fdr"),
      geary_p_adj = p.adjust(geary_p, method = "fdr"),
      gi_p_adj = p.adjust(gi_p, method = "fdr"),
      moran_flag = if_else(
        moran_p_adj <= config$alpha &
          (!config$require_moran_cluster_direction | moran_similarity_direction),
        "Significant spatial cluster",
        "Not significant"
      ),
      geary_flag = if_else(geary_p_adj <= config$alpha, "Significant spatial cluster", "Not significant"),
      gi_flag = if_else(gi_p_adj <= config$alpha, "Significant spatial cluster", "Not significant")
    )

  result_sf <- bind_cols(dgp_obj$sf, out |> select(-id))
  readr::write_csv(
    sf::st_drop_geometry(result_sf),
    file.path(dirs$tables, paste0(dgp_obj$id, "_univariate_lisa_results.csv"))
  )

  out
}

save_map_grid <- function(maps, file_path) {
  if (!length(maps)) {
    return(invisible(NULL))
  }

  ncol <- min(3, length(maps))
  fig <- patchwork::wrap_plots(maps, ncol = ncol, guides = "collect") &
    theme(legend.position = "bottom")
  ggsave(
    file_path,
    fig,
    width = 5 * ncol,
    height = 4 * ceiling(length(maps) / ncol),
    dpi = 300,
    bg = "white"
  )
}

save_comparison_figures <- function(dgp_obj, mlsas_results, geary_results, dirs, config) {
  comparison_maps <- list()
  for (res in geary_results) {
    comparison_maps[[paste0("geary_", res$set_name)]] <- res$map
  }
  for (res in mlsas_results) {
    comparison_maps[[paste0("mlsas_strict_", res$set_name)]] <- res$strict_map
  }

  save_map_grid(
    comparison_maps,
    file.path(dirs$figures, paste0(dgp_obj$id, "_comparison_maps.png"))
  )
  if (isTRUE(config$save_diagnostic_figures)) {
    diagnostic_maps <- list()
    for (res in geary_results) {
      diagnostic_maps[[paste0("geary_", res$set_name)]] <- res$map
    }
    for (res in mlsas_results) {
      diagnostic_maps[[paste0("mlsas_exploratory_", res$set_name)]] <- res$exploratory_map
      diagnostic_maps[[paste0("mlsas_strict_", res$set_name)]] <- res$strict_map
    }
    save_map_grid(
      diagnostic_maps,
      file.path(dirs$figures, paste0(dgp_obj$id, "_diagnostic_maps.png"))
    )
  }

  invisible(NULL)
}

run_dgp_analysis <- function(config) {
  set.seed(config$seed)
  if (config$workers > 1) {
    future::plan(future::multisession, workers = config$workers)
  } else {
    future::plan(future::sequential)
  }
  on.exit(future::plan(future::sequential), add = TRUE)

  dgp_obj <- load_dgp(config$dgp_id, input_dir = config$input_dir)
  dirs <- create_run_dirs(config)
  write_run_metadata(config, dgp_obj, dirs)
  nb_obj <- make_nb(dgp_obj$sf, queen = config$queen)
  saveRDS(nb_obj, file.path(dirs$objects, paste0(config$dgp_id, "_nb_object.rds")))

  continuous_sets <- dgp_obj$spec$continuous_sets
  mixed_sets <- dgp_obj$spec$mixed_sets
  if (!is.null(config$continuous_variables)) {
    missing_variables <- setdiff(config$continuous_variables, names(dgp_obj$sf))
    if (length(missing_variables)) {
      stop(
        "DGP에 존재하지 않는 연속형 변수가 지정됨: ",
        paste(missing_variables, collapse = ", "),
        call. = FALSE
      )
    }
    continuous_sets <- stats::setNames(
      list(config$continuous_variables),
      config$continuous_set_name
    )
  }
  if (!isTRUE(config$include_mixed_sets)) {
    mixed_sets <- list()
  }

  geary_results <- purrr::imap(
    continuous_sets,
    ~ run_multigeary_for_variable_set(dgp_obj, .x, .y, config, dirs)
  )
  mlsas_cont_results <- purrr::imap(
    continuous_sets,
    ~ run_mlsas_for_variable_set(dgp_obj, nb_obj, .x, .y, config, dirs)
  )
  mlsas_mixed_results <- purrr::imap(
    mixed_sets,
    ~ run_mlsas_for_variable_set(dgp_obj, nb_obj, .x, .y, config, dirs)
  )
  lisa_results <- run_univariate_lisa_comparison(dgp_obj, nb_obj, config, dirs)

  detection_summary <- bind_rows(c(
    purrr::map(geary_results, "summary"),
    purrr::map(mlsas_cont_results, "exploratory_summary"),
    purrr::map(mlsas_cont_results, "strict_summary"),
    purrr::map(mlsas_mixed_results, "exploratory_summary"),
    purrr::map(mlsas_mixed_results, "strict_summary")
  ))

  readr::write_csv(detection_summary, file.path(dirs$tables, paste0(config$dgp_id, "_detection_summary.csv")))
  save_comparison_figures(
    dgp_obj,
    c(mlsas_cont_results, mlsas_mixed_results),
    geary_results,
    dirs,
    config
  )

  result <- list(
    config = config,
    dirs = dirs,
    dgp = dgp_obj,
    geary_results = geary_results,
    mlsas_results = c(mlsas_cont_results, mlsas_mixed_results),
    univariate_lisa = lisa_results,
    detection_summary = detection_summary
  )
  saveRDS(result, file.path(dirs$objects, paste0(config$dgp_id, "_analysis_result.rds")))
  result
}

run_professor_dgp_analysis <- run_dgp_analysis

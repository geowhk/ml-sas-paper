args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L)) {
  stop("Usage: Rscript 02_revision-simulation/R/run_dgp43_lovo.R REPO_DIR RUN_DIR [CONFIG_FILE]")
}

repo_dir <- normalizePath(args[1], mustWork = TRUE)
run_dir <- normalizePath(args[2], mustWork = FALSE)
if (file.exists(run_dir)) stop("Output must not exist: ", run_dir)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)

config_file <- if (length(args) == 3L) {
  normalizePath(args[3], mustWork = TRUE)
} else {
  file.path(repo_dir, "02_revision-simulation", "config", "dgp43_lovo_revision_v2.R")
}
config <- source(config_file)$value
config$workers <- as.integer(Sys.getenv("MLSAS_WORKERS", as.character(config$workers)))
source(file.path(repo_dir, "02_revision-simulation", "R", "mlsas_helpers.R"))
source(file.path(repo_dir, "02_revision-simulation", "R", "mlsas_optimized.R"))
source(file.path(repo_dir, "02_revision-simulation", "R", "mlsas_inference_optimized.R"))
source(file.path(repo_dir, "02_revision-simulation", "R", "dgp.R"))

suppressPackageStartupMessages({
  library(sf)
  library(spdep)
  library(randomForest)
  library(Matrix)
  library(future)
  library(future.apply)
  library(dplyr)
  library(ggplot2)
})

dirs <- file.path(run_dir, c("tables", "figures", "objects", "metadata", "logs"))
names(dirs) <- c("tables", "figures", "objects", "metadata", "logs")
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
progress_file <- file.path(run_dir, "progress.csv")

write_status <- function(stage, detail = "") {
  row <- data.frame(
    time_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stage = stage,
    detail = detail
  )
  utils::write.table(
    row, progress_file, sep = ",", row.names = FALSE,
    col.names = !file.exists(progress_file), append = file.exists(progress_file), quote = TRUE
  )
  message(sprintf("[%s] %s", stage, detail))
}

auc_rank <- function(scores, truth) {
  truth <- as.logical(truth)
  n1 <- sum(truth)
  n0 <- sum(!truth)
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

make_addcl1 <- function(x, seed) {
  set.seed(seed)
  as.data.frame(lapply(x, function(col) sample(col, replace = TRUE)))
}

fit_similarity <- function(x, synthetic, seed, ntree, mtry, nodesize, validate = FALSE) {
  set.seed(seed)
  compute_path_proximity_optimized(
    data = x,
    n_tree = ntree,
    mtry = mtry,
    nodesize = nodesize,
    workers = config$workers,
    validate = validate,
    synthetic_data = synthetic
  )
}

summarise_variant <- function(name, variables, observed, local_inference = NULL) {
  out <- data.frame(
    specification = name,
    variables = paste(variables, collapse = ";"),
    n_variables = length(variables),
    global_mlsas = observed$global,
    true_cluster_mean = mean(observed$local[truth]),
    background_mean = mean(observed$local[!truth]),
    true_minus_background = mean(observed$local[truth]) - mean(observed$local[!truth]),
    local_stat_auc = auc_rank(observed$local, truth)
  )
  if (!is.null(local_inference)) {
    cls <- local_inference$classification
    out$true_similarity_detected <- sum(cls[truth] == "Similarity cluster")
    out$true_any_significant <- sum(cls[truth] != "Not significant")
    out$background_similarity_false_positive <- sum(cls[!truth] == "Similarity cluster")
    out$background_any_significant <- sum(cls[!truth] != "Not significant")
    out$true_similarity_detection_rate <- mean(cls[truth] == "Similarity cluster")
    out$background_similarity_false_positive_rate <- mean(cls[!truth] == "Similarity cluster")
  }
  out
}

write_status("start", config$run_label)
saveRDS(config, file.path(dirs["metadata"], "config.rds"))
dput(config, file = file.path(dirs["metadata"], "config.R"))
capture.output(sessionInfo(), file = file.path(dirs["metadata"], "sessionInfo.txt"))

dgp <- load_dgp(config$dgp_id, file.path(repo_dir, "02_revision-simulation", "sim_dgp"))
if (!identical(config$variables, c("x1", "x2", "cat1"))) {
  stop("DGP 4.3 LOVO expects x1, x2, and cat1 in that order.")
}
data_full <- prepare_attribute_data(dgp$sf, config$variables)
truth <- as.logical(dgp$sf$true_cluster)
nb <- make_nb(dgp$sf, queen = config$queen)
variants <- list(
  full = config$variables,
  leave_x1_out = c("x2", "cat1"),
  leave_x2_out = c("x1", "cat1"),
  leave_cat1_out = c("x1", "x2")
)
saveRDS(list(dgp = dgp$sf, nb = nb, variants = variants), file.path(dirs["objects"], "analysis_inputs.rds"))

permutations <- NULL
if (isTRUE(config$run_primary_inference)) {
  write_status("permutations", sprintf("n=%d nperm=%d seed=%d", nrow(data_full), config$nperm, config$permutation_seed))
  permutations <- make_permutation_indices(nrow(data_full), config$nperm, config$permutation_seed)
  saveRDS(permutations, file.path(dirs["objects"], "permutation_indices.rds"), compress = "xz")
} else {
  write_status("permutations_skipped", "surface-sensitivity run; no local significance claims")
}

primary_synthetic <- make_addcl1(data_full, config$primary_synthetic_seed)
primary_objects <- list()
primary_summary <- list()
primary_local <- list()

for (name in names(variants)) {
  variables <- variants[[name]]
  write_status("primary_fit", sprintf("%s: %s", name, paste(variables, collapse = ",")))
  similarity <- fit_similarity(
    data_full[, variables, drop = FALSE],
    primary_synthetic[, variables, drop = FALSE],
    config$primary_rf_seed,
    config$rf_ntree,
    config$rf_mtry,
    config$rf_nodesize,
    validate = TRUE
  )
  saveRDS(similarity, file.path(dirs["objects"], paste0("similarity_", name, ".rds")))
  observed <- mlsas_observed_optimized(nb, similarity, config$style)
  inference <- NULL
  if (isTRUE(config$run_primary_inference)) {
    write_status("primary_inference", sprintf("%s: %d permutations", name, config$nperm))
    inference <- mlsas_permutation_inference_optimized(
      nb = nb,
      similarity = similarity,
      permutations = permutations,
      workers = config$workers,
      style = config$style,
      alpha = config$alpha,
      p_adjust_method = config$p_adjust_method
    )
  }
  local_table <- cbind(
    sf::st_drop_geometry(dgp$sf)[, c("id", "true_cluster", "cluster_type")],
    specification = name,
    local_mlsas = observed$local
  )
  if (!is.null(inference)) {
    local_table <- cbind(local_table, inference$local[, setdiff(names(inference$local), "local_mlsas"), drop = FALSE])
  }
  utils::write.csv(local_table, file.path(dirs["tables"], paste0("primary_local_", name, ".csv")), row.names = FALSE)
  variant_summary <- summarise_variant(
    name,
    variables,
    observed,
    if (is.null(inference)) NULL else inference$local
  )
  if (!is.null(inference)) {
    variant_summary <- cbind(variant_summary, inference$global)
  }
  primary_summary[[name]] <- variant_summary
  primary_local[[name]] <- observed$local
  primary_objects[[name]] <- list(
    observed = observed,
    global = if (is.null(inference)) NULL else inference$global,
    local = if (is.null(inference)) NULL else inference$local
  )
  rm(similarity)
  invisible(gc())
}

primary_summary <- bind_rows(primary_summary)
utils::write.csv(primary_summary, file.path(dirs["tables"], "primary_inference_summary.csv"), row.names = FALSE)
saveRDS(primary_objects, file.path(dirs["objects"], "primary_inference_results.rds"), compress = "xz")

full_primary <- primary_objects$full
primary_lovo <- lapply(setdiff(names(variants), "full"), function(name) {
  reduced <- primary_objects[[name]]
  local_delta <- full_primary$observed$local - reduced$observed$local
  row <- data.frame(
    removed_variable = sub("leave_(.*)_out", "\\1", name),
    reduced_specification = name,
    full_global_mlsas = full_primary$observed$global,
    reduced_global_mlsas = reduced$observed$global,
    global_change_when_included = full_primary$observed$global - reduced$observed$global,
    mean_local_change_when_included = mean(local_delta),
    mean_absolute_local_change = mean(abs(local_delta)),
    true_cluster_mean_change = mean(local_delta[truth]),
    background_mean_change = mean(local_delta[!truth]),
    true_cluster_mean_absolute_change = mean(abs(local_delta[truth])),
    background_mean_absolute_change = mean(abs(local_delta[!truth])),
    full_auc = auc_rank(full_primary$observed$local, truth),
    reduced_auc = auc_rank(reduced$observed$local, truth),
    auc_change_when_included = auc_rank(full_primary$observed$local, truth) - auc_rank(reduced$observed$local, truth)
  )
  if (!is.null(full_primary$local) && !is.null(reduced$local)) {
    row$classification_changes <- sum(full_primary$local$classification != reduced$local$classification)
    row$true_similarity_detection_change <- sum(full_primary$local$classification[truth] == "Similarity cluster") -
      sum(reduced$local$classification[truth] == "Similarity cluster")
    row$background_similarity_false_positive_change <- sum(full_primary$local$classification[!truth] == "Similarity cluster") -
      sum(reduced$local$classification[!truth] == "Similarity cluster")
  }
  row
}) |>
  bind_rows()
utils::write.csv(primary_lovo, file.path(dirs["tables"], "primary_lovo_summary.csv"), row.names = FALSE)

write_status("stability_start", sprintf("%d paired RF seeds; ntree=%d", config$stability_reps, config$stability_ntree))
stability_rows <- list()
for (rep_id in seq_len(config$stability_reps)) {
  seed <- config$stability_seed_start + rep_id - 1L
  synthetic <- make_addcl1(data_full, seed + 100000L)
  observed_rep <- list()
  for (name in names(variants)) {
    variables <- variants[[name]]
    similarity <- fit_similarity(
      data_full[, variables, drop = FALSE],
      synthetic[, variables, drop = FALSE],
      seed,
      config$stability_ntree,
      config$rf_mtry,
      config$rf_nodesize
    )
    observed_rep[[name]] <- mlsas_observed_optimized(nb, similarity, config$style)
    rm(similarity)
  }
  full_obs <- observed_rep$full
  for (name in setdiff(names(variants), "full")) {
    reduced <- observed_rep[[name]]
    local_delta <- full_obs$local - reduced$local
    stability_rows[[length(stability_rows) + 1L]] <- data.frame(
      repetition = rep_id,
      seed = seed,
      removed_variable = sub("leave_(.*)_out", "\\1", name),
      full_global_mlsas = full_obs$global,
      reduced_global_mlsas = reduced$global,
      global_change_when_included = full_obs$global - reduced$global,
      mean_local_change_when_included = mean(local_delta),
      mean_absolute_local_change = mean(abs(local_delta)),
      true_cluster_mean_change = mean(local_delta[truth]),
      background_mean_change = mean(local_delta[!truth]),
      true_cluster_mean_absolute_change = mean(abs(local_delta[truth])),
      background_mean_absolute_change = mean(abs(local_delta[!truth])),
      full_auc = auc_rank(full_obs$local, truth),
      reduced_auc = auc_rank(reduced$local, truth),
      auc_change_when_included = auc_rank(full_obs$local, truth) - auc_rank(reduced$local, truth)
    )
  }
  utils::write.csv(bind_rows(stability_rows), file.path(dirs["tables"], "stability_by_repetition.csv"), row.names = FALSE)
  write_status("stability_progress", sprintf("%d/%d seed=%d", rep_id, config$stability_reps, seed))
  rm(synthetic, observed_rep)
  invisible(gc())
}

stability <- bind_rows(stability_rows)
stability_summary <- stability |>
  # Rename the repetition-level metric before summarise(). Otherwise dplyr's
  # sequential evaluation makes sd(mean_absolute_local_change) refer to the
  # scalar mean created on the preceding line, yielding NA.
  rename(repetition_absolute_local_change = mean_absolute_local_change) |>
  group_by(removed_variable) |>
  summarise(
    repetitions = n(),
    mean_global_change = mean(global_change_when_included),
    sd_global_change = sd(global_change_when_included),
    mean_absolute_local_change = mean(repetition_absolute_local_change),
    sd_absolute_local_change = sd(repetition_absolute_local_change),
    mean_true_cluster_change = mean(true_cluster_mean_change),
    sd_true_cluster_change = sd(true_cluster_mean_change),
    mean_auc_change = mean(auc_change_when_included),
    sd_auc_change = sd(auc_change_when_included),
    auc_change_positive_repetitions = sum(auc_change_when_included > 0),
    true_cluster_change_positive_repetitions = sum(true_cluster_mean_change > 0),
    .groups = "drop"
  ) |>
  arrange(desc(mean_absolute_local_change))
utils::write.csv(stability_summary, file.path(dirs["tables"], "stability_summary.csv"), row.names = FALSE)

plot_data <- stability |>
  mutate(removed_variable = factor(removed_variable, levels = c("x1", "x2", "cat1")))
p1 <- ggplot(plot_data, aes(removed_variable, mean_absolute_local_change, fill = removed_variable)) +
  geom_boxplot(width = 0.68, outlier.size = 0.8) +
  scale_fill_manual(values = c(x1 = "#4477AA", x2 = "#66CCEE", cat1 = "#EE7733")) +
  labs(x = "Removed variable", y = "Mean absolute local ML-SAS change", title = "DGP 4.3 paired LOVO sensitivity") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
ggsave(file.path(dirs["figures"], "lovo_local_change.png"), p1, width = 7, height = 5, dpi = 300, bg = "white")

p2 <- ggplot(plot_data, aes(removed_variable, auc_change_when_included, fill = removed_variable)) +
  geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
  geom_boxplot(width = 0.68, outlier.size = 0.8) +
  scale_fill_manual(values = c(x1 = "#4477AA", x2 = "#66CCEE", cat1 = "#EE7733")) +
  labs(x = "Removed variable", y = "Full minus reduced local-statistic AUC", title = "DGP 4.3 LOVO change in true-cluster separation") +
  theme_bw(base_size = 11) + theme(legend.position = "none")
ggsave(file.path(dirs["figures"], "lovo_auc_change.png"), p2, width = 7, height = 5, dpi = 300, bg = "white")

if (isTRUE(config$run_primary_inference)) {
  map_levels <- c("Similarity cluster", "Dissimilarity/outlier", "Not significant")
  map_colours <- c("Similarity cluster" = "#D73027", "Dissimilarity/outlier" = "#4575B4", "Not significant" = "#E6E6E6")
  map_data <- bind_rows(lapply(c("full", "leave_cat1_out"), function(name) {
    out <- dgp$sf
    out$classification <- factor(primary_objects[[name]]$local$classification, levels = map_levels)
    out$specification <- ifelse(name == "full", "Full: x1 + x2 + cat1", "LOVO: x1 + x2")
    out
  }))
  p3 <- ggplot(map_data) +
    geom_sf(aes(fill = classification), colour = "white", linewidth = 0.08) +
    facet_wrap(~ specification) +
    scale_fill_manual(values = map_colours, drop = FALSE) +
    labs(fill = NULL, title = "DGP 4.3 significance after removing cat1") +
    theme_void(base_size = 10) + theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
  ggsave(file.path(dirs["figures"], "full_vs_leave_cat1_significance.png"), p3, width = 10, height = 5.5, dpi = 300, bg = "white")
}

summary_lines <- c(
  "# DGP 4.3 paired LOVO analysis",
  "",
  sprintf("- Full variables: `%s`", paste(config$variables, collapse = "`, `")),
  sprintf("- Confirmatory RF: %d trees, mtry %d, nodesize %d, seed %d", config$rf_ntree, config$rf_mtry, config$rf_nodesize, config$primary_rf_seed),
  if (isTRUE(config$run_primary_inference)) {
    sprintf("- Confirmatory inference: %s shared label permutations, BH-adjusted local p-values", format(config$nperm, big.mark = ","))
  } else {
    "- Confirmatory inference: not run; this output is limited to paired surface sensitivity"
  },
  sprintf("- Stability: %d paired RF seeds with %d trees", config$stability_reps, config$stability_ntree),
  "- Pairing: every full/reduced comparison reuses the same RF seed and the same Addcl1 values for shared columns.",
  "",
  "## Primary LOVO results",
  "",
  paste(capture.output(print(primary_lovo, row.names = FALSE)), collapse = "\n"),
  "",
  "## Stability summary",
  "",
  paste(capture.output(print(stability_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Interpretation boundary",
  "",
  "LOVO measures conditional sensitivity of DGP 4.3 ML-SAS to removing one input. Surface-change and AUC quantities are not causal variable importance measures. Detection-change claims require a separate approved local-inference run."
)
writeLines(summary_lines, file.path(run_dir, "DGP43_LOVO_REPORT.md"))

manifest_files <- list.files(run_dir, recursive = TRUE, full.names = TRUE)
manifest_files <- manifest_files[file.info(manifest_files)$isdir %in% FALSE]
manifest <- data.frame(
  file = substring(manifest_files, nchar(run_dir) + 2L),
  bytes = file.info(manifest_files)$size,
  md5 = unname(tools::md5sum(manifest_files))
)
utils::write.csv(manifest, file.path(dirs["metadata"], "manifest.csv"), row.names = FALSE)
write_status("complete", sprintf("outputs=%d", nrow(manifest) + 1L))

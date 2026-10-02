args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 05_empirical-study/R/run_reduced_urbanity_analysis.R REPO_DIR RUN_DIR")
}

repo_dir <- normalizePath(args[1], mustWork = TRUE)
run_dir <- normalizePath(args[2], mustWork = FALSE)
if (file.exists(run_dir)) stop("Output must not exist: ", run_dir)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
study_dir <- file.path(repo_dir, "05_empirical-study")

config <- source(file.path(study_dir, "config", "reduced_urbanity_production.R"))$value
source(file.path(repo_dir, "02_revision-simulation", "R", "mlsas_helpers.R"))
source(file.path(repo_dir, "02_revision-simulation", "R", "mlsas_optimized.R"))
source(file.path(study_dir, "R", "mlsas_inference.R"))

suppressPackageStartupMessages({
  library(tmap)
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

make_addcl1 <- function(x, seed) {
  set.seed(seed)
  as.data.frame(lapply(x, function(col) sample(col, replace = TRUE)))
}

fit_similarity <- function(x, synthetic, seed, ntree, validate = FALSE) {
  set.seed(seed)
  compute_path_proximity_optimized(
    data = x,
    n_tree = ntree,
    workers = config$workers,
    # Preserve the original 48-worker tree summation groups at any worker count.
    chunk_size = max(1L, ceiling(ntree / 96L)),
    validate = validate,
    synthetic_data = synthetic
  )
}

write_status("start", config$run_label)
saveRDS(config, file.path(dirs["metadata"], "config.rds"))
dput(config, file = file.path(dirs["metadata"], "config.R"))
capture.output(sessionInfo(), file = file.path(dirs["metadata"], "sessionInfo.txt"))

source(file.path(study_dir,"R","prepare_inputs.R"))
prepared<-prepare_empirical_inputs();check_empirical_inputs(prepared,repo_dir)
main_sf<-prepared$sf;queen_nb<-prepared$queen
continuous_data<-prepared$continuous;mixed_data<-prepared$mixed
mixed_vars<-c(config$continuous_vars,config$categorical_var)
saveRDS(list(sf = main_sf, queen = queen_nb, continuous = continuous_data, mixed = mixed_data),
        file.path(dirs["objects"], "analysis_inputs.rds"))
write_status("data_prepared", sprintf("n=%d continuous=%s categorical=%s", nrow(main_sf), paste(config$continuous_vars, collapse = ";"), config$categorical_var))

primary_synthetic_mixed <- make_addcl1(mixed_data, config$primary_synthetic_seed)
primary_synthetic_continuous <- primary_synthetic_mixed[, config$continuous_vars, drop = FALSE]

write_status("primary_fit", sprintf("continuous ntree=%d", config$main_ntree))
similarity_continuous <- fit_similarity(
  continuous_data, primary_synthetic_continuous,
  config$primary_seed, config$main_ntree, validate = TRUE
)
write_status("primary_fit", sprintf("mixed ntree=%d", config$main_ntree))
similarity_mixed <- fit_similarity(
  mixed_data, primary_synthetic_mixed,
  config$primary_seed, config$main_ntree, validate = TRUE
)
saveRDS(similarity_continuous, file.path(dirs["objects"], "similarity_continuous.rds"))
saveRDS(similarity_mixed, file.path(dirs["objects"], "similarity_mixed.rds"))

write_status("permutations", sprintf("nperm=%d seed=%d", config$nperm, config$permutation_seed))
permutations <- make_permutations(nrow(main_sf), config$nperm, config$permutation_seed)
saveRDS(permutations, file.path(dirs["objects"], "permutation_indices.rds"), compress = "xz")

write_status("inference", "continuous")
inference_continuous <- mlsas_permutation_inference(
  queen_nb, similarity_continuous, permutations, config$workers,
  alpha = config$alpha, p_adjust_method = config$p_adjust_method
)
write_status("inference", "mixed")
inference_mixed <- mlsas_permutation_inference(
  queen_nb, similarity_mixed, permutations, config$workers,
  alpha = config$alpha, p_adjust_method = config$p_adjust_method
)

id_cols <- sf::st_drop_geometry(main_sf)[, c("code", "name", "province")]
local_continuous <- cbind(id_cols, specification = "Three continuous variables", inference_continuous$local)
local_mixed <- cbind(id_cols, specification = "Three continuous variables + urbanity", inference_mixed$local)
utils::write.csv(local_continuous, file.path(dirs["tables"], "local_mlsas_continuous.csv"), row.names = FALSE)
utils::write.csv(local_mixed, file.path(dirs["tables"], "local_mlsas_mixed.csv"), row.names = FALSE)

global_results <- bind_rows(
  cbind(specification = "Three continuous variables", inference_continuous$global),
  cbind(specification = "Three continuous variables + urbanity", inference_mixed$global)
)
utils::write.csv(global_results, file.path(dirs["tables"], "global_mlsas.csv"), row.names = FALSE)

comparison <- id_cols |>
  mutate(
    continuous_local_mlsas = inference_continuous$local$local_mlsas,
    mixed_local_mlsas = inference_mixed$local$local_mlsas,
    delta_local_mlsas = mixed_local_mlsas - continuous_local_mlsas,
    continuous_p_adjusted = inference_continuous$local$p_adjusted,
    mixed_p_adjusted = inference_mixed$local$p_adjusted,
    continuous_class = inference_continuous$local$classification,
    mixed_class = inference_mixed$local$classification,
    class_changed = continuous_class != mixed_class
  )
utils::write.csv(comparison, file.path(dirs["tables"], "continuous_vs_mixed_local.csv"), row.names = FALSE)

classification_summary <- bind_rows(
  data.frame(specification = "Three continuous variables", classification = inference_continuous$local$classification),
  data.frame(specification = "Three continuous variables + urbanity", classification = inference_mixed$local$classification)
) |>
  count(specification, classification, name = "municipalities")
utils::write.csv(classification_summary, file.path(dirs["tables"], "classification_summary.csv"), row.names = FALSE)

write_status("ablation", sprintf("%d variables x %d paired repetitions", length(mixed_vars), config$ablation_reps))
ablation_rows <- list()
for (rep_id in seq_len(config$ablation_reps)) {
  seed <- config$ablation_seeds[rep_id]
  paired_synthetic <- make_addcl1(mixed_data, seed + 200000L)
  full_similarity <- fit_similarity(mixed_data, paired_synthetic, seed, config$ablation_ntree)
  full_observed <- mlsas_observed(queen_nb, full_similarity)
  for (variable in mixed_vars) {
    reduced_vars <- setdiff(mixed_vars, variable)
    reduced_similarity <- fit_similarity(
      mixed_data[, reduced_vars, drop = FALSE],
      paired_synthetic[, reduced_vars, drop = FALSE],
      seed, config$ablation_ntree
    )
    reduced_observed <- mlsas_observed(queen_nb, reduced_similarity)
    local_delta <- full_observed$local - reduced_observed$local
    ablation_rows[[length(ablation_rows) + 1L]] <- data.frame(
      repetition = rep_id,
      seed = seed,
      removed_variable = variable,
      full_global_mlsas = full_observed$global,
      reduced_global_mlsas = reduced_observed$global,
      global_change_when_included = full_observed$global - reduced_observed$global,
      mean_local_change_when_included = mean(local_delta),
      mean_absolute_local_change = mean(abs(local_delta)),
      max_absolute_local_change = max(abs(local_delta))
    )
  }
  utils::write.csv(bind_rows(ablation_rows), file.path(dirs["tables"], "ablation_by_repetition.csv"), row.names = FALSE)
  write_status("ablation_progress", sprintf("%d/%d seed=%d", rep_id, config$ablation_reps, seed))
}

ablation_by_rep <- bind_rows(ablation_rows)
ablation_summary <- ablation_by_rep |>
  group_by(removed_variable) |>
  summarise(
    repetitions = n(),
    mean_global_change_when_included = mean(global_change_when_included),
    sd_global_change_when_included = sd(global_change_when_included),
    mean_abs_local_change = mean(mean_absolute_local_change),
    sd_abs_local_change = sd(mean_absolute_local_change),
    max_absolute_local_change = max(max_absolute_local_change),
    .groups = "drop"
  ) |>
  rename(
    mean_absolute_local_change = mean_abs_local_change,
    sd_absolute_local_change = sd_abs_local_change
  ) |>
  arrange(desc(mean_absolute_local_change))
utils::write.csv(ablation_summary, file.path(dirs["tables"], "ablation_summary.csv"), row.names = FALSE)

class_levels <- c("Similarity cluster", "Dissimilarity/outlier", "Not significant")
class_colours <- c(
  "Similarity cluster" = "#D73027",
  "Dissimilarity/outlier" = "#4575B4",
  "Not significant" = "#E6E6E6"
)
map_data <- bind_rows(lapply(c("continuous", "mixed"), function(spec) {
  out <- main_sf
  if (spec == "continuous") {
    out$classification <- factor(inference_continuous$local$classification, levels = class_levels)
    out$specification <- "population + dwelling value + employment rate"
  } else {
    out$classification <- factor(inference_mixed$local$classification, levels = class_levels)
    out$specification <- "Three continuous variables + ordered urbanity"
  }
  out
}))
p_map <- ggplot(map_data) +
  geom_sf(aes(fill = classification), colour = "white", linewidth = 0.06) +
  facet_wrap(~ specification) +
  scale_fill_manual(values = class_colours, drop = FALSE) +
  labs(fill = NULL, title = "Local ML-SAS significance with and without urbanity", subtitle = "BH-adjusted q < .05; queen contiguity") +
  theme_void(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"))
ggsave(file.path(dirs["figures"], "significance_with_without_urbanity.png"), p_map, width = 12, height = 6.2, dpi = 300, bg = "white")

ablation_plot <- ablation_summary |>
  mutate(
    variable_type = ifelse(removed_variable == config$categorical_var, "Ordered categorical", "Continuous"),
    removed_variable = reorder(removed_variable, mean_absolute_local_change),
    lower = pmax(0, mean_absolute_local_change - sd_absolute_local_change),
    upper = mean_absolute_local_change + sd_absolute_local_change
  )
p_lovo <- ggplot(ablation_plot, aes(removed_variable, mean_absolute_local_change)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0, colour = "grey45") +
  geom_col(aes(fill = variable_type), width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Continuous" = "#4477AA", "Ordered categorical" = "#EE7733")) +
  labs(x = NULL, y = "Mean absolute change in local ML-SAS", fill = NULL, title = "Paired LOVO in the reduced mixed specification", subtitle = sprintf("%d paired RF seeds; %d trees per fit", config$ablation_reps, config$ablation_ntree)) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(dirs["figures"], "reduced_model_lovo.png"), p_lovo, width = 8, height = 5.5, dpi = 300, bg = "white")

report <- c(
  "# Reduced-variable NLD_muni analysis with ordered urbanity",
  "",
  sprintf("- Continuous variables: `%s`", paste(config$continuous_vars, collapse = "`, `")),
  sprintf("- Ordered categorical variable: `%s`", config$categorical_var),
  sprintf("- Primary analysis: %d trees and %s shared label permutations", config$main_ntree, format(config$nperm, big.mark = ",")),
  sprintf("- LOVO: %d paired RF seeds and %d trees per fit", config$ablation_reps, config$ablation_ntree),
  "",
  "## Global inference",
  "",
  paste(capture.output(print(global_results, row.names = FALSE)), collapse = "\n"),
  "",
  "## Local classification counts",
  "",
  paste(capture.output(print(classification_summary, row.names = FALSE)), collapse = "\n"),
  "",
  sprintf("- Municipalities whose classification changed after adding urbanity: %d", sum(comparison$class_changed)),
  "",
  "## Paired LOVO",
  "",
  paste(capture.output(print(ablation_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Interpretation boundary",
  "",
  "LOVO is feature-removal sensitivity of ML-SAS, not causal importance. This reduced specification is a targeted supplementary analysis and does not replace the full 14-continuous-variable empirical model."
)
writeLines(report, file.path(run_dir, "REDUCED_URBANITY_REPORT.md"))

manifest_files <- list.files(run_dir, recursive = TRUE, full.names = TRUE)
manifest_files <- manifest_files[file.info(manifest_files)$isdir %in% FALSE]
manifest <- data.frame(
  file = substring(manifest_files, nchar(run_dir) + 2L),
  bytes = file.info(manifest_files)$size,
  md5 = unname(tools::md5sum(manifest_files))
)
utils::write.csv(manifest, file.path(dirs["metadata"], "manifest.csv"), row.names = FALSE)
write_status("complete", sprintf("outputs=%d", nrow(manifest) + 1L))

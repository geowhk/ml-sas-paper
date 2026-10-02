args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript rerun_conditional_local_inference.R REPO_DIR SOURCE_RUN_DIR OUTPUT_RUN_DIR BASE_RUN_DIR")
}

repo_dir <- normalizePath(args[1], mustWork = TRUE)
source_run_dir <- normalizePath(args[2], mustWork = TRUE)
output_run_dir <- normalizePath(args[3], mustWork = FALSE)
if (dir.exists(output_run_dir) && length(list.files(output_run_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Output run directory already exists and is not empty: ", output_run_dir)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(Matrix)
  library(sf)
  library(spdep)
})

source(file.path(repo_dir, "02_revision-simulation", "R", "mlsas_helpers.R"))
source(file.path(repo_dir, "05_empirical-study", "R", "mlsas_inference.R"))

dirs <- file.path(output_run_dir, c("inputs", "tables", "objects", "figures", "metadata", "logs"))
names(dirs) <- c("inputs", "tables", "objects", "figures", "metadata", "logs")
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
progress_file <- file.path(output_run_dir, "progress.csv")

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

nperm <- 99999L
alpha <- 0.05
p_adjust_method <- "BH"

source_config <- readRDS(file.path(source_run_dir, "metadata", "config.rds"))
base_run_dir <- normalizePath(args[4], mustWork = TRUE)
input_files <- c(
  analysis_inputs = file.path(base_run_dir, "objects", "analysis_inputs.rds"),
  similarity_continuous = file.path(base_run_dir, "objects", "similarity_continuous.rds"),
  similarity_mixed = file.path(base_run_dir, "objects", "similarity_mixed.rds"),
  permutation_indices = file.path(source_run_dir, "objects", "permutation_indices.rds"),
  old_global = file.path(source_run_dir, "tables", "global_mlsas.csv"),
  old_local_continuous = file.path(source_run_dir, "tables", "local_mlsas_continuous.csv"),
  old_local_mixed = file.path(source_run_dir, "tables", "local_mlsas_mixed.csv")
)
stopifnot(all(file.exists(input_files)))

source_manifest <- data.frame(
  role = names(input_files),
  path = unname(normalizePath(input_files)),
  bytes = unname(file.info(input_files)$size),
  md5 = unname(tools::md5sum(input_files)),
  stringsAsFactors = FALSE
)
utils::write.csv(source_manifest, file.path(dirs["inputs"], "SOURCE_MANIFEST.csv"), row.names = FALSE)

write_status("start", sprintf("conditional local inference; nperm=%d", nperm))
inputs <- readRDS(input_files["analysis_inputs"])
similarity_continuous <- readRDS(input_files["similarity_continuous"])
similarity_mixed <- readRDS(input_files["similarity_mixed"])
old_continuous <- utils::read.csv(input_files["old_local_continuous"], stringsAsFactors = FALSE)
old_mixed <- utils::read.csv(input_files["old_local_mixed"], stringsAsFactors = FALSE)
old_global <- utils::read.csv(input_files["old_global"], stringsAsFactors = FALSE)
permutations <- readRDS(input_files["permutation_indices"])

stopifnot(
  inherits(inputs$sf, "sf"),
  nrow(inputs$sf) == 337L,
  identical(dim(similarity_continuous), c(337L, 337L)),
  identical(dim(similarity_mixed), c(337L, 337L)),
  identical(dim(permutations), c(337L, nperm)),
  nrow(old_continuous) == 337L,
  nrow(old_mixed) == 337L
)

conditional_summary <- function(similarity) {
  observed <- mlsas_observed(inputs$queen, similarity, style = "W")
  n <- length(observed$local)
  mc <- matrix(NA_real_, nrow = n, ncol = 4L)
  colnames(mc) <- c("permutation_mean", "permutation_sd", "count_high", "count_low")

  for (i in seq_len(n)) {
    neighbours <- which(observed$P[i, ] != 0)
    assigned <- permutations[neighbours, , drop = FALSE]
    focal_displaced <- assigned == i
    if (any(focal_displaced)) {
      focal_labels <- matrix(
        permutations[i, ],
        nrow = length(neighbours), ncol = nperm,
        byrow = TRUE
      )
      assigned[focal_displaced] <- focal_labels[focal_displaced]
    }
    q_values <- matrix(
      observed$Q[i, assigned],
      nrow = length(neighbours), ncol = nperm
    )
    permuted_local <- n * colSums(
      q_values * observed$P[i, neighbours]
    )
    mc[i, ] <- c(
      mean(permuted_local),
      stats::sd(permuted_local),
      sum(permuted_local >= observed$local[i]),
      sum(permuted_local <= observed$local[i])
    )
  }

  analytic_mean <- analytic_variance <- numeric(n)
  for (i in seq_len(n)) {
    p_i <- observed$P[i, -i] * n
    q_i <- observed$Q[i, -i]
    analytic_mean[i] <- sum(p_i) * mean(q_i)
    analytic_variance[i] <- sum((p_i - mean(p_i))^2) *
      sum((q_i - mean(q_i))^2) / (n - 2)
  }

  p_high <- (mc[, "count_high"] + 1) / (nperm + 1)
  p_low <- (mc[, "count_low"] + 1) / (nperm + 1)
  p_two <- pmin(1, 2 * pmin(p_high, p_low))
  p_adjusted <- stats::p.adjust(p_two, method = p_adjust_method)
  direction <- ifelse(
    observed$local >= mc[, "permutation_mean"],
    "Similarity cluster",
    "Dissimilarity/outlier"
  )

  data.frame(
    local_mlsas = observed$local,
    permutation_mean = mc[, "permutation_mean"],
    permutation_sd = mc[, "permutation_sd"],
    analytic_mean = analytic_mean,
    analytic_sd = sqrt(pmax(0, analytic_variance)),
    mc_mean_se_z = (mc[, "permutation_mean"] - analytic_mean) /
      (sqrt(pmax(0, analytic_variance)) / sqrt(nperm)),
    p_high = p_high,
    p_low = p_low,
    p_two_sided = p_two,
    p_adjusted = p_adjusted,
    classification = ifelse(p_adjusted <= alpha, direction, "Not significant"),
    stringsAsFactors = FALSE
  )
}

write_status("inference", "continuous specification")
new_continuous_core <- conditional_summary(similarity_continuous)
write_status("inference", "mixed specification")
new_mixed_core <- conditional_summary(similarity_mixed)

id_cols <- sf::st_drop_geometry(inputs$sf)[, c("code", "name", "province")]
new_continuous <- cbind(
  id_cols,
  specification = "Three continuous variables",
  new_continuous_core
)
new_mixed <- cbind(
  id_cols,
  specification = "Three continuous variables + urbanity",
  new_mixed_core
)

new_comparison <- id_cols |>
  mutate(
    continuous_local_mlsas = new_continuous$local_mlsas,
    mixed_local_mlsas = new_mixed$local_mlsas,
    delta_local_mlsas = mixed_local_mlsas - continuous_local_mlsas,
    continuous_p_adjusted = new_continuous$p_adjusted,
    mixed_p_adjusted = new_mixed$p_adjusted,
    continuous_class = new_continuous$classification,
    mixed_class = new_mixed$classification,
    class_changed = continuous_class != mixed_class
  )

classification_summary <- bind_rows(
  data.frame(specification = "Three continuous variables", classification = new_continuous$classification),
  data.frame(specification = "Three continuous variables + urbanity", classification = new_mixed$classification)
) |>
  count(specification, classification, name = "municipalities")

old_new_comparison <- bind_rows(
  data.frame(
    id_cols,
    specification = "Three continuous variables",
    local_mlsas = new_continuous$local_mlsas,
    old_local_mlsas = old_continuous$local_mlsas,
    old_permutation_mean = old_continuous$permutation_mean,
    new_permutation_mean = new_continuous$permutation_mean,
    old_p_adjusted = old_continuous$p_adjusted,
    new_p_adjusted = new_continuous$p_adjusted,
    old_classification = old_continuous$classification,
    new_classification = new_continuous$classification
  ),
  data.frame(
    id_cols,
    specification = "Three continuous variables + urbanity",
    local_mlsas = new_mixed$local_mlsas,
    old_local_mlsas = old_mixed$local_mlsas,
    old_permutation_mean = old_mixed$permutation_mean,
    new_permutation_mean = new_mixed$permutation_mean,
    old_p_adjusted = old_mixed$p_adjusted,
    new_p_adjusted = new_mixed$p_adjusted,
    old_classification = old_mixed$classification,
    new_classification = new_mixed$classification
  )
) |>
  mutate(
    permutation_mean_change = new_permutation_mean - old_permutation_mean,
    adjusted_p_change = new_p_adjusted - old_p_adjusted,
    classification_changed = old_classification != new_classification
  )

old_new_summary <- old_new_comparison |>
  group_by(specification) |>
  summarise(
    municipalities = n(),
    unchanged_observed_statistic = max(abs(local_mlsas - old_local_mlsas)) < 1e-12,
    mean_old_permutation_mean = mean(old_permutation_mean),
    mean_new_permutation_mean = mean(new_permutation_mean),
    mean_absolute_permutation_mean_change = mean(abs(permutation_mean_change)),
    median_absolute_adjusted_p_change = median(abs(adjusted_p_change)),
    classification_changes = sum(classification_changed),
    .groups = "drop"
  )

utils::write.csv(old_global, file.path(dirs["tables"], "global_mlsas.csv"), row.names = FALSE)
utils::write.csv(new_continuous, file.path(dirs["tables"], "local_mlsas_continuous.csv"), row.names = FALSE)
utils::write.csv(new_mixed, file.path(dirs["tables"], "local_mlsas_mixed.csv"), row.names = FALSE)
utils::write.csv(new_comparison, file.path(dirs["tables"], "continuous_vs_mixed_local.csv"), row.names = FALSE)
utils::write.csv(classification_summary, file.path(dirs["tables"], "classification_summary.csv"), row.names = FALSE)
utils::write.csv(old_new_comparison, file.path(dirs["tables"], "total_vs_conditional_local.csv"), row.names = FALSE)
utils::write.csv(old_new_summary, file.path(dirs["tables"], "total_vs_conditional_summary.csv"), row.names = FALSE)

validation <- data.frame(
  specification = c("Three continuous variables", "Three continuous variables + urbanity"),
  max_abs_mc_mean_se_z = c(
    max(abs(new_continuous$mc_mean_se_z)),
    max(abs(new_mixed$mc_mean_se_z))
  ),
  mean_abs_mc_minus_analytic_mean = c(
    mean(abs(new_continuous$permutation_mean - new_continuous$analytic_mean)),
    mean(abs(new_mixed$permutation_mean - new_mixed$analytic_mean))
  ),
  max_relative_sd_error = c(
    max(abs(new_continuous$permutation_sd - new_continuous$analytic_sd) / new_continuous$analytic_sd),
    max(abs(new_mixed$permutation_sd - new_mixed$analytic_sd) / new_mixed$analytic_sd)
  )
)
utils::write.csv(validation, file.path(dirs["metadata"], "MONTE_CARLO_VALIDATION.csv"), row.names = FALSE)
if (any(validation$max_abs_mc_mean_se_z > 6.5) || any(validation$max_relative_sd_error > 0.03)) {
  stop("Monte Carlo results failed the prespecified analytical-moment validation thresholds.")
}

config <- list(
  run_label = basename(output_run_dir),
  source_total_run_dir = source_run_dir,
  source_base_run_dir = base_run_dir,
  reused_similarity_matrices = TRUE,
  global_inference = "unchanged total-label randomization from source run",
  local_inference = "focal-profile conditional randomization",
  local_sampling = "each stored total permutation transformed by swapping the focal profile into the focal position; the remaining n-1 profiles are thereby uniformly permuted",
  tails = "two-sided; equality counted in both one-sided tails",
  p_value = "plus-one Monte Carlo",
  multiplicity = "BH across 337 municipalities separately by specification",
  nperm = nperm,
  source_permutation_seed = source_config$permutation_seed,
  shared_transformed_draws_across_specifications = TRUE,
  alpha = alpha,
  p_adjust_method = p_adjust_method
)
saveRDS(config, file.path(dirs["metadata"], "config.rds"))
dput(config, file = file.path(dirs["metadata"], "config.R"))
capture.output(sessionInfo(), file = file.path(dirs["metadata"], "sessionInfo.txt"))

transition_table <- old_new_comparison |>
  count(specification, old_classification, new_classification, name = "municipalities")
utils::write.csv(transition_table, file.path(dirs["tables"], "total_vs_conditional_transitions.csv"), row.names = FALSE)

report <- c(
  "# Conditional local empirical inference rerun",
  "",
  sprintf("- Immutable source total-randomization run: `%s`", source_run_dir),
  sprintf("- Immutable source forest/input run: `%s`", base_run_dir),
  "- Observed 1,000-tree path-proximity matrices and global inference were held fixed.",
  sprintf("- New local inference: focal-profile conditional randomization, %s draws per municipality.", format(nperm, big.mark = ",")),
  "- Local tests: two-sided plus-one Monte Carlo p-values; BH adjustment across 337 municipalities within each specification.",
  "",
  "## New local classification counts",
  "",
  paste(capture.output(print(classification_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Total versus conditional summary",
  "",
  paste(capture.output(print(old_new_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Classification transitions from old total to new conditional inference",
  "",
  paste(capture.output(print(transition_table, row.names = FALSE)), collapse = "\n"),
  "",
  "## Monte Carlo validation against analytical conditional moments",
  "",
  paste(capture.output(print(validation, row.names = FALSE)), collapse = "\n")
)
writeLines(report, file.path(output_run_dir, "COMPARISON_REPORT.md"))

manifest_files <- list.files(output_run_dir, recursive = TRUE, full.names = TRUE)
manifest_files <- manifest_files[file.info(manifest_files)$isdir %in% FALSE]
manifest_path <- file.path(dirs["metadata"], "manifest.csv")
manifest_files <- setdiff(manifest_files, manifest_path)
manifest <- data.frame(
  path = substring(manifest_files, nchar(output_run_dir) + 2L),
  bytes = file.info(manifest_files)$size,
  md5 = unname(tools::md5sum(manifest_files))
)
utils::write.csv(manifest, manifest_path, row.names = FALSE)
write_status("complete", sprintf("outputs=%d", nrow(manifest) + 1L))

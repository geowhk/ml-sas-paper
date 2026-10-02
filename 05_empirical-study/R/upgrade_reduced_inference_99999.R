args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript upgrade_reduced_inference_99999.R REPO_DIR SOURCE_RUN_DIR OUTPUT_RUN_DIR")
}

repo_dir <- normalizePath(args[1], mustWork = TRUE)
source_run_dir <- normalizePath(args[2], mustWork = TRUE)
output_run_dir <- normalizePath(args[3], mustWork = FALSE)
if (file.exists(output_run_dir)) stop("Output must not exist: ", output_run_dir)
dir.create(output_run_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(future)
  library(future.apply)
  library(Matrix)
  library(randomForest)
  library(sf)
  library(spdep)
})

source(file.path(repo_dir, "02_revision-simulation", "R", "mlsas_helpers.R"))
source(file.path(repo_dir, "05_empirical-study", "R", "mlsas_inference.R"))

dirs <- file.path(output_run_dir, c("tables", "objects", "metadata", "logs"))
names(dirs) <- c("tables", "objects", "metadata", "logs")
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
permutation_seed <- 20260729L
detected_cores <- parallel::detectCores(logical = TRUE)
workers <- as.integer(Sys.getenv("REDUCED_WORKERS", "1"))
stopifnot(is.finite(workers),workers>=1L)

input_file <- file.path(source_run_dir, "objects", "analysis_inputs.rds")
continuous_similarity_file <- file.path(source_run_dir, "objects", "similarity_continuous.rds")
mixed_similarity_file <- file.path(source_run_dir, "objects", "similarity_mixed.rds")

write_status("start", sprintf("nperm=%d workers=%d", nperm, workers))
inputs <- readRDS(input_file)
similarity_continuous <- readRDS(continuous_similarity_file)
similarity_mixed <- readRDS(mixed_similarity_file)

stopifnot(
  inherits(inputs$sf, "sf"),
  nrow(inputs$sf) == 337L,
  nrow(similarity_continuous) == 337L,
  nrow(similarity_mixed) == 337L
)

write_status("permutations", sprintf("generating %d permutations with seed %d", nperm, permutation_seed))
permutations <- make_permutations(nrow(inputs$sf), nperm, permutation_seed)
saveRDS(permutations, file.path(dirs["objects"], "permutation_indices.rds"), compress = "xz")

write_status("inference", "continuous")
inference_continuous <- mlsas_permutation_inference(
  inputs$queen, similarity_continuous, permutations,
  workers = workers, alpha = 0.05, p_adjust_method = "BH"
)

write_status("inference", "mixed")
inference_mixed <- mlsas_permutation_inference(
  inputs$queen, similarity_mixed, permutations,
  workers = workers, alpha = 0.05, p_adjust_method = "BH"
)

id_cols <- sf::st_drop_geometry(inputs$sf)[, c("code", "name", "province")]
local_continuous <- cbind(
  id_cols,
  specification = "Three continuous variables",
  inference_continuous$local
)
local_mixed <- cbind(
  id_cols,
  specification = "Three continuous variables + urbanity",
  inference_mixed$local
)
global_results <- bind_rows(
  cbind(specification = "Three continuous variables", inference_continuous$global),
  cbind(specification = "Three continuous variables + urbanity", inference_mixed$global)
)
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
classification_summary <- bind_rows(
  data.frame(
    specification = "Three continuous variables",
    classification = inference_continuous$local$classification
  ),
  data.frame(
    specification = "Three continuous variables + urbanity",
    classification = inference_mixed$local$classification
  )
) |>
  count(specification, classification, name = "municipalities")

utils::write.csv(global_results, file.path(dirs["tables"], "global_mlsas.csv"), row.names = FALSE)
utils::write.csv(local_continuous, file.path(dirs["tables"], "local_mlsas_continuous.csv"), row.names = FALSE)
utils::write.csv(local_mixed, file.path(dirs["tables"], "local_mlsas_mixed.csv"), row.names = FALSE)
utils::write.csv(comparison, file.path(dirs["tables"], "continuous_vs_mixed_local.csv"), row.names = FALSE)
utils::write.csv(classification_summary, file.path(dirs["tables"], "classification_summary.csv"), row.names = FALSE)

config <- list(
  run_label = "NLD_muni_reduced_urbanity_perm99999",
  source_run_dir = source_run_dir,
  reused_similarity_matrices = TRUE,
  continuous_vars = c("population", "dwelling_value", "employment_rate"),
  categorical_var = "urbanity",
  rf_seed = 20260714L,
  synthetic_seed = 20260814L,
  permutation_seed = permutation_seed,
  nperm = nperm,
  workers = workers,
  alpha = 0.05,
  p_adjust_method = "BH"
)
saveRDS(config, file.path(dirs["metadata"], "config.rds"))
dput(config, file = file.path(dirs["metadata"], "config.R"))
capture.output(sessionInfo(), file = file.path(dirs["metadata"], "sessionInfo.txt"))

report <- c(
  "# Reduced-variable NLD_muni 99,999-permutation inference",
  "",
  sprintf("- Source similarity matrices: `%s`", source_run_dir),
  "- The fitted 1,000-tree RF similarity matrices were held fixed.",
  sprintf("- Shared label permutations: %s", format(nperm, big.mark = ",")),
  "- Inference: two-sided finite-sample permutation p-values; BH adjustment across 337 municipalities.",
  "",
  "## Global inference",
  "",
  paste(capture.output(print(global_results, row.names = FALSE)), collapse = "\n"),
  "",
  "## Local classification counts",
  "",
  paste(capture.output(print(classification_summary, row.names = FALSE)), collapse = "\n"),
  "",
  sprintf("- Municipalities whose classification changed after adding urbanity: %d", sum(comparison$class_changed))
)
writeLines(report, file.path(output_run_dir, "PERM99999_REPORT.md"))

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

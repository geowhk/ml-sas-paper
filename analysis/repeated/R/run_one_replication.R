#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
this_file <- normalizePath(sub("^--file=", "", file_arg[[1]]))
source(file.path(dirname(this_file), "production_common.R"))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) != 1L) stop("Expected --", name, "=...", call. = FALSE)
  sub(paste0("^--", name, "="), "", hit[[1]])
}

dgp <- as.integer(get_arg("dgp"))
replication_id <- as.integer(get_arg("replication"))
run_dir <- normalizePath(get_arg("run_dir"), mustWork = TRUE)
paths <- resolve_production_paths(this_file)
spec <- load_production_sources(paths)

max_rep <- switch(as.character(dgp), `2` = 100L, `3` = 100L, `4` = 100L, stop("Unsupported DGP", call. = FALSE))
if (replication_id < 1L || replication_id > max_rep) stop("Replication outside frozen range", call. = FALSE)

dgp_dir <- file.path(run_dir, paste0("dgp", dgp))
final_dir <- file.path(dgp_dir, sprintf("rep_%04d", replication_id))
temp_dir <- paste0(final_dir, ".partial")
if (dir.exists(final_dir) || dir.exists(temp_dir)) stop("Replication output already exists: ", final_dir, call. = FALSE)
dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

failure_file <- file.path(temp_dir, "FAILED.json")
tryCatch({
  dataset_seed <- professor_dataset_seed(dgp, replication_id)
  sf_obj <- generate_professor_dgp(dgp, replication_id)
  if (nrow(sf_obj) != 400L) stop("Generated dataset does not contain 400 cells", call. = FALSE)
  if (dgp > 2L && sum(sf_obj$true_cluster) != 60L) stop("Target mask does not contain 60 cells", call. = FALSE)

  results <- list()
  if (dgp == 2L) {
    results$path_continuous <- run_mlsas_method(sf_obj, dgp, replication_id, "path_continuous", c("x1", "x2", "x3"), 1L, "path", spec)
    results$path_mixed <- run_mlsas_method(sf_obj, dgp, replication_id, "path_mixed", c("x1", "x2", "x3", "cat1"), 2L, "path", spec)
  } else if (dgp == 3L) {
    results$local_geary <- run_local_geary_method(sf_obj, dgp, replication_id, "local_geary", c("x1", "x2", "x3"), 3L, spec)
    results$path_continuous <- run_mlsas_method(sf_obj, dgp, replication_id, "path_continuous", c("x1", "x2", "x3"), 4L, "path", spec)
  } else if (dgp == 4L) {
    results$local_geary <- run_local_geary_method(sf_obj, dgp, replication_id, "local_geary", c("x1", "x2"), 3L, spec)
    results$path_continuous <- run_mlsas_method(sf_obj, dgp, replication_id, "path_continuous", c("x1", "x2"), 4L, "path", spec)
    results$path_mixed <- run_mlsas_method(sf_obj, dgp, replication_id, "path_mixed", c("x1", "x2", "cat1"), 5L, "path", spec)
    results$terminal_mixed <- run_mlsas_method(sf_obj, dgp, replication_id, "terminal_mixed", c("x1", "x2", "cat1"), 6L, "terminal", spec)
  }

  method_names <- names(results)
  if (!identical(sort(method_names), sort(expected_methods(dgp)))) stop("Method set mismatch", call. = FALSE)
  for (method in method_names) {
    utils::write.csv(results[[method]]$local, file.path(temp_dir, paste0(method, "_local.csv")), row.names = FALSE)
  }
  metrics <- do.call(rbind, lapply(results, `[[`, "metrics"))
  utils::write.csv(metrics, file.path(temp_dir, "metrics.csv"), row.names = FALSE, na = "")

  metadata <- list(
    status = "COMPLETED",
    freeze_id = spec$freeze_id,
    dgp_source = "professor_original_DGP_generation.R",
    dgp_source_sha256 = spec$source$sha256,
    dgp = dgp,
    replication_id = replication_id,
    dataset_seed = dataset_seed,
    methods = method_names,
    started_and_completed_in_single_attempt = TRUE,
    completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    settings = spec$inference
  )
  jsonlite::write_json(metadata, file.path(temp_dir, "metadata.json"), pretty = TRUE, auto_unbox = TRUE)
  writeLines("COMPLETED", file.path(temp_dir, "COMPLETE"))
  if (!file.rename(temp_dir, final_dir)) stop("Atomic finalization failed", call. = FALSE)
  cat("COMPLETED dgp=", dgp, " replication=", replication_id, " permutations=99999\n", sep = "")
}, error = function(e) {
  failure <- list(
    status = "FAILED",
    freeze_id = spec$freeze_id,
    dgp = dgp,
    replication_id = replication_id,
    failed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    message = conditionMessage(e)
  )
  jsonlite::write_json(failure, failure_file, pretty = TRUE, auto_unbox = TRUE)
  message("FAILED dgp=", dgp, " replication=", replication_id, ": ", conditionMessage(e))
  quit(status = 1L, save = "no")
})

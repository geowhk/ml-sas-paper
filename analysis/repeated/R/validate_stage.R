#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(jsonlite))

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) != 1L) stop("Expected --", name, "=...", call. = FALSE)
  sub(paste0("^--", name, "="), "", hit[[1]])
}

dgp <- as.integer(get_arg("dgp"))
run_dir <- normalizePath(get_arg("run_dir"), mustWork = TRUE)
expected_reps <- switch(as.character(dgp), `2` = 100L, `3` = 100L, `4` = 100L, stop("Unsupported DGP", call. = FALSE))
expected_methods <- switch(as.character(dgp), `2` = 2L, `3` = 2L, `4` = 4L)
freeze_id <- "S4-REPEATED-PROFESSOR-DGP-v2"
source_hash <- "8599bf88fca254343c62486b5840597b8c5a5933e71c0ffeeec484b83fe0e90c"

dgp_dir <- file.path(run_dir, paste0("dgp", dgp))
partials <- list.files(dgp_dir, pattern = "\\.partial$", full.names = TRUE)
if (length(partials)) stop("Partial/failed replication directories present: ", paste(basename(partials), collapse = ", "), call. = FALSE)
rep_dirs <- sort(list.files(dgp_dir, pattern = "^rep_[0-9]{4}$", full.names = TRUE))
if (length(rep_dirs) != expected_reps) stop("Expected ", expected_reps, " completed replications; found ", length(rep_dirs), call. = FALSE)

all_metrics <- list()
for (path in rep_dirs) {
  if (!file.exists(file.path(path, "COMPLETE"))) stop("Missing COMPLETE: ", path, call. = FALSE)
  metadata <- jsonlite::read_json(file.path(path, "metadata.json"), simplifyVector = TRUE)
  if (!identical(metadata$status, "COMPLETED")) stop("Metadata status mismatch: ", path, call. = FALSE)
  if (!identical(metadata$freeze_id, freeze_id)) stop("Freeze ID mismatch: ", path, call. = FALSE)
  if (!identical(metadata$dgp_source_sha256, source_hash)) stop("DGP source hash mismatch: ", path, call. = FALSE)
  if (as.integer(metadata$settings$permutations) != 99999L) stop("Permutation count mismatch: ", path, call. = FALSE)
  if (!identical(metadata$settings$mlsas_alternative, "greater")) stop("ML-SAS tail mismatch: ", path, call. = FALSE)
  if (!identical(metadata$settings$p_adjust, "BH")) stop("Multiplicity rule mismatch: ", path, call. = FALSE)

  metrics <- read.csv(file.path(path, "metrics.csv"), stringsAsFactors = FALSE)
  if (nrow(metrics) != expected_methods) stop("Method count mismatch: ", path, call. = FALSE)
  local_files <- list.files(path, pattern = "_local\\.csv$", full.names = TRUE)
  if (length(local_files) != expected_methods) stop("Local-file count mismatch: ", path, call. = FALSE)
  for (file in local_files) {
    local <- read.csv(file, stringsAsFactors = FALSE)
    if (nrow(local) != 400L) stop("Local row-count mismatch: ", file, call. = FALSE)
    if (grepl("path_|terminal_", basename(file)) && "p_upper" %in% names(local)) {
      if (any(local$p_upper < 1 / 100000 - 1e-15) || any(local$p_upper > 1)) {
        stop("ML-SAS plus-one p-value grid mismatch: ", file, call. = FALSE)
      }
    }
  }
  all_metrics[[length(all_metrics) + 1L]] <- metrics
}

metrics <- do.call(rbind, all_metrics)
if (nrow(metrics) != expected_reps * expected_methods) stop("Stage metric row count mismatch", call. = FALSE)
if (anyDuplicated(metrics[c("dgp", "replication_id", "method")])) stop("Duplicate method/replication rows", call. = FALSE)
utils::write.csv(metrics, file.path(dgp_dir, paste0("dgp", dgp, "_stage_metrics.csv")), row.names = FALSE, na = "")
writeLines(
  c(
    paste0("dgp=", dgp),
    paste0("replications=", expected_reps),
    paste0("methods_per_replication=", expected_methods),
    paste0("metric_rows=", nrow(metrics)),
    "permutations=99999",
    paste0("freeze_id=", freeze_id),
    paste0("dgp_source_sha256=", source_hash),
    "status=VALIDATED"
  ),
  file.path(dgp_dir, "STAGE_VALIDATED.txt")
)
cat("VALIDATED dgp=", dgp, " replications=", expected_reps, " metric_rows=", nrow(metrics), " permutations=99999\n", sep = "")

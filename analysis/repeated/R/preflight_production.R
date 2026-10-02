#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
this_file <- normalizePath(sub("^--file=", "", file_arg[[1]]))
source(file.path(dirname(this_file), "production_common.R"))

paths <- resolve_production_paths(this_file)
spec <- load_production_sources(paths)
generator_output <- system2(
  "Rscript",
  shQuote(file.path(paths$protocol_dir, "R", "validate_professor_generators.R")),
  stdout = TRUE,
  stderr = TRUE
)
generator_status <- attr(generator_output, "status")
if (is.null(generator_status)) generator_status <- 0L
if (generator_status != 0L) stop(paste(generator_output, collapse = " | "), call. = FALSE)

stopifnot(
  identical(spec$freeze_id, "S4-REPEATED-PROFESSOR-DGP-v2"),
  identical(spec$status, "validated_ready_for_production"),
  as.integer(spec$inference$rf_ntree) == 100L,
  as.integer(spec$inference$rf_mtry) == 1L,
  as.integer(spec$inference$rf_nodesize) == 1L,
  as.integer(spec$inference$permutations) == 99999L,
  as.integer(spec$dgp2$replications) == 100L,
  as.integer(spec$dgp3$replications) == 100L,
  as.integer(spec$dgp4$replications) == 100L,
  as.integer(spec$production_gate$expected_metric_rows) == 800L,
  identical(spec$inference$mlsas_alternative, "greater"),
  identical(spec$inference$p_adjust, "BH"),
  is.function(run_mlsas_method),
  is.function(run_local_geary_method),
  is.function(generate_professor_dgp)
)

cat(paste(generator_output, collapse = "\n"), "\n")
cat("PASS: production dependencies loaded\n")
cat("PASS: all repeated DGPs use 100 datasets and the expected metric total is 800\n")
cat("PASS: frozen analysis settings include 99,999 permutations\n")

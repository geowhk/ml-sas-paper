#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
this_file <- normalizePath(sub("^--file=", "", file_arg[[1]]))
protocol_dir <- normalizePath(file.path(dirname(this_file), ".."))
repo_dir <- normalizePath(file.path(protocol_dir, "../.."))
source_snapshot <- file.path(protocol_dir, "source_snapshot", "DGP_generation.R")
expected_source_hash <- "8599bf88fca254343c62486b5840597b8c5a5933e71c0ffeeec484b83fe0e90c"

hash <- unname(tools::md5sum(source_snapshot))
sha_output <- system2("shasum", c("-a", "256", shQuote(source_snapshot)), stdout = TRUE)
actual_source_hash <- strsplit(sha_output[[1]], "[[:space:]]+")[[1]][[1]]
if (!identical(actual_source_hash, expected_source_hash)) {
  stop("Professor source SHA-256 mismatch", call. = FALSE)
}

source(file.path(protocol_dir, "R", "generate_professor_dgps.R"))

# Load the professor's unmodified script in an isolated environment. The script
# also creates the four fixed datasets; no files are written by sourcing it.
professor_env <- new.env(parent = globalenv())
sys.source(source_snapshot, envir = professor_env)

grid_new <- make_professor_grid_data(20, 20, "hex", TRUE, 1)
grid_original <- professor_env$make_grid_data(20, 20, "hex", TRUE, 1)

new_objects <- list(
  dgp1 = simulate_professor_dgp1(grid_new, seed = 2026L),
  dgp2 = simulate_professor_dgp2(grid_new, seed = 2027L),
  dgp3 = simulate_professor_dgp3(grid_new, seed = 2028L),
  dgp4 = simulate_professor_dgp4(grid_new, seed = 2029L)
)
original_objects <- list(
  dgp1 = professor_env$simulate_dgp1_univariate_geary_clusters(grid_original, seed = 2026L),
  dgp2 = professor_env$simulate_dgp2_null_mixed(grid_original, seed = 2027L),
  dgp3 = professor_env$simulate_dgp3_linear_multivariate_continuous(grid_original, seed = 2028L),
  dgp4 = professor_env$simulate_dgp4_nonlinear_mixed_interaction(grid_original, seed = 2029L)
)

legacy_files <- c(
  dgp1 = "sim_dgp1.rds",
  dgp2 = "sim_dgp2.rds",
  dgp3 = "sim_dgp3.1.rds",
  dgp4 = "sim_dgp4.3.rds"
)

assert_close_objects <- function(actual, expected, label, tolerance = 5e-15) {
  common <- intersect(names(actual), names(expected))
  for (name in setdiff(common, "geometry")) {
    a <- actual[[name]]
    b <- expected[[name]]
    if (is.numeric(a) && is.numeric(b)) {
      if (length(a) != length(b) || max(abs(a - b)) > tolerance) {
        stop(label, " numeric mismatch: ", name, call. = FALSE)
      }
    } else if (!identical(a, b)) {
      stop(label, " exact mismatch: ", name, call. = FALSE)
    }
  }
  geometry_match <- sf::st_equals_exact(actual, expected, par = 1e-12)
  if (!all(lengths(geometry_match) == 1L)) {
    stop(label, " geometry mismatch", call. = FALSE)
  }
  invisible(TRUE)
}

for (name in names(new_objects)) {
  assert_close_objects(new_objects[[name]], original_objects[[name]], paste(name, "extracted-vs-source"))
  legacy <- readRDS(file.path(repo_dir, "02_revision-simulation", "sim_dgp", legacy_files[[name]]))
  assert_close_objects(new_objects[[name]], legacy, paste(name, "extracted-vs-fixed-RDS"))
}

if (!identical(as.character(new_objects$dgp4$cat1), as.character(original_objects$dgp4$cat1))) {
  stop("DGP 4 categorical RNG order does not match professor source", call. = FALSE)
}
if (sum(new_objects$dgp2$true_cluster) != 0L) stop("DGP 2 truth mismatch", call. = FALSE)
if (sum(new_objects$dgp3$true_cluster) != 60L) stop("DGP 3 target count mismatch", call. = FALSE)
if (sum(new_objects$dgp4$true_cluster) != 60L) stop("DGP 4 target count mismatch", call. = FALSE)

same_seed_a <- generate_professor_dgp(4L, 1L)
same_seed_b <- generate_professor_dgp(4L, 1L)
different_seed <- generate_professor_dgp(4L, 2L)
if (!identical(same_seed_a, same_seed_b)) stop("Same-seed determinism failed", call. = FALSE)
if (identical(same_seed_a$x1, different_seed$x1)) stop("Different-seed variation failed", call. = FALSE)

cat("PASS: professor source SHA-256\n")
cat("PASS: extracted DGP 1-4 functions match professor source\n")
cat("PASS: fixed DGP 1-4 RDS reconstruction, including DGP 4 cat1\n")
cat("PASS: target counts and seed determinism\n")

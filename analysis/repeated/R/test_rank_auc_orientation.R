#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
this_file <- normalizePath(sub("^--file=", "", file_arg[[1]]))
source(file.path(dirname(this_file), "production_common.R"))

truth <- c(TRUE, TRUE, FALSE, FALSE)
mlsas_statistic <- c(4, 3, 2, 1)
local_geary_statistic <- c(1, 2, 3, 4)

stopifnot(
  identical(orient_rank_auc_score(mlsas_statistic, "higher"), mlsas_statistic),
  identical(orient_rank_auc_score(local_geary_statistic, "lower"), -local_geary_statistic),
  identical(safe_rank_auc(orient_rank_auc_score(mlsas_statistic, "higher"), truth), 1),
  identical(safe_rank_auc(orient_rank_auc_score(local_geary_statistic, "lower"), truth), 1)
)

cat("PASS: ML-SAS and Local Geary rank-AUC orientations are higher-is-target-like.\n")

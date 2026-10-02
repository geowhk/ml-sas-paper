#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
hit <- grep("^--run_dir=", args, value = TRUE)
if (length(hit) != 1L) stop("Expected --run_dir=...", call. = FALSE)
run_dir <- normalizePath(sub("^--run_dir=", "", hit[[1]]), mustWork = TRUE)
summary_hit <- grep("^--summary_dir=", args, value = TRUE)
if (length(summary_hit) > 1L) stop("Expected at most one --summary_dir=...", call. = FALSE)
summary_dir <- if (length(summary_hit) == 1L) {
  sub("^--summary_dir=", "", summary_hit[[1]])
} else {
  file.path(run_dir, "summary")
}
if (!grepl("^/", summary_dir)) summary_dir <- file.path(run_dir, summary_dir)
correction_hit <- grep("^--correct_legacy_local_geary_auc=", args, value = TRUE)
if (length(correction_hit) > 1L) {
  stop("Expected at most one --correct_legacy_local_geary_auc=...", call. = FALSE)
}
if (length(correction_hit)) stop("Legacy AUC correction is forbidden for v6 outputs")
correct_legacy_local_geary_auc <- FALSE
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

for (dgp in 2:4) {
  marker <- file.path(run_dir, paste0("dgp", dgp), "STAGE_VALIDATED.txt")
  if (!file.exists(marker)) stop("Stage validation marker missing for DGP ", dgp, call. = FALSE)
}

metrics <- dplyr::bind_rows(lapply(2:4, function(dgp) {
  read.csv(
    file.path(run_dir, paste0("dgp", dgp), paste0("dgp", dgp, "_stage_metrics.csv")),
    stringsAsFactors = FALSE
  )
}))
expected_rows <- 100L * (2L + 2L + 4L)
if (nrow(metrics) != expected_rows) {
  stop("Expected ", expected_rows, " method-replication rows; found ", nrow(metrics), call. = FALSE)
}
if (anyDuplicated(metrics[c("dgp", "replication_id", "method")])) {
  stop("Duplicate method-replication keys", call. = FALSE)
}
for (dgp in 2:4) {
  if (length(unique(metrics$replication_id[metrics$dgp == dgp])) != 100L) {
    stop("DGP ", dgp, " does not contain exactly 100 datasets", call. = FALSE)
  }
}
write.csv(metrics, file.path(summary_dir, "all_replication_metrics.csv"), row.names = FALSE, na = "")

exact_binomial <- function(successes, trials) {
  ci <- stats::binom.test(successes, trials, conf.level = 0.95)$conf.int
  c(estimate = successes / trials, lower = ci[[1]], upper = ci[[2]])
}

null_summary <- metrics |>
  dplyr::filter(dgp == 2) |>
  dplyr::group_by(method) |>
  dplyr::group_modify(function(data, key) {
    interval <- exact_binomial(sum(data$any_rejection), nrow(data))
    data.frame(
      replications = nrow(data),
      datasets_with_any_rejection = sum(data$any_rejection),
      prob_any_rejection = interval[["estimate"]],
      prob_any_rejection_ci_lower = interval[["lower"]],
      prob_any_rejection_ci_upper = interval[["upper"]],
      empirical_fdr_fwer = interval[["estimate"]],
      mean_rejections = mean(data$discoveries),
      sd_rejections = stats::sd(data$discoveries),
      mean_cell_false_positive_proportion = mean(data$discoveries / 400)
    )
  }) |>
  dplyr::ungroup()
write.csv(null_summary, file.path(summary_dir, "dgp2_null_operating_characteristics.csv"), row.names = FALSE)

nonnull <- metrics |> dplyr::filter(dgp %in% c(3, 4))
nonnull_summary <- nonnull |>
  dplyr::group_by(dgp, method) |>
  dplyr::summarise(
    replications = dplyr::n(),
    mean_recall = mean(recall),
    sd_recall = stats::sd(recall),
    mean_background_fpr = mean(background_fpr),
    sd_background_fpr = stats::sd(background_fpr),
    precision_defined_reps = sum(!is.na(precision)),
    conditional_mean_precision = mean(precision, na.rm = TRUE),
    conditional_sd_precision = stats::sd(precision, na.rm = TRUE),
    no_discovery_reps = sum(no_discovery),
    no_discovery_frequency = mean(no_discovery),
    mean_f1 = mean(f1),
    f1_defined_reps = sum(is.finite(f1)),
    mean_jaccard = mean(jaccard, na.rm = TRUE),
    mean_discoveries = mean(discoveries),
    mean_rank_auc = mean(rank_auc, na.rm = TRUE),
    .groups = "drop"
  )
no_disc_intervals <- nonnull |>
  dplyr::group_by(dgp, method) |>
  dplyr::group_modify(function(data, key) {
    interval <- exact_binomial(sum(data$no_discovery), nrow(data))
    data.frame(
      no_discovery_ci_lower = interval[["lower"]],
      no_discovery_ci_upper = interval[["upper"]]
    )
  }) |>
  dplyr::ungroup()
nonnull_summary <- dplyr::left_join(nonnull_summary, no_disc_intervals, by = c("dgp", "method"))
write.csv(nonnull_summary, file.path(summary_dir, "dgp3_dgp4_method_summaries.csv"), row.names = FALSE)

contrast_definitions <- data.frame(
  dgp = c(3L, 4L, 4L, 4L),
  lhs = c("path_continuous", "path_mixed", "path_mixed", "path_continuous"),
  rhs = c("local_geary", "path_continuous", "terminal_mixed", "local_geary"),
  contrast = c(
    "dgp3_path_continuous_minus_local_geary",
    "dgp4_path_mixed_minus_path_continuous",
    "dgp4_path_mixed_minus_terminal_mixed",
    "dgp4_path_continuous_minus_local_geary"
  ),
  stringsAsFactors = FALSE
)
contrast_metrics <- c("recall", "background_fpr", "precision", "f1", "jaccard", "discoveries", "rank_auc")

set.seed(5200001L)
contrast_rows <- list()
contrast_rep_rows <- list()
for (i in seq_len(nrow(contrast_definitions))) {
  def <- contrast_definitions[i, ]
  wide <- metrics |>
    dplyr::filter(dgp == def$dgp, method %in% c(def$lhs, def$rhs)) |>
    dplyr::select(replication_id, method, dplyr::all_of(contrast_metrics)) |>
    tidyr::pivot_wider(names_from = method, values_from = dplyr::all_of(contrast_metrics), names_sep = "__") |>
    dplyr::arrange(replication_id)
  if (nrow(wide) != 100L) stop("Paired contrast is incomplete: ", def$contrast, call. = FALSE)

  rep_output <- data.frame(replication_id = wide$replication_id)
  for (metric in contrast_metrics) {
    difference <- wide[[paste0(metric, "__", def$lhs)]] - wide[[paste0(metric, "__", def$rhs)]]
    rep_output[[metric]] <- difference
    values <- difference[is.finite(difference)]
    if (!length(values)) {
      estimate <- lower <- upper <- NA_real_
      n_pairs <- 0L
    } else {
      estimate <- mean(values)
      before <- .Random.seed
      boot_means <- replicate(9999L, mean(sample(values, length(values), replace = TRUE)))
      if (metric == "f1") {
        # v6 used the old finite-F1 stream for subsequent metrics only.
        legacy_f1 <- function(method) {
          z <- metrics[metrics$dgp == def$dgp & metrics$method == method, ]
          z <- z[order(z$replication_id), ]
          ifelse(is.finite(z$precision) & is.finite(z$recall) &
                 z$precision + z$recall > 0,
                 2*z$precision*z$recall/(z$precision+z$recall), NA_real_)
        }
        legacy <- legacy_f1(def$lhs) - legacy_f1(def$rhs)
        legacy <- legacy[is.finite(legacy)]
        .Random.seed <- before
        if (length(legacy)) invisible(replicate(9999L, mean(sample(legacy, length(legacy), replace=TRUE))))
      }
      interval <- stats::quantile(boot_means, probs = c(0.025, 0.975), names = FALSE, type = 7)
      lower <- interval[[1]]
      upper <- interval[[2]]
      n_pairs <- length(values)
    }
    contrast_rows[[length(contrast_rows) + 1L]] <- data.frame(
      dgp = def$dgp,
      contrast = def$contrast,
      lhs = def$lhs,
      rhs = def$rhs,
      metric = metric,
      n_paired_datasets = n_pairs,
      mean_difference = estimate,
      bootstrap_ci_lower = lower,
      bootstrap_ci_upper = upper,
      stringsAsFactors = FALSE
    )
  }
  rep_output$dgp <- def$dgp
  rep_output$contrast <- def$contrast
  contrast_rep_rows[[length(contrast_rep_rows) + 1L]] <- rep_output
}

contrast_summary <- dplyr::bind_rows(contrast_rows)
contrast_replications <- dplyr::bind_rows(contrast_rep_rows)
write.csv(contrast_summary, file.path(summary_dir, "prespecified_paired_contrasts.csv"), row.names = FALSE, na = "")
write.csv(contrast_replications, file.path(summary_dir, "paired_differences_by_replication.csv"), row.names = FALSE, na = "")

manifest <- list(
  status = "SUMMARIZED",
  freeze_id = "S4-REPEATED-PROFESSOR-DGP-v2",
  dgp_source = "professor_original_DGP_generation.R",
  dgp_source_sha256 = "8599bf88fca254343c62486b5840597b8c5a5933e71c0ffeeec484b83fe0e90c",
  summarized_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  production_rows = nrow(metrics),
  dgp2_replications = 100L,
  dgp3_replications = 100L,
  dgp4_replications = 100L,
  permutations_per_dataset_method = 99999L,
  bootstrap_resamples = 9999L,
  bootstrap_seed = 5200001L,
  rank_auc_orientation = "higher_values_indicate_target_evidence",
  local_geary_rank_auc_correction = correct_legacy_local_geary_auc,
  correction_scope = if (correct_legacy_local_geary_auc) {
    "derived_rank_auc_and_dependent_summaries_only_no_permutation_rerun"
  } else {
    "none"
  },
  source_excludes_smoke_and_fixed_realizations = TRUE
)
jsonlite::write_json(manifest, file.path(summary_dir, "summary_metadata.json"), pretty = TRUE, auto_unbox = TRUE)
writeLines(capture.output(sessionInfo()), file.path(summary_dir, "sessionInfo.txt"))
cat("SUMMARIZED production_rows=", nrow(metrics), " contrasts=", nrow(contrast_summary), "\n", sep = "")

# Optional v6 audit interval, separate from the paper's mean-only table.
set.seed(5200002L)
v <- metrics$f1[metrics$dgp == 3 & metrics$method == "path_continuous"]
stopifnot(length(v)==100L,all(is.finite(v)))
ci <- quantile(replicate(9999L,mean(sample(v,length(v),replace=TRUE))),c(.025,.975))
write.csv(data.frame(n=100,mean=mean(v),lower=ci[1],upper=ci[2],bootstrap_seed=5200002),
          file.path(summary_dir,"dgp3_f1_mean_audit.csv"),row.names=FALSE)

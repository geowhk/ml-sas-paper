suppressPackageStartupMessages({
  library(yaml)
})

source("02_revision-simulation/R/compare_analysis.R")

parse_cli_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key, call. = FALSE)
    }
    key <- sub("^--", "", key)
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      out[[key]] <- TRUE
      i <- i + 1
    } else {
      out[[key]] <- args[[i + 1]]
      i <- i + 2
    }
  }
  out
}

coerce_config_value <- function(value) {
  if (is.logical(value) || is.numeric(value) || is.integer(value)) {
    return(value)
  }
  if (!is.character(value) || length(value) != 1) {
    return(value)
  }
  if (tolower(value) %in% c("true", "false")) {
    return(tolower(value) == "true")
  }
  numeric_value <- suppressWarnings(as.numeric(value))
  if (!is.na(numeric_value) && grepl("^-?[0-9]+([.][0-9]+)?$", value)) {
    return(numeric_value)
  }
  value
}

args <- parse_cli_args(commandArgs(trailingOnly = TRUE))
config_values <- list()

if (!is.null(args$config)) {
  config_values <- yaml::read_yaml(args$config)
}

if (is.null(config_values$dgp_ids)) {
  dgp_ids <- c("dgp1", "dgp2", "dgp3", "dgp4")
} else {
  dgp_ids <- config_values$dgp_ids
}
config_values$dgp_ids <- NULL

for (nm in setdiff(names(args), "config")) {
  if (nm == "dgp_ids") {
    dgp_ids <- strsplit(args[[nm]], ",", fixed = TRUE)[[1]]
  } else {
    config_values[[nm]] <- coerce_config_value(args[[nm]])
  }
}

dgp_ids <- trimws(dgp_ids)
run_results <- vector("list", length(dgp_ids))
names(run_results) <- dgp_ids

for (dgp_id in dgp_ids) {
  config_values$dgp_id <- dgp_id
  config <- do.call(default_analysis_config, config_values)

  message("============================================================")
  message("DGP batch run")
  message("DGP: ", config$dgp_id)
  message("RF trees: ", config$rf_ntree)
  message("RF mtry: ", config$rf_mtry)
  message("RF nodesize: ", config$rf_nodesize)
  message("Permutations: ", config$perm_nsim)
  message("ML-SAS alternative: ", config$mlsas_alternative)
  message("Workers: ", config$workers)
  message("============================================================")

  run_results[[dgp_id]] <- run_dgp_analysis(config)
  message("Completed ", dgp_id, ". Run directory: ", run_results[[dgp_id]]$dirs$run)
}

message("============================================================")
message("Batch completed.")
message("Run directories:")
for (dgp_id in names(run_results)) {
  message(dgp_id, ": ", run_results[[dgp_id]]$dirs$run)
}
message("============================================================")

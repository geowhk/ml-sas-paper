library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

repo_root <- normalizePath(".")

source(file.path(repo_root, "03_hyperparameter-sensitivity/R/sensitivity_helpers.R"))

theme_revision <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 10),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "vertical",
      panel.grid.minor = element_blank()
    )
}

add_global_observed <- function(results, n_true, n_total) {
  n_false <- n_total - n_true
  results |>
    mutate(global_mlsas_observed = (true_local_mean * n_true + false_local_mean * n_false) / n_total)
}

summarise_global <- function(results, group_vars) {
  results |>
    group_by(across(all_of(group_vars))) |>
    summarise(
      global_mean = mean(global_mlsas_observed, na.rm = TRUE),
      global_sd = sd(global_mlsas_observed, na.rm = TRUE),
      global_q25 = quantile(global_mlsas_observed, 0.25, na.rm = TRUE),
      global_q50 = quantile(global_mlsas_observed, 0.50, na.rm = TRUE),
      global_q75 = quantile(global_mlsas_observed, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

render_main_global <- function(n_true, n_total) {
  base <- file.path(repo_root, "03_hyperparameter-sensitivity/runs/mixed_hyperparameter_path")
  obj <- readRDS(file.path(base, "objects/sensitivity_results.rds"))
  results <- add_global_observed(obj$results, n_true, n_total)

  by_condition <- summarise_global(results, c("proximity_type", "ntree", "mtry", "nodesize"))
  by_ntree <- summarise_global(results, "ntree")
  by_mtry <- summarise_global(results, "mtry")
  by_nodesize <- summarise_global(results, "nodesize")

  write_csv(results, file.path(base, "tables/global_mlsas_observed_raw.csv"))
  write_csv(by_condition, file.path(base, "tables/global_mlsas_observed_by_condition.csv"))
  write_csv(by_ntree, file.path(base, "tables/global_mlsas_observed_by_ntree.csv"))
  write_csv(by_mtry, file.path(base, "tables/global_mlsas_observed_by_mtry.csv"))
  write_csv(by_nodesize, file.path(base, "tables/global_mlsas_observed_by_nodesize.csv"))

  p1 <- by_condition |>
    mutate(mtry = factor(mtry), nodesize = factor(nodesize)) |>
    ggplot(aes(ntree, global_mean, color = mtry, linetype = nodesize, group = interaction(mtry, nodesize))) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = sort(unique(by_condition$ntree))) +
    labs(
      x = "ntree",
      y = "Observed global ML-SAS",
      color = "mtry",
      linetype = "nodesize",
      title = "Observed Global ML-SAS Across Hyperparameter Settings",
      subtitle = "DGP 4.3 mixed condition, path proximity"
    ) +
    theme_revision()

  ggsave(file.path(base, "figures/revision_global_mlsas_by_hyperparameter.png"), p1, width = 10, height = 6, dpi = 300)

  p2 <- bind_rows(
    by_ntree |> mutate(parameter = "ntree", value = as.character(ntree)),
    by_mtry |> mutate(parameter = "mtry", value = as.character(mtry)),
    by_nodesize |> mutate(parameter = "nodesize", value = as.character(nodesize))
  ) |>
    mutate(value = factor(value, levels = unique(value))) |>
    ggplot(aes(value, global_mean, ymin = global_q25, ymax = global_q75)) +
    geom_pointrange(color = "#4c78a8", linewidth = 0.5) +
    facet_wrap(~ parameter, scales = "free_x", nrow = 1) +
    labs(
      x = NULL,
      y = "Observed global ML-SAS",
      title = "Marginal Variation in Observed Global ML-SAS",
      subtitle = "Points are means; vertical ranges show interquartile ranges across replications and other hyperparameters"
    ) +
    theme_revision()

  ggsave(file.path(base, "figures/revision_global_mlsas_marginal_summary.png"), p2, width = 10, height = 4.8, dpi = 300)

  p3 <- results |>
    mutate(
      ntree = factor(ntree, levels = sort(unique(ntree))),
      mtry = factor(mtry, levels = sort(unique(mtry))),
      nodesize = factor(nodesize, levels = sort(unique(nodesize)))
    ) |>
    ggplot(aes(ntree, global_mlsas_observed)) +
    geom_boxplot(width = 0.62, outlier.size = 0.7, linewidth = 0.35, fill = "#d8e7f5", color = "#355c7d") +
    facet_grid(mtry ~ nodesize, labeller = labeller(
      mtry = function(x) paste0("mtry = ", x),
      nodesize = function(x) paste0("min node size = ", x)
    )) +
    labs(
      x = "Number of trees",
      y = "Observed global ML-SAS",
      title = "Observed Global ML-SAS by Hyperparameter Combination",
      subtitle = "Each panel fixes mtry and minimum node size; boxes show the distribution across 100 replications"
    ) +
    theme_revision() +
    theme(
      axis.text.x = element_text(size = 8),
      panel.border = element_rect(color = "grey45", fill = NA, linewidth = 0.5),
      panel.spacing = unit(0.8, "lines")
    )

  ggsave(file.path(base, "figures/revision_global_mlsas_facet_boxplot.png"), p3, width = 11.2, height = 7.8, dpi = 300)

  invisible(list(results = results, by_condition = by_condition, by_ntree = by_ntree, by_mtry = by_mtry, by_nodesize = by_nodesize))
}

render_path_terminal_global <- function(n_true, n_total) {
  base <- file.path(repo_root, "03_hyperparameter-sensitivity/runs/mixed_path_vs_terminal")
  obj <- readRDS(file.path(base, "objects/sensitivity_results.rds"))
  results <- add_global_observed(obj$results, n_true, n_total)

  by_condition <- summarise_global(results, c("proximity_type", "ntree", "mtry", "nodesize"))
  by_proximity <- summarise_global(results, "proximity_type")
  by_proximity_mtry <- summarise_global(results, c("proximity_type", "mtry"))
  by_proximity_nodesize <- summarise_global(results, c("proximity_type", "nodesize"))

  write_csv(results, file.path(base, "tables/global_mlsas_observed_raw.csv"))
  write_csv(by_condition, file.path(base, "tables/global_mlsas_observed_by_condition.csv"))
  write_csv(by_proximity, file.path(base, "tables/global_mlsas_observed_by_proximity.csv"))
  write_csv(by_proximity_mtry, file.path(base, "tables/global_mlsas_observed_by_proximity_mtry.csv"))
  write_csv(by_proximity_nodesize, file.path(base, "tables/global_mlsas_observed_by_proximity_nodesize.csv"))

  p1 <- by_condition |>
    mutate(mtry = factor(mtry), proximity_type = factor(proximity_type, levels = c("path", "terminal"))) |>
    ggplot(aes(nodesize, global_mean, color = proximity_type, linetype = mtry, group = interaction(proximity_type, mtry))) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = sort(unique(by_condition$nodesize))) +
    labs(
      x = "nodesize",
      y = "Observed global ML-SAS",
      color = "Proximity",
      linetype = "mtry",
      title = "Observed Global ML-SAS: Path vs Terminal Proximity",
      subtitle = "DGP 4.3 mixed condition, ntree = 100"
    ) +
    theme_revision()

  ggsave(file.path(base, "figures/revision_global_mlsas_path_vs_terminal.png"), p1, width = 9, height = 5.5, dpi = 300, bg = "white")

  p2 <- results |>
    mutate(
      proximity_type = factor(proximity_type, levels = c("path", "terminal")),
      mtry = factor(mtry, levels = sort(unique(mtry))),
      nodesize = factor(nodesize, levels = sort(unique(nodesize)))
    ) |>
    ggplot(aes(nodesize, global_mlsas_observed, fill = proximity_type)) +
    geom_boxplot(width = 0.68, outlier.size = 0.7, linewidth = 0.35, position = position_dodge(width = 0.75)) +
    facet_wrap(~ mtry, labeller = labeller(mtry = function(x) paste0("mtry = ", x)), nrow = 1) +
    scale_fill_manual(values = c(path = "#1b9e77", terminal = "#d95f02")) +
    labs(
      x = "Minimum node size",
      y = "Observed global ML-SAS",
      fill = "Proximity",
      title = "Observed Global ML-SAS by Proximity Definition",
      subtitle = "DGP 4.3 mixed condition, ntree = 100; boxes show 100 replications"
    ) +
    theme_revision() +
    theme(
      panel.border = element_rect(color = "grey45", fill = NA, linewidth = 0.5),
      panel.spacing = unit(0.9, "lines")
    )

  ggsave(file.path(base, "figures/revision_global_mlsas_path_terminal_boxplot.png"), p2, width = 10.4, height = 5.2, dpi = 300, bg = "white")

  invisible(list(results = results, by_condition = by_condition, by_proximity = by_proximity, by_proximity_mtry = by_proximity_mtry))
}

inputs <- load_mixed_dgp43_inputs()
n_total <- length(inputs$true_cluster)
n_true <- sum(inputs$true_cluster)

main <- render_main_global(n_true, n_total)
path_terminal <- render_path_terminal_global(n_true, n_total)

cat("Global ML-SAS summaries saved.\n")
cat("Main analysis range:", range(main$results$global_mlsas_observed), "\n")
cat("Path-terminal analysis range:", range(path_terminal$results$global_mlsas_observed), "\n")

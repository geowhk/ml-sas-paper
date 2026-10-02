args <- commandArgs(trailingOnly = TRUE)
repo_dir <- if (length(args) >= 1L) normalizePath(args[1], mustWork = TRUE) else normalizePath(".")
if (!length(args) %in% c(2L,3L)) stop("Usage: Rscript scripts/render_section4.R REPO_DIR NEW_OUTPUT_DIR")
output_dir <- args[2]
if (file.exists(output_dir)) stop("Output must not exist")
repo_dir <- if(length(args)==3L) normalizePath(args[3]) else file.path(repo_dir,"reference","figure_inputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(sf)
  library(tidyr)
})

font_family <- "Arial"
colour_blue <- "#0072B2"
colour_orange <- "#D55E00"
colour_green <- "#009E73"
colour_grey <- "#E6E6E6"
colour_dark <- "#333333"

theme_publication <- function(base_size = 9.5) {
  theme_classic(base_size = base_size, base_family = font_family) +
    theme(
      text = element_text(colour = colour_dark),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1, colour = colour_dark),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "plain", margin = margin(3, 3, 3, 3)),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(0.55, "cm"),
      panel.spacing = unit(0.65, "lines"),
      plot.margin = margin(5, 6, 5, 5)
    )
}

theme_map <- function(base_size = 9.5) {
  theme_void(base_size = base_size, base_family = font_family) +
    theme(
      plot.title = element_text(
        size = base_size,
        face = "plain",
        hjust = 0.5,
        margin = margin(0, 0, 4, 0)
      ),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 1),
      legend.key.height = unit(0.35, "cm"),
      legend.key.width = unit(0.55, "cm"),
      plot.margin = margin(2, 2, 2, 2)
    )
}

save_publication <- function(plot, filename, width, height) {
  png_file <- file.path(output_dir, paste0(filename, ".png"))
  pdf_file <- file.path(output_dir, paste0(filename, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, units = "in", dpi = 600, bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, units = "in", device = cairo_pdf, bg = "white")
  c(png = png_file, pdf = pdf_file)
}

make_binary_map <- function(sf_data, flag_column, panel_label) {
  levels <- c("Significant spatial cluster", "Not significant")
  sf_data |>
    mutate(display = factor(.data[[flag_column]], levels = levels)) |>
    ggplot() +
    geom_sf(aes(fill = display), colour = "white", linewidth = 0.075) +
    scale_fill_manual(
      values = c("Significant spatial cluster" = colour_orange, "Not significant" = colour_grey),
      breaks = levels,
      drop = FALSE
    ) +
    coord_sf(expand = FALSE) +
    labs(title = panel_label, fill = NULL) +
    theme_map()
}

manifest_rows <- list()
record_figure <- function(number, filename, role, source, width, height, checks) {
  manifest_rows[[length(manifest_rows) + 1L]] <<- data.frame(
    figure = number,
    png = paste0(filename, ".png"),
    pdf = paste0(filename, ".pdf"),
    role = role,
    approved_source = source,
    width_in = width,
    height_in = height,
    png_dpi = 600,
    title = "none",
    subtitle = "none",
    validation_checks = checks
  )
}

# Figure 1 -----------------------------------------------------------------
fig1_source <- file.path(
  repo_dir,
  "03_hyperparameter-sensitivity/runs/mixed_hyperparameter_path/tables/global_mlsas_observed_raw.csv"
)
fig1_data <- read_csv(fig1_source, show_col_types = FALSE) |>
  mutate(
    ntree = factor(ntree, levels = c(25, 50, 100, 250, 500)),
    mtry = factor(mtry, levels = c(1, 2, 3)),
    nodesize = factor(nodesize, levels = c(1, 5, 10))
  )
stopifnot(nrow(fig1_data) == 4500L, all(table(fig1_data$mtry, fig1_data$nodesize, fig1_data$ntree) == 100L))
p1 <- ggplot(fig1_data, aes(ntree, global_mlsas_observed)) +
  geom_boxplot(
    width = 0.62,
    outlier.size = 0.35,
    outlier.alpha = 0.5,
    linewidth = 0.3,
    fill = "#D9EEF7",
    colour = colour_blue
  ) +
  facet_grid(
    mtry ~ nodesize,
    labeller = labeller(
      mtry = function(x) paste0("mtry = ", x),
      nodesize = function(x) paste0("min. node size = ", x)
    )
  ) +
  labs(x = "Number of trees", y = "Observed global ML-SAS") +
  theme_publication(8.5) +
  theme(axis.text.x = element_text(size = 6.8))
save_publication(p1, "Figure_01_hyperparameter_robustness", 7.4, 6.2)
record_figure(
  1, "Figure_01_hyperparameter_robustness", "Path-proximity hyperparameter robustness",
  fig1_source, 7.4, 6.2, "4500 rows; 100 RF realizations per 45 conditions"
)

# Figure 2 -----------------------------------------------------------------
fig2_source <- file.path(
  repo_dir,
  "03_hyperparameter-sensitivity/runs/mixed_path_vs_terminal/tables/global_mlsas_observed_raw.csv"
)
fig2_data <- read_csv(fig2_source, show_col_types = FALSE) |>
  mutate(
    proximity_type = factor(proximity_type, levels = c("path", "terminal"), labels = c("Path", "Terminal-node")),
    nodesize = factor(nodesize, levels = c(1, 5, 10)),
    mtry = factor(mtry, levels = c(1, 2, 3))
  )
stopifnot(nrow(fig2_data) == 1800L, all(table(fig2_data$proximity_type, fig2_data$mtry, fig2_data$nodesize) == 100L))
p2 <- ggplot(fig2_data, aes(nodesize, global_mlsas_observed, fill = proximity_type)) +
  geom_boxplot(
    width = 0.66,
    position = position_dodge(width = 0.76),
    outlier.size = 0.45,
    outlier.alpha = 0.5,
    linewidth = 0.3
  ) +
  facet_wrap(~mtry, nrow = 1, labeller = labeller(mtry = function(x) paste0("mtry = ", x))) +
  scale_fill_manual(values = c("Path" = colour_blue, "Terminal-node" = colour_orange)) +
  labs(x = "Minimum node size", y = "Observed global ML-SAS", fill = NULL) +
  theme_publication()
save_publication(p2, "Figure_02_path_vs_terminal", 7.4, 3.7)
record_figure(
  2, "Figure_02_path_vs_terminal", "Path versus terminal-node proximity",
  fig2_source, 7.4, 3.7, "1800 rows; 100 independent RF realizations per 18 conditions"
)

# Figures 3 and 4 ----------------------------------------------------------
global_source <- file.path(repo_dir, "04_current-statistics-comparison/runs/revision_v2")
stats <- read_csv(file.path(global_source, "tables/permutation_statistics.csv"), show_col_types = FALSE)
correlations <- read_csv(file.path(global_source, "tables/correlation_summary.csv"), show_col_types = FALSE)
stopifnot(nrow(stats) == 1000L)

continuous_long <- stats |>
  select(continuous_mlsas, moran_i, geary_c) |>
  pivot_longer(c(moran_i, geary_c), names_to = "comparison_key", values_to = "statistic") |>
  mutate(
    panel = factor(
      comparison_key,
      levels = c("moran_i", "geary_c"),
      labels = c("(a) Moran's I", "(b) Geary's c")
    )
  )
continuous_labels <- correlations |>
  filter(comparison %in% c("Moran's I", "Geary's c")) |>
  mutate(
    panel = factor(comparison, levels = c("Moran's I", "Geary's c"), labels = c("(a) Moran's I", "(b) Geary's c")),
    label = sprintf("Pearson r = %.3f", pearson_r)
  )
stopifnot(
  abs(cor(stats$continuous_mlsas, stats$moran_i) - correlations$pearson_r[correlations$comparison == "Moran's I"]) < 1e-12,
  abs(cor(stats$continuous_mlsas, stats$geary_c) - correlations$pearson_r[correlations$comparison == "Geary's c"]) < 1e-12
)
p3 <- ggplot(continuous_long, aes(continuous_mlsas, statistic)) +
  geom_point(size = 0.65, alpha = 0.42, colour = colour_dark) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = colour_blue, fill = "#BFDDEE", linewidth = 0.65) +
  geom_text(
    data = continuous_labels,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.08,
    vjust = 1.35,
    inherit.aes = FALSE,
    family = font_family,
    size = 3.0
  ) +
  facet_wrap(~panel, scales = "free_y", nrow = 1) +
  labs(x = "Global ML-SAS", y = NULL) +
  theme_publication() +
  theme(panel.border = element_rect(colour = "#666666", fill = NA, linewidth = 0.35))
save_publication(p3, "Figure_03_continuous_benchmark", 7.4, 3.5)
record_figure(
  3, "Figure_03_continuous_benchmark", "Continuous global benchmark",
  global_source, 7.4, 3.5, "n=1000; r(Moran)=0.472681; r(Geary)=-0.468342"
)

bb_row <- correlations |> filter(comparison == "BB join-count")
stopifnot(nrow(bb_row) == 1L, abs(cor(stats$categorical_mlsas, stats$bb_join_count) - bb_row$pearson_r) < 1e-12)
p4 <- ggplot(stats, aes(categorical_mlsas, bb_join_count)) +
  geom_point(size = 0.7, alpha = 0.42, colour = colour_dark) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = colour_blue, fill = "#BFDDEE", linewidth = 0.65) +
  annotate(
    "text", x = Inf, y = Inf, label = sprintf("Pearson r = %.3f", bb_row$pearson_r),
    hjust = 1.08, vjust = 1.35, family = font_family, size = 3.0
  ) +
  labs(x = "Global ML-SAS", y = "BB join-count") +
  theme_publication() +
  theme(panel.border = element_rect(colour = "#666666", fill = NA, linewidth = 0.35))
save_publication(p4, "Figure_04_bb_join_count", 5.1, 4.0)
record_figure(
  4, "Figure_04_bb_join_count", "Binary global benchmark against BB join-count",
  global_source, 5.1, 4.0, "n=1000; BB uses symmetric binary adjacency; r=0.874269"
)

# Figure 5 -----------------------------------------------------------------
dgp1_file <- file.path(repo_dir, "02_revision-simulation/sim_dgp/sim_dgp1.rds")
dgp1 <- readRDS(dgp1_file)
fig5_prod <- file.path(
  repo_dir,
  "02_revision-simulation/runs/revision_v2/production/dgp1_section4_v2_ntree100_perm99999_20260718_133406"
)
fig5_correction <- file.path(
  repo_dir,
  "02_revision-simulation/runs/revision_v2/corrections/dgp1_local_geary_greater_perm99999_20260719"
)
lisa1 <- read_csv(file.path(fig5_correction, "tables/dgp1_univariate_lisa_results_corrected.csv"), show_col_types = FALSE)
mlsas1 <- read_csv(file.path(fig5_prod, "tables/dgp1_univariate_mlsas_local_results.csv"), show_col_types = FALSE)
flags1 <- tibble(
  id = lisa1$id,
  moran_flag = lisa1$moran_flag,
  geary_flag = lisa1$geary_flag,
  gi_flag = lisa1$gi_flag,
  mlsas_flag = mlsas1$cluster_flag_strict
)
map1 <- left_join(dgp1, flags1, by = "id")
fig5_counts <- c(
  sum(map1$moran_flag == "Significant spatial cluster"),
  sum(map1$geary_flag == "Significant spatial cluster"),
  sum(map1$gi_flag == "Significant spatial cluster"),
  sum(map1$mlsas_flag == "Significant spatial cluster")
)
stopifnot(identical(as.integer(fig5_counts), c(61L, 83L, 66L, 69L)))
p5 <- (
  make_binary_map(map1, "moran_flag", "(a) Local Moran's I") +
    make_binary_map(map1, "geary_flag", "(b) Local Geary's c") +
    make_binary_map(map1, "gi_flag", "(c) Getis-Ord Gi*") +
    make_binary_map(map1, "mlsas_flag", "(d) ML-SAS")
) + plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "bottom")
save_publication(p5, "Figure_05_dgp1_local_comparison", 7.4, 6.0)
record_figure(
  5, "Figure_05_dgp1_local_comparison", "DGP 1 local comparison",
  paste(fig5_prod, fig5_correction, sep = " + "), 7.4, 6.0,
  "significant cells: Moran=61; Geary=83; Gi*=66; ML-SAS=69"
)

# Figure 6 -----------------------------------------------------------------
dgp31_file <- file.path(repo_dir, "02_revision-simulation/sim_dgp/sim_dgp3.1.rds")
dgp31 <- readRDS(dgp31_file)
fig6_source <- file.path(
  repo_dir,
  "02_revision-simulation/runs/revision_v2/production/dgp3.1_section4_v2_ntree100_perm99999_20260718_152439"
)
mg31 <- read_csv(file.path(fig6_source, "tables/dgp3.1_continuous_multigeary_local_results.csv"), show_col_types = FALSE)
ml31 <- read_csv(file.path(fig6_source, "tables/dgp3.1_continuous_mlsas_local_results.csv"), show_col_types = FALSE)
map31 <- left_join(
  dgp31,
  tibble(id = mg31$id, multigeary_flag = mg31$cluster_flag, mlsas_flag = ml31$cluster_flag_strict),
  by = "id"
)
stopifnot(
  sum(map31$multigeary_flag == "Significant spatial cluster") == 33L,
  sum(map31$mlsas_flag == "Significant spatial cluster") == 31L
)
p6 <- (
  make_binary_map(map31, "multigeary_flag", "(a) Multivariate Local Geary") +
    make_binary_map(map31, "mlsas_flag", "(b) ML-SAS")
) + plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "bottom")
save_publication(p6, "Figure_06_dgp3_local_comparison", 7.4, 3.45)
record_figure(
  6, "Figure_06_dgp3_local_comparison", "DGP 3 continuous local comparison",
  fig6_source, 7.4, 3.45, "significant cells: Multivariate Local Geary=33; ML-SAS=31"
)

# Figure 7 -----------------------------------------------------------------
dgp43_file <- file.path(repo_dir, "02_revision-simulation/sim_dgp/sim_dgp4.3.rds")
dgp43 <- readRDS(dgp43_file)
fig7_source <- file.path(
  repo_dir,
  "02_revision-simulation/runs/revision_v2/production/dgp4.3_section4_v2_ntree100_perm99999_20260718_155829"
)
mg43 <- read_csv(file.path(fig7_source, "tables/dgp4.3_continuous_multigeary_local_results.csv"), show_col_types = FALSE)
ml43c <- read_csv(file.path(fig7_source, "tables/dgp4.3_continuous_mlsas_local_results.csv"), show_col_types = FALSE)
ml43m <- read_csv(file.path(fig7_source, "tables/dgp4.3_mixed_mlsas_local_results.csv"), show_col_types = FALSE)
map43 <- left_join(
  dgp43,
  tibble(
    id = mg43$id,
    multigeary_flag = mg43$cluster_flag,
    mlsas_continuous_flag = ml43c$cluster_flag_strict,
    mlsas_mixed_flag = ml43m$cluster_flag_strict
  ),
  by = "id"
)
stopifnot(
  sum(map43$multigeary_flag == "Significant spatial cluster") == 42L,
  sum(map43$mlsas_continuous_flag == "Significant spatial cluster") == 19L,
  sum(map43$mlsas_mixed_flag == "Significant spatial cluster") == 46L
)
p7 <- (
  make_binary_map(map43, "multigeary_flag", "(a) Multivariate Local Geary") +
    make_binary_map(map43, "mlsas_continuous_flag", "(b) ML-SAS: continuous") +
    make_binary_map(map43, "mlsas_mixed_flag", "(c) ML-SAS: mixed")
) + plot_layout(ncol = 3, guides = "collect") & theme(legend.position = "bottom")
save_publication(p7, "Figure_07_dgp4_mixed_comparison", 7.4, 3.0)
record_figure(
  7, "Figure_07_dgp4_mixed_comparison", "DGP 4 nonlinear mixed-type comparison",
  fig7_source, 7.4, 3.0, "significant cells: Multivariate Local Geary=42; ML-SAS continuous=19; ML-SAS mixed=46"
)

# Manifest and deterministic checks ----------------------------------------
manifest <- bind_rows(manifest_rows)
manifest$png_bytes <- file.info(file.path(output_dir, manifest$png))$size
manifest$pdf_bytes <- file.info(file.path(output_dir, manifest$pdf))$size
manifest$png_md5 <- unname(tools::md5sum(file.path(output_dir, manifest$png)))
manifest$pdf_md5 <- unname(tools::md5sum(file.path(output_dir, manifest$pdf)))
stopifnot(
  nrow(manifest) == 7L,
  all(manifest$png_bytes > 0),
  all(manifest$pdf_bytes > 0),
  !anyDuplicated(manifest$png),
  !anyDuplicated(manifest$pdf)
)
write_csv(manifest, file.path(output_dir, "figure_manifest.csv"))

writeLines(
  c(
    "Section 4 publication figure render completed.",
    paste0("created_at_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    "figures=7",
    "png_dpi=600",
    "pdf_device=cairo_pdf",
    "overall_titles=none",
    "subtitles=none",
    "font=Arial",
    "palette=colourblind-safe blue/orange plus neutral grey"
  ),
  file.path(output_dir, "RENDER_METADATA.txt")
)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
message("Stage 4 figures written to: ", output_dir)

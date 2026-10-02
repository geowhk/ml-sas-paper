args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) stop("Usage: Rscript scripts/render_empirical_results.R CONDITIONAL_RUN_DIR ANALYSIS_INPUTS_RDS NEW_OUTPUT_DIR")
run_dir <- normalizePath(args[1],mustWork=TRUE)
input_rds <- normalizePath(args[2],mustWork=TRUE)
output_dir <- args[3]
if (file.exists(output_dir)) stop("Output must not exist")
dir.create(output_dir,recursive=TRUE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
  library(patchwork)
  library(ggrepel)
})

comparison_csv <- file.path(run_dir, "tables", "continuous_vs_mixed_local.csv")
global_csv <- file.path(run_dir, "tables", "global_mlsas.csv")
inputs <- readRDS(input_rds)
comparison <- utils::read.csv(comparison_csv, stringsAsFactors = FALSE)
global_results <- utils::read.csv(global_csv, stringsAsFactors = FALSE)

stopifnot(
  inherits(inputs$sf, "sf"),
  nrow(inputs$sf) == 337L,
  nrow(comparison) == 337L,
  !anyDuplicated(inputs$sf$code),
  !anyDuplicated(comparison$code),
  setequal(inputs$sf$code, comparison$code),
  all(c("continuous_class", "mixed_class", "class_changed") %in% names(comparison)),
  identical(unique(global_results$nperm), 99999L)
)

map_data <- inputs$sf
matched_rows <- match(map_data$code, comparison$code)
stopifnot(!anyNA(matched_rows))
map_data$continuous_class <- comparison$continuous_class[matched_rows]
map_data$mixed_class <- comparison$mixed_class[matched_rows]

transition_levels <- c(
  "Similarity in both",
  "New similarity in mixed",
  "Similarity lost in mixed",
  "Dissimilarity/outlier in both",
  "Nonsignificant in both"
)

map_data$transition <- with(
  map_data,
  ifelse(
    continuous_class == "Similarity cluster" & mixed_class == "Similarity cluster",
    transition_levels[1],
    ifelse(
      continuous_class == "Not significant" & mixed_class == "Similarity cluster",
      transition_levels[2],
      ifelse(
        continuous_class == "Similarity cluster" & mixed_class == "Not significant",
        transition_levels[3],
        ifelse(
          continuous_class == "Dissimilarity/outlier" & mixed_class == "Dissimilarity/outlier",
          transition_levels[4],
          ifelse(
            continuous_class == "Not significant" & mixed_class == "Not significant",
            transition_levels[5],
            NA_character_
          )
        )
      )
    )
  )
)

stopifnot(!anyNA(map_data$transition))
map_data$transition <- factor(map_data$transition, levels = transition_levels)
transition_counts <- table(map_data$transition)
stopifnot(sum(transition_counts) == 337L)

transition_colours <- c(
  "Similarity in both" = "#D55E00",
  "New similarity in mixed" = "#E69F00",
  "Similarity lost in mixed" = "#CC79A7",
  "Dissimilarity/outlier in both" = "#0072B2",
  "Nonsignificant in both" = "#E6E6E6"
)

transition_labels <- sprintf(
  "%s (n = %d)",
  transition_levels,
  as.integer(transition_counts)
)
names(transition_labels) <- transition_levels

label_lookup <- c(
  "GM0163" = "Hellendoorn",
  "GM0297" = "Zaltbommel",
  "GM0335" = "Montfoort",
  "GM0394" = "Haarlemmermeer",
  "GM0590" = "Papendrecht",
  "GM1842" = "Midden-Delfland",
  "GM1901" = "Bodegraven-Reeuwijk",
  "GM1916" = "Leidschendam-Voorburg"
)
stopifnot(all(names(label_lookup) %in% map_data$code))

make_label_data <- function(x) {
  points <- suppressWarnings(st_point_on_surface(x))
  coords <- st_coordinates(points)
  data.frame(
    code = x$code,
    label = unname(label_lookup[x$code]),
    X = coords[, "X"],
    Y = coords[, "Y"],
    stringsAsFactors = FALSE
  )
}

main_labels <- make_label_data(map_data[map_data$code %in% c("GM0163", "GM0297"), ])
west_codes <- c("GM0335", "GM0394", "GM0590", "GM1842", "GM1901", "GM1916")
west_focus <- map_data[map_data$code %in% west_codes, ]
west_labels <- make_label_data(west_focus)
west_bbox <- st_bbox(st_buffer(st_union(west_focus), dist = 17000))
west_extent <- st_as_sfc(west_bbox)
national_bbox <- st_bbox(map_data)
national_width <- unname(national_bbox["xmax"] - national_bbox["xmin"])
main_xlim <- c(
  unname(national_bbox["xmin"] - 0.72 * national_width),
  unname(national_bbox["xmax"])
)

font_family <- "Arial"

base_map <- ggplot(map_data) +
  geom_sf(aes(fill = transition), colour = "white", linewidth = 0.075) +
  geom_sf(
    data = west_extent,
    inherit.aes = FALSE,
    fill = NA,
    colour = "#333333",
    linewidth = 0.48
  ) +
  ggrepel::geom_label_repel(
    data = main_labels,
    aes(x = X, y = Y, label = label),
    inherit.aes = FALSE,
    size = 2.4,
    family = font_family,
    colour = "#222222",
    fill = "white",
    label.size = 0.15,
    box.padding = 0.25,
    point.padding = 0.2,
    min.segment.length = 0,
    segment.size = 0.25,
    seed = 20260720,
    max.overlaps = Inf
  ) +
  coord_sf(
    xlim = main_xlim,
    ylim = c(national_bbox["ymin"], national_bbox["ymax"]),
    expand = FALSE,
    datum = NA
  ) +
  scale_fill_manual(
    values = transition_colours,
    breaks = transition_levels,
    labels = transition_labels,
    drop = FALSE
  ) +
  labs(fill = NULL) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  theme_void(base_size = 9.5, base_family = font_family) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 7.6, colour = "#333333"),
    legend.key.height = grid::unit(0.32, "cm"),
    legend.key.width = grid::unit(0.50, "cm"),
    legend.spacing.x = grid::unit(0.12, "cm"),
    legend.spacing.y = grid::unit(0.03, "cm"),
    plot.margin = margin(2, 2, 2, 2),
    plot.background = element_rect(fill = "white", colour = NA)
  )

inset_map <- ggplot(map_data) +
  geom_sf(aes(fill = transition), colour = "white", linewidth = 0.16) +
  ggrepel::geom_label_repel(
    data = west_labels,
    aes(x = X, y = Y, label = label),
    inherit.aes = FALSE,
    size = 2.35,
    family = font_family,
    colour = "#222222",
    fill = "white",
    label.size = 0.12,
    box.padding = 0.28,
    point.padding = 0.18,
    min.segment.length = 0,
    segment.size = 0.24,
    seed = 20260720,
    force = 10,
    max.overlaps = Inf
  ) +
  coord_sf(
    xlim = c(west_bbox["xmin"], west_bbox["xmax"]),
    ylim = c(west_bbox["ymin"], west_bbox["ymax"]),
    expand = FALSE,
    datum = NA
  ) +
  scale_fill_manual(values = transition_colours, drop = FALSE) +
  annotate(
    "label",
    x = west_bbox["xmin"] + 0.02 * (west_bbox["xmax"] - west_bbox["xmin"]),
    y = west_bbox["ymax"] - 0.03 * (west_bbox["ymax"] - west_bbox["ymin"]),
    label = "Western Netherlands",
    hjust = 0,
    vjust = 1,
    family = font_family,
    size = 2.4,
    label.size = 0,
    fill = "white",
    colour = "#333333"
  ) +
  theme_void(base_size = 8.5, base_family = font_family) +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "#555555", fill = NA, linewidth = 0.35),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

p <- base_map + patchwork::inset_element(
  inset_map,
  left = 0.012,
  bottom = 0.590,
  right = 0.405,
  top = 0.988,
  align_to = "full"
)

filename <- "Figure_08_NLD_muni_local_clusters"
png_file <- file.path(output_dir, paste0(filename, ".png"))
pdf_file <- file.path(output_dir, paste0(filename, ".pdf"))
width_in <- 7.2
height_in <- 6.3

ggsave(
  png_file,
  p,
  width = width_in,
  height = height_in,
  units = "in",
  dpi = 600,
  bg = "white"
)
ggsave(
  pdf_file,
  p,
  width = width_in,
  height = height_in,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

manifest <- data.frame(
  figure = 8L,
  png = basename(png_file),
  pdf = basename(pdf_file),
  role = "Mixed-specification local ML-SAS classifications and transitions from the continuous-only specification",
  approved_source_geometry = normalizePath(input_rds),
  approved_source_comparison = normalizePath(comparison_csv),
  variables = "population; dwelling_value; employment_rate; ordered urbanity",
  sample_n = nrow(map_data),
  queen_contiguity = TRUE,
  permutations = unique(global_results$nperm),
  inference = "two-sided",
  multiplicity_adjustment = "Benjamini-Hochberg",
  alpha = 0.05,
  similarity_both = unname(transition_counts[transition_levels[1]]),
  new_similarity_mixed = unname(transition_counts[transition_levels[2]]),
  similarity_lost_mixed = unname(transition_counts[transition_levels[3]]),
  dissimilarity_both = unname(transition_counts[transition_levels[4]]),
  nonsignificant_both = unname(transition_counts[transition_levels[5]]),
  labelled_municipalities = paste(unname(label_lookup), collapse = "; "),
  inset_location = "upper-left outside national map",
  inset_extent_box = TRUE,
  width_in = width_in,
  height_in = height_in,
  png_dpi = 600L,
  title = "none",
  subtitle = "none",
  png_bytes = file.info(png_file)$size,
  pdf_bytes = file.info(pdf_file)$size,
  png_md5 = unname(tools::md5sum(png_file)),
  pdf_md5 = unname(tools::md5sum(pdf_file))
)
utils::write.csv(manifest, file.path(output_dir, "figure_manifest.csv"), row.names = FALSE)

writeLines(
  c(
    "Section 5 transition-aware publication figure render completed.",
    paste0("created_at_utc=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    "figures=1",
    "png_dpi=600",
    "pdf_device=cairo_pdf",
    "overall_titles=none",
    "subtitles=none",
    "font=Arial",
    "palette=colourblind-safe Okabe-Ito plus neutral grey",
    paste0("local_inference=two-sided permutation p-values with BH adjustment; nperm=", manifest$permutations),
    sprintf(
      "transition_counts=%d similarity in both; %d new similarity; %d lost similarity; %d dissimilarity/outlier in both; %d nonsignificant in both",
      manifest$similarity_both,
      manifest$new_similarity_mixed,
      manifest$similarity_lost_mixed,
      manifest$dissimilarity_both,
      manifest$nonsignificant_both
    ),
    paste0("labelled_municipalities=", manifest$labelled_municipalities)
  ),
  file.path(output_dir, "RENDER_METADATA.txt")
)
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))

message("Transition-aware Section 5 figure written to: ", output_dir)

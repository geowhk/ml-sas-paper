env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) as.integer(value) else as.integer(default)
}

list(
  run_label = "dgp4_lovo_revision_v2",
  dgp_id = "dgp4.3",
  variables = c("x1", "x2", "cat1"),
  primary_rf_seed = 123L,
  primary_synthetic_seed = 100123L,
  permutation_seed = 123L,
  rf_ntree = env_int("DGP43_LOVO_NTREE", 100L),
  rf_mtry = env_int("DGP43_LOVO_MTRY", 1L),
  rf_nodesize = env_int("DGP43_LOVO_NODESIZE", 1L),
  nperm = env_int("DGP43_LOVO_NPERM", 99999L),
  workers = env_int("DGP43_LOVO_WORKERS", 24L),
  stability_reps = env_int("DGP43_LOVO_REPS", 100L),
  stability_seed_start = 20260715L,
  stability_ntree = env_int("DGP43_LOVO_STABILITY_NTREE", 100L),
  run_primary_inference = FALSE,
  alpha = 0.05,
  p_adjust_method = "BH",
  style = "W",
  queen = TRUE
)


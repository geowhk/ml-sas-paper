env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) as.integer(value) else as.integer(default)
}

list(
  run_label = "NLD_muni_reduced_urbanity",
  continuous_vars = c("population", "dwelling_value", "employment_rate"),
  categorical_var = "urbanity",
  primary_seed = 20260714L,
  primary_synthetic_seed = 20260814L,
  permutation_seed = 20260729L,
  main_ntree = env_int("REDUCED_MAIN_NTREE", 1000L),
  nperm = env_int("REDUCED_NPERM", 4999L),
  workers = env_int("REDUCED_WORKERS", 1L),
  ablation_reps = env_int("REDUCED_ABLATION_REPS", 5L),
  ablation_ntree = env_int("REDUCED_ABLATION_NTREE", 500L),
  ablation_seeds = 20260714L + 0:4,
  alpha = 0.05,
  p_adjust_method = "BH"
)

# Reproduction guide

Run all commands from the repository root.

## 1. Install

Use R 4.5.2 and Python 3.13 for the checked environment. R package versions are
locked in `renv.lock`; Python versions are in `environment/requirements-figures.txt`.
From the repository root:

```sh
Rscript scripts/setup.R
python3 -m venv .venv
.venv/bin/python -m pip install -r environment/requirements-figures.txt
```

The R setup creates a project library and activates it for subsequent Rscript calls.
It does not silently upgrade unavailable locked packages. Source installations need
compilers and spatial system libraries; read [installation troubleshooting](INSTALL.md).
First-time package download/compilation is separate from analysis runtime.

## 2. Reproduce saved results — the quick check

```sh
.venv/bin/python scripts/reproduce.py --mode saved --out runs/saved
```

This regenerates final F1 summaries and bootstrap intervals, Appendix C diagnostics,
and all ten figures. Empirical inputs come from the pinned `tmap::NLD_muni` dataset;
their attributes, geometry and queen adjacency must match v6 fingerprints for all
337 analytical units. Synthetic local evidence and derived empirical result tables
are included; no manuscript, review history, original municipal RDS/GeoJSON or fonts
are bundled. Saved mode does **not** refit models.

## 3. Recompute models, inference, tables and figures — slow

```sh
.venv/bin/python scripts/reproduce.py --mode full --out runs/full --workers 1 --parallelism 1
```

Full mode runs representative DGP analyses, sensitivity analyses, the path/terminal
comparison, the existing-statistics comparison, the empirical base/global/conditional
pipeline, and all 300 repeated datasets (800 method-dataset results). It then builds
the final summaries, Appendix C and figures from these new results. It compares the
final F1 CSVs, empirical results and Appendix C CSVs to v6 at absolute tolerance 1e-12.
A mismatch stops the run; saved result values are never inserted to make a new run pass.

The original research computations used AWS. Full recomputation is intended for a
server or cloud environment with sufficient CPU, memory and uninterrupted runtime;
saved mode provides a practical local check.

**Plan for hours to days, especially with one worker.** The archived repeated-DGP batch alone took 19 h 32 min;
that is an old-run measurement, not a prediction for your machine. Read the runtime
notes before starting. `--parallelism` controls independent repeated datasets;
`--workers` controls workers inside the other analyses. Both default to one to avoid
oversubscribing a laptop. Seeds, tree counts and permutation counts stay fixed.
The existing-statistics comparison intentionally uses its original 1,000 permutations;
final repeated/empirical inference uses 99,999.

Each stage writes a log. `reproduction.json` records commands, exit status and elapsed
seconds. It says PASS only when the requested mode completes, including comparisons.
To resume after an interruption:

```sh
.venv/bin/python scripts/reproduce.py --mode full --out runs/full --workers 1 --parallelism 1 --resume
```

Resumption rejects changed code, reference inputs or mode. Completed stages are retained.
Sensitivity and path/terminal analyses retain completed condition checkpoints.
Other incomplete non-repeated stage outputs move to `incomplete/` and that stage restarts.
The repeated runner retains completed replications, archives partial replications,
reruns those with the same seeds, and validates stages again. Do not edit a running output.
To rerun after a code change, choose a new output directory.

## Outputs and paper mapping

| Paper item | Output in a full run |
|---|---|
| Table 1 design | Frozen specification in `analysis/repeated/repeated_simulation_spec.yml`; generator preflight |
| Table 2 | `repeated/summary/dgp2_null_operating_characteristics.csv` |
| Tables 3–4 | `repeated/summary/dgp3_dgp4_method_summaries.csv` |
| Paired differences/intervals | `repeated/summary/prespecified_paired_contrasts.csv` |
| Tables C1–C3 | `appendix/dgp3_pvalue_tail_rates.csv`, `category_diagnostics_summary.csv`, `category_pair_auc_summary.csv` |
| Figures 1–4 | `section4/`: historical filenames 03, 04, 01, 02, respectively |
| Figures 5–7, B1 | `dgp-figures/` |
| Figure 8 | `figure8/` |
| Figure 9 | `figure9/` |
| Section 5 numerical results | `empirical/conditional/tables/`, `empirical/global/tables/` |

Extra historical DGP plots emitted by the R Section 4 renderer are not the final
Figures 5–7; use the Python outputs. The published 600 dpi assets were used for visual
comparison. R/Python generate PDF and PNG; font and renderer versions can change layout.

Individual analyses remain available via `python3 scripts/run_analysis.py --help`.
The legacy LOVO experiment is available separately; it is not a required dependency of
the final figure/table mapping. Quick regression checks:

```sh
python3 tests/verify_pipeline.py
Rscript tests/check_empirical_workers.R .
python3 tests/check_resume.py
```

## Final calculation rules

- Non-null F1 is `2 TP / (2 TP + FP + FN)`. No discoveries contribute zero.
- Precision is undefined without discoveries. Archived DGP 3 path-continuous precision
  uses 99 datasets; F1, recall, FPR and rank AUC use all 100.
- Null DGP 2 F1 remains undefined. This is separate from the non-null no-discovery case.
- Local Geary rank AUC uses the negated raw statistic once. Fresh summaries reject the
  historical `1-AUC` correction option.
- Bootstrap F1 uses all pairs while retaining v6's historical RNG progression for
  unrelated contrasts. Historical F1 replay is isolated from fresh analyses.
- Empirical input paths are explicit; archived absolute paths are never followed.

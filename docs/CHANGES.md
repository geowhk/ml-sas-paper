# Changes from frozen sources

1. Repeated protocol moved to `analysis/repeated`; root resolution rebased. Original DGP
   source bytes and its asserted SHA-256 are preserved.
2. `make_method_metrics` computes non-null F1 directly from TP/FP/FN. Undefined null-target
   F1 and conditional precision retain their historical meaning.
3. Fresh summaries reject legacy AUC transformation. The corrected score orientation is
   used at computation time. No old postprocessing runs after a new experiment.
4. Fresh summaries include `f1_defined_reps` and the optional F1 mean audit interval.
   For each F1 contrast, bootstrap estimates use the corrected values from the entry RNG
   state. The subsequent RNG state advances as in the historical finite-F1 sample stream.
   This reproduces v6's policy of retaining unrelated archived intervals, without substituting
   archived intervals into a new run. Seeds remain 5200001 and 5200002; 9,999 resamples.
5. Historical replay now takes explicit repository/output arguments and refuses overwrites.
6. Conditional empirical analysis receives BASE_RUN_DIR explicitly. Figure 9 takes its
   input RDS explicitly; neither follows the old config's absolute source path.
7. A thin analysis launcher supplies output paths without changing analytical parameters.
   Sensitivity/comparison accept an optional output override. New output guards were added.
8. Figure scripts use selected inputs and explicit outputs. Figure 8 rendering was extracted
   from the original script without its duplicate historical DGP rendering stages.
9. Historical code comments/config source paths identifying a local Downloads folder were
   replaced by the preserved relative source path. No original preservation files were edited.

A package API and a new execution of every full-sized experiment remain outside the
completed validation scope. An isolated R-library restore and Appendix C diagnostic replay
have now passed. This preparation does not certify the paper's mathematical proofs.

## Complete execution paths

- Added `reproduce.py` with saved and full modes, per-stage logs/timings, strict comparisons,
  code/input-bound resumption, and retained incomplete outputs.
- Added Appendix C export/analysis using the same saved-local diagnostic formulas, with
  the obsolete diagnostic F1 expression replaced by the v6 count-based definition.
- Added package-sourced empirical preparation with frozen attribute, geometry and adjacency
  fingerprints. Only fingerprints and derived result tables are bundled, not source geometry.
- Added fresh-result figure assembly. It derives sensitivity global values from the original
  weighted local means and maps newly produced representative outputs to renderer filenames.
- Fixed the original empirical inference `cut(..., breaks=1)` failure for a single worker or
  a single permutation. A regression check compares serial and parallel results.
- Fixed empirical tree chunk grouping at its original 48-worker value while allowing fewer
  executing workers. This preserves the floating-point accumulation order; mixed proximity
  reproduced the saved matrix with maximum absolute difference zero in the checked rerun.
- Worker overrides change scheduling only. No permutation counts, seeds, sample sizes,
  tree counts, targets or inferential rules were reduced for speed.
- Sensitivity condition evaluation now consistently uses `future_lapply(future.seed=TRUE)`,
  including a sequential plan. The original parallel configuration used L'Ecuyer-CMRG;
  falling back to plain `map` at one worker would change the RNG kind despite identical
  per-condition seeds. The worker-invariance test anchors its first condition to v6.

- Sensitivity/path-terminal runs checkpoint each completed condition atomically and reject
  resumption if code or configuration changes. Repeated runs resume at replication boundaries.

# Release status — 2026-09-30

Repository: `geowhk/ml-sas-paper`. DGP publication permission was confirmed by
the maintainer on 2026-10-02. R package code is maintained separately in
[ml-sas-r](https://github.com/geowhk/ml-sas-r).

The validation results below describe release preparation on 2026-09-30;
publication checks on the separated repository are recorded in VALIDATION.json.

| Check | Result and scope |
|---|---|
| Frozen preservation integrity | PASS, 3,965 read-only checks; frozen source not modified |
| Saved F1 and aggregation | PASS, five final CSVs match v6 at absolute tolerance 1e-12 |
| Metric edge cases | PASS, 600 non-null confusion matrices; conditional precision; AUC orientation |
| Generator/preflight | PASS, source checksum, fixed DGP 1–4 values/geometry and seed determinism |
| R dependency restoration | PASS, isolated library with 169 pinned packages; host native compilers/libraries required; not a clean OS installation |
| Saved end-to-end mode | PASS, F1 bootstrap, C1–C3, empirical input reconstruction and all ten figures |
| Appendix C | PASS, all three CSVs regenerated from 240,000 local evidence rows and compared with v6 |
| Empirical full model/inference | PASS, fresh models plus 99,999 global and conditional permutations; all eight final CSVs match v6 numerically with maximum observed difference zero; 14 min 25 s with one worker |
| Fresh repeated simulation | PASS for DGP 3 replication 91, both methods, 99,999 permutations; all metric columns match v6 within 1e-12; 27 min 32 s |
| Worker handling | PASS, empirical serial/parallel regression and sensitivity first-condition serial/parallel equality anchored to v6 |
| Resume handling | PASS, overwrite/config guards, completed-output integrity checks, and sensitivity condition checkpoint reuse |
| Fresh-input figure assembly | PASS against a fixture assembled from archived analysis outputs; not evidence that every model has been rerun |
| Figure appearance | Ten figures inspected alongside frozen 600 dpi TIFFs; no pixel-identity claim |
| All 300 repeated datasets recomputed | NOT RUN in this preparation |
| Entire sensitivity grids and complete full mode | NOT RUN in this preparation; full execution route is implemented but not yet validated end to end at paper scale |
| Other operating systems | NOT TESTED |

See `VALIDATION.json`, `EMPIRICAL_RERUN.json`, `REPEATED_SAMPLE_RERUN.json`,
`SENSITIVITY_CHECKPOINT_VALIDATION.json`, and `../environment/restore-validation.json`.
The `full` command must finish its comparisons before its output is marked PASS.
Full-scale or cross-platform discrepancies remain possible; a mismatch stops execution.

Publication inventory is `PUBLISH_FILES.csv`. The candidate excludes manuscripts, reviewer
history, raw empirical RDS/GeoJSON, private environments, temporary runs and credentials.

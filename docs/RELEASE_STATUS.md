# Validation history and current scope

Repository: `geowhk/ml-sas-paper`. DGP publication permission was confirmed by
the maintainer on 2026-10-02. R package code is maintained separately in
[ml-sas-r](https://github.com/geowhk/ml-sas-r).

## Recorded validation

The checks below were performed on 30 September 2026. Saved-result reproduction
was checked again on the separated repository on 2 October 2026. These are dated
validation records, not a claim that every analysis has been rerun after each
documentation update. Machine-readable records are in [VALIDATION.json](VALIDATION.json).

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

See `VALIDATION.json`, `EMPIRICAL_RERUN.json`, `REPEATED_SAMPLE_RERUN.json`,
`SENSITIVITY_CHECKPOINT_VALIDATION.json`, and `../environment/restore-validation.json`.
## Current validation limits

- All 300 repeated datasets have not been recomputed during repository validation.
- The complete sensitivity grids and full workflow have not been validated end
  to end at paper scale. The full execution route is provided.
- The paper workflow has been checked on macOS, not Windows or Linux. Checks of
  the separate R package do not validate this research workflow.

The `full` command must finish its comparisons before its output is marked PASS.
Full-scale or cross-platform discrepancies remain possible; a mismatch stops execution.

## Distributed files

[PUBLISH_FILES.csv](PUBLISH_FILES.csv) records distributed file checksums. The
repository includes source code, selected reference results, documentation and
README figures. Municipal inputs are obtained from the pinned public `tmap`
package; local environments and generated runs are not distributed.

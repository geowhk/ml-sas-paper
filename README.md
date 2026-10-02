# ml-sas-paper

Reproduction code for **ML-SAS: A Machine-Learning-Based Spatial Autocorrelation
Statistic for Mixed-Type Data**, by **Woohyung Kim and Sang-Il Lee**.

**Status:** accepted for publication in *Geographical Analysis*. The DOI and
published article link will be added when available.

## About the paper

ML-SAS combines spatial weights with attribute similarities learned from
unsupervised random-forest decision paths. It defines global and local spatial
autocorrelation statistics for continuous, categorical and mixed-type data,
with permutation tests for inference. The paper examines simulated spatial
processes, forest settings and proximity definitions, and applies the method
to Dutch municipalities.

This repository contains the paper's simulation and empirical-analysis workflows,
including the DGP code shared with permission, selected reference results,
and scripts for regenerating tables and figures from the final v6 analysis.
For applying ML-SAS to your own data, see the separate
[ml-sas-r package](https://github.com/geowhk/ml-sas-r) and its
[usage tutorial](https://geowhk.github.io/ml-sas-r/getting-started.html).

## Getting started

The checked environment uses **R 4.5.2 and Python 3.13**. R dependencies are pinned
in `renv.lock`; Python dependencies are listed in
`environment/requirements-figures.txt`. The commands below use a macOS/Linux shell.
The paper workflow has been checked on macOS; the R package's cross-platform
checks do not cover this separate reproduction workflow.

```sh
git clone https://github.com/geowhk/ml-sas-paper.git
cd ml-sas-paper
Rscript scripts/setup.R
python3 -m venv .venv
.venv/bin/python -m pip install -r environment/requirements-figures.txt
```

See [installation requirements](docs/INSTALL.md) for compilers and spatial libraries.

## Reproduce the results

**Quick check using saved results:** regenerates final summaries, Appendix C
diagnostics and ten figures without refitting the models.

```sh
.venv/bin/python scripts/reproduce.py --mode saved --out runs/saved
```

**Full recomputation:** reruns models, inference, simulations and sensitivity
analyses before assembling and comparing the final results.

```sh
.venv/bin/python scripts/reproduce.py --mode full --out runs/full --workers 1 --parallelism 1
```

The original research used **AWS**. Full recomputation can take hours to days
and requires adequate memory and uninterrupted runtime; start with saved mode
on a personal computer. See [runtime and resource guidance](docs/RUNTIME.md).
Add `--resume` to the same full command to resume an interrupted run.

The [reproduction guide](docs/REPRODUCING.md) gives stage details, output locations
and mappings to the paper's tables and figures.

## Validation

Saved-result reproduction, the full empirical analysis and one repeated-simulation
replication have been checked against the final v6 results. The entire 300-dataset
simulation batch and complete sensitivity grids have **not** been rerun during
this release preparation. See [validation scope](docs/RELEASE_STATUS.md).

## Authors and license

Woohyung Kim and Sang-Il Lee. Contact: whk3423@snu.ac.kr.

Project code is available under the [MIT license](LICENSE). Upstream dependencies
and data retain their own terms; see [sources and attribution](docs/THIRD_PARTY.md).

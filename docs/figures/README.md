# Figures in the README

These PNGs were regenerated from the repository's saved results on 2 October
2026. They correspond to Figures 7 and 9 in Kim and Lee (2026),
https://doi.org/10.1111/gean.70062. They are analysis-rendered images, not extracts
of the publisher's typeset PDF. Fonts and layout can differ across environments.

After installing the dependencies, run from the repository root:

```sh
.venv/bin/python scripts/reproduce.py --mode saved --out runs/saved
```

| README image | Generated output | Meaning |
|---|---|---|
| `figure-07.png` | `runs/saved/dgp-figures/Figure_07_dgp4_inputs_results.png` | Representative DGP 4 inputs and local detections |
| `figure-09.png` | `runs/saved/figure9/Figure_08_NLD_muni_local_clusters.png` | Continuous-to-mixed municipal classification transitions |

The Figure 9 renderer retains the historical filename `Figure_08`; the `figure9/`
output directory maps to Figure 9 in the paper. Use the DGP output above for
Figure 7, rather than the older plot emitted under `section4/`.

Figure 9 uses Statistics Netherlands/CBS municipal data and boundaries supplied
through `tmap::NLD_muni`; see [sources and rights](../THIRD_PARTY.md).

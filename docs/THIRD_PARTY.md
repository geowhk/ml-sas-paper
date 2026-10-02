# Sources and rights

The supplied DGP source and adaptations are included with permission, confirmed
by the maintainer on 2026-10-02. Project code is distributed under the MIT license
in `../LICENSE`; upstream dependencies and data retain their own terms.

The supplied original file remains byte-for-byte preserved at
`analysis/repeated/source_snapshot/DGP_generation.R`, SHA-256
`8599bf88fca254343c62486b5840597b8c5a5933e71c0ffeeec484b83fe0e90c`.
The relative frozen-source paths and source hashes are in `SOURCE_MANIFEST.csv`.

R/Python dependencies are not vendored. Direct R license metadata is recorded in
`environment/R-packages.csv`; consult upstream packages for complete terms.

Empirical inputs are loaded from tmap 4.2's `NLD_muni` (Statistics Netherlands/CBS,
2022 demographic data). The paper uses the largest queen-connected component:
337 of the package's 345 municipalities. Original geometry and demographic attributes
are not distributed here. Fingerprints verify them after acquisition through the package.
Derived ML-SAS result tables are included for saved-result reproduction.

Source/version references:
- https://r-tmap.github.io/tmap/articles/22_examples_choro_NLD.html
- https://www.cbs.nl/nl-nl/maatwerk/2024/11/kerncijfers-wijken-en-buurten-2022

Older tmap map documentation names both CBS and Kadaster and requires attribution:
https://r-tmap.github.io/tmap/reference/Shapes.html . Older and newer NLD_muni datasets
must not be substituted for each other. A fingerprint mismatch stops preparation.

Publication PDFs/TIFFs were used locally for visual comparison only and are not included.
Manuscripts, responses to reviewers, internal history, credentials and private source
snapshots are excluded. The paper's DOI and publication link will be added when available.

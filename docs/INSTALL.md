# Installation and environment checks

R 4.5.2 and Python 3.13 are the tested runtime versions. Install using the root README.
The lockfile pins R packages and their recursive dependencies; it does not install R,
compilers, operating-system libraries or fonts. Packages absent from current CRAN binary
repositories may be compiled from their archived source versions.

## macOS

Use the Command Line Tools and the GNU Fortran matching your CRAN R installation.
The authoritative toolchain instructions are https://mac.r-project.org/tools/ .
Do not replace R's toolchain globally as part of this project setup.

During the isolated-library test on this host, C++ source builds initially failed with
`cmath file not found`, and the R-configured `/opt/gfortran/bin/gfortran` was absent.
The C++ headers existed under the active SDK. A command-scoped diagnostic/workaround was:

```sh
CPLUS_INCLUDE_PATH="$(xcrun --show-sdk-path)/usr/include/c++/v1" Rscript scripts/setup.R
```

Use this only when that directory exists and the error matches. The local test also used
an explicit project-only Makevars pointing to the host's existing GNU Fortran 14.2 and
its libraries. This is a recorded local toolchain choice, not a portable setting to copy:
prefer the matching R-provided Fortran distribution. No user-wide Makevars was edited.

## Linux

A working C/C++/Fortran toolchain, CMake, and development libraries for GDAL, GEOS, PROJ,
udunits2, curl, OpenSSL, XML, fonts and text rendering are needed if R packages build
from source. Package names vary by distribution. Linux clean installation has not been
validated for this repository; do not interpret the macOS checks as Linux certification.

## What setup does

`scripts/setup.R` loads an isolated renv project, restores the lockfile with `retry=FALSE`,
and activates it. It does not silently replace a pinned version with a newer package.
See https://rstudio.github.io/renv/reference/restore.html .

Local renv libraries, caches, generated runs and virtual environments are excluded
from version control by `.gitignore`. `docs/PUBLISH_FILES.csv` records checksums
for distributed files (excluding the inventory itself).

If an input fingerprint, numeric comparison or package installation fails, retain the
log and fix the cause. Do not loosen numerical tolerances or bypass fingerprints just
to get PASS. Rendered image bytes may differ even when numeric inputs match.

The isolated test ultimately restored 169 packages after additionally supplying gettext
headers and building Unidata udunits 2.2.28 into a work-local prefix. For an existing custom
udunits installation, setup accepts `MLSAS_UDUNITS_PREFIX=/absolute/prefix`; this supplies
its `include/` and `lib/` directories to the units package configure step. This option is
not needed when the normal system development package is installed. The upstream source
is https://downloads.unidata.ucar.edu/udunits/ . A missing `libintl.h` refers to gettext
headers; it is an installation requirement, not a reason to change the analysis packages.

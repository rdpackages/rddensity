# Changelog

Notable project changes are listed from newest to oldest.

## 2026-05-21 - Python and R Modernization

- Reworked Python `rddensity`/`rdbwdensity` data preparation to use NumPy arrays internally while preserving the public Pandas result objects and numerical baseline.
- Added `scripts/benchmark-runtime.py` for repeatable Python, R, and optional Stata runtime benchmarks.
- Hardened the Stata benchmark runner for multi-size runs and missing-output diagnostics.
- Aligned Python package metadata and documentation version with the R/Stata `2.6` release line.
- Added cached kernel moment matrices in Python and R, matching the new `lpdensity` implementation style and avoiding repeated numerical integration.
- Reworked Python's core `rddensity` estimator to use NumPy arrays in the inner linear algebra routine.
- Fixed Python estimated-bandwidth indexing for current Pandas releases and added a regression test for that path.
- Removed Stata `preserve`/`restore` from `rddensityEST` and `rdbwdensity`; both now operate on typed temporary variables and sort inside Mata.

## 2026-05-21 - Stata Precision Baselines

- Added Stata option `precision(single|double)` to `rddensity` and `rdbwdensity`; `double` is now the default, while `single` preserves a legacy float-storage path.
- Routed Stata bandwidth selection, binomial testing, plotting helper variables, and generated variables through the selected precision mode.
- Extended the Stata numerical baseline to record both single- and double-precision results.
- Added a cross-language numerical comparison script for Python, R, Stata single precision, and Stata double precision.
- Regenerated Stata help PDFs documenting the new precision option.

## 2026-05-21 - Repository Setup

- Modernized repository infrastructure following the RDROBUST and LPDENSITY setup, including GitHub Actions CI, Dependabot configuration, issue templates, pull request template, security policy, and Python PyPI publishing through GitHub trusted publishing.
- Moved the repository from `master` to `main`, updated Stata installation URLs, removed the stale remote `master` branch, and added branch-protection settings for the `main` branch.
- Reorganized the top-level and Python READMEs, refreshed package metadata, and standardized the GPL-3 license notice.
- Updated maintainer and contributor metadata across Python, R, and Stata package files, keeping Matias D. Cattaneo as maintainer and Ricardo Masini as a Python contributor.
- Added Python, R, and Stata numerical baseline fixtures/scripts to protect current numerical behavior during future refactors.
- Added Stata help PDF generation scripts and regenerated help PDFs for `rddensity` and `rdbwdensity`.
- Added local repository notes in `AGENTS.md`, ignored from version control.
- Removed generated Python/R build artifacts and cached files from the tracked public package surfaces.
- Fixed the CI repository layout check to tolerate whitespace in the Stata package index.

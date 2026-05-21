# Changelog

Notable project changes are listed from newest to oldest.

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

# Changelog

## In Progress: May 22, 2026

- Began audit of the 0.2.0 R/Stata update and new Python package.
- Updated repository checks and numerical-baseline tooling for the new
  R vignette/data layout and pyproject-based Python package.
- Fixed initial R packaging warnings found by `R CMD check`: ASCII-only
  warning text, missing `stats::coef` import, and documentation for
  `rdhte_contrast()`.
- Cleaned the Python README encoding artifacts and development install path.
- Fixed Stata baseline failures caused by matrix elements embedded directly in
  data expressions and by CR2/CR3 aliases being passed to `regress` through
  unsupported HC-with-cluster syntax.
- Verified the new Python package tests and local no-isolation build; noted
  that an isolated build needs network access to refresh build dependencies.
- Ran first-pass benchmarks against commit `d4a05c7`; automatic-bandwidth R
  examples are now roughly 5x-9x faster on the local Windows setup.
- Regenerated R and Stata `*-next.json` baselines and compared them with the
  saved pre-update baselines; numerical changes remain and are being audited.
- Switched the Stata help PDF build list to the native Stata conversion path
  for all four help files, including `rdhte_plot`.
- Added Stata `precision(double|single)` to `rdhte` and `rdbwhte`; the default
  generated-variable path is now double precision, with `single` retained for
  backward numerical compatibility.
- Extended the Stata numerical-baseline runner with an optional precision
  argument so double- and single-precision compatibility checks can be tracked.
- Expanded the R illustration with `rdbwhte`, `covs.eff`, plotting, and table
  export examples, and added a visible Python `rdhte_illustration.py` entrypoint.
- Updated the top-level README to point to the new Python illustration and
  the Stata `rdhte_plot` PDF help file.
- Moved the Python illustration entrypoint to `Python/rdhte_illustration.py`
  and updated wrapper paths/documentation.
- Updated Matias, Max, Filippo, and Rocio contact emails across R, Stata, and
  Python package metadata/help files, and aligned the Python package license
  notice with the root `LICENSE.md`.
- Switched R, Python, and Stata package references for the two unpublished
  RDHTE papers to canonical arXiv links while retaining rdpackages links for
  background references.
- Brought Stata help-file orientation into the R and Python documentation:
  companion-command notes, empirical examples, option/return-value summaries,
  and fuller reference lists.
- Removed the R vignette directory and stale R session/build debris; updated
  repository checks for the no-vignette package layout.
- Made the Python illustration self-contained, including local source-path
  resolution and a headless plotting backend for saved figures.
- Removed cross-platform/platform-specific prose from public R, Python, and
  Stata help text while keeping Stata-only `precision()` documented in Stata.
- Added Python numerical-baseline generation and Stata cross-platform
  comparison tooling; current R/Python results agree at `1e-5`, and
  R/Stata/Python agree at `2e-5` on the shared illustration baseline fields.
- Aligned the Stata illustration with the cluster default variance path used
  by the R and Python illustrations.
- Refactored the R and Python illustrations into simple sequential
  replication scripts that mirror the Stata illustration: load the package,
  read the bundled data, and run the examples directly.
- Cleaned ignored verification debris (logs, caches, build artifacts,
  `R CMD check` output, generated plot/table files, and local baseline
  outputs) while keeping local-only `AGENTS.md`.
- Tightened the PyPI publishing workflow by separating distribution build from
  the OIDC-powered publish job and binding publication to the `pypi`
  environment.

## Modernization Summary: May 22, 2026

- Reorganized the public README to follow the RDROBUST multi-language layout, including Python, R, Stata, references, replication data, and funding sections.
- Added GitHub-facing project infrastructure: CI for repository layout and R checks, Python package checks that activate when `Python/rdhte/` is added, Python publishing to PyPI, Dependabot configuration, issue and pull request templates, and a security reporting policy.
- Replaced the placeholder license with the RDROBUST-style GPL-3.0 license notice adapted for `rdhte`.
- Added local and repository guardrails through `.gitignore`, `.gitattributes`, and local-only `AGENTS.md` notes.
- Added Stata help PDF generation scripts for `rdhte`, `rdbwhte`, and `rdhte_lincom`.
- Added numerical baseline scripts for R and Stata so future upgrades can compare results against working-tree reference outputs.
- Fixed Stata command compatibility issues uncovered by baseline generation: bandwidth matrix elements are copied before use in data expressions, optional efficiency covariate terms are omitted when empty, polynomial factor-variable syntax is normalized, and clustered variance options are routed to `regress` as `vce(cluster ...)`.

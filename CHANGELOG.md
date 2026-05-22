# Changelog

## Modernization Summary: May 22, 2026

- Reorganized the public README to follow the RDROBUST multi-language layout, including Python, R, Stata, references, replication data, and funding sections.
- Added GitHub-facing project infrastructure: CI for repository layout and R checks, Python package checks that activate when `Python/rdhte/` is added, Python publishing to PyPI, Dependabot configuration, issue and pull request templates, and a security reporting policy.
- Replaced the placeholder license with the RDROBUST-style GPL-3.0 license notice adapted for `rdhte`.
- Added local and repository guardrails through `.gitignore`, `.gitattributes`, and local-only `AGENTS.md` notes.
- Added Stata help PDF generation scripts for `rdhte`, `rdbwhte`, and `rdhte_lincom`.
- Added numerical baseline scripts for R and Stata so future upgrades can compare results against working-tree reference outputs.
- Fixed Stata command compatibility issues uncovered by baseline generation: bandwidth matrix elements are copied before use in data expressions, optional efficiency covariate terms are omitted when empty, polynomial factor-variable syntax is normalized, and clustered variance options are routed to `regress` as `vce(cluster ...)`.

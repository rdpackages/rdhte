# rdhte Numerical Baselines

Local benchmark harness for checking numerical stability across the R and Stata
implementations. Add the Python implementation here after `Python/rdhte/` lands.

The scripts write JSON files under `docs/audit/baselines/` by default. That
directory is local-only under the top-level `.gitignore`.

## Cases

- `binary_left`: illustration data with `w_left` as a binary subgroup variable.
- `binary_left_joint`: same subgroup variable with a common bandwidth.
- `binary_left_bw`: bandwidth selection for the binary subgroup case.
- `categorical_ideology`: illustration data with unordered ideology categories.
- `continuous_strength`: illustration data with `w_strength` as a continuous heterogeneity variable.
- `interaction_strength`: fixed-bandwidth binary-by-continuous interaction.
- `average_manual`: fixed-bandwidth average RD treatment effect, used to compare against `rdrobust`.

## Commands

From the repository root:

```powershell
Rscript scripts/numerical-baselines/run_r_baseline.R
```

If Stata is available:

```powershell
& "C:\Program Files\StataNow19\StataMP-64.exe" /e do scripts/numerical-baselines/run_stata_baseline.do
```

Then compare two JSON files:

```powershell
python scripts/numerical-baselines/compare_baselines.py docs/audit/baselines/r-current.json docs/audit/baselines/r-next.json
```

Use tight tolerances for baseline-to-baseline checks within the same language.
Use looser tolerances for cross-language checks because object structures and
printed precision differ between R and Stata.

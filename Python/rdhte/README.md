# rdhte (Python)

Heterogeneous treatment effects in regression discontinuity (RD) designs.
Point estimates, robust bias-corrected confidence intervals, and
heteroskedastic or cluster-robust standard errors for conditional RD
treatment effects.

The theoretical framework is developed in Calonico, Cattaneo, Farrell,
Palomba, and Titiunik (2025).
[arXiv](https://arxiv.org/abs/2503.13696)

## Install (development)

```bash
cd Python/rdhte
pip install -e ".[plot,test]"
```

Hard dependencies: `numpy`, `scipy`, `pandas`, `statsmodels`, `rdrobust`.
Plotting (`rdhte.plot.plot`) additionally requires `plotnine`.

## Quick start

```python
import numpy as np
from rdhte import rdhte, rdbwhte, rdhte_lincom

rng = np.random.default_rng(0)
n = 5000
x = rng.uniform(-1.0, 1.0, n)
W = np.where(rng.random(n) < 0.5, "A", "B")              # categorical W
y = ((x >= 0) * (1.0 + 0.5 * (W == "B"))                 # CATE = 1.0 (A), 1.5 (B)
     + 0.3 * x + rng.normal(0, 0.3, n))

# Per-cell CATEs with rdrobust-style bandwidth selection (default: mserd).
m = rdhte(y=y, x=x, c=0.0, covs_hte=W)
print(m)                 # short header
print(m.summary())       # full table with Nh/h columns
print(m.tidy())          # broom-style DataFrame: term, estimate, std.error, ...
print(m.glance())        # one-row metadata DataFrame

# Data-driven bandwidth per cell, without running the full estimator.
bw = rdbwhte(y=y, x=x, c=0.0, covs_hte=W)
print(bw.h)              # (n_lev, 2) array of left/right bandwidths

# Test H0: CATE(A) == CATE(B), and separately H0: CATE(A) == 0.
res = rdhte_lincom(m, linfct=["A - B = 0", "A = 0"])
print(res["individual"])
print(res["joint"])
```

## Empirical illustration

The package includes a full empirical illustration using the bundled
Granzier, Pons, and Tricaud election dataset.
From the `Python` folder:

```bash
python rdhte_illustration.py
```

## Command Reference

The companion commands are:

| Command | Purpose |
|---------|---------|
| `rdhte()` | Estimate conditional RD treatment effects and robust bias-corrected inference. |
| `rdbwhte()` | Compute data-driven left/right bandwidths without running the full estimator. |
| `rdhte_lincom()` | Test linear restrictions of the fitted CATE vector. |

Main estimation options:

| Option | Notes |
|--------|-------|
| `covs_hte` | Heterogeneity covariates. Accepts no covariate, one categorical/continuous variable, a list/DataFrame/2-D array for Cartesian-product cells, or a formula string with `data=`. |
| `c`, `p`, `q` | RD cutoff and local-polynomial orders. Default `p=1`; `q=None` resolves to `p + 1`. |
| `h`, `h_l`, `h_r` | Manual common or side-specific bandwidths. If omitted, `rdrobust.rdbwselect` is called using `bwselect`. |
| `kernel` | `triangular`, `epanechnikov`, or `uniform` (`tri`, `epa`, `uni` also accepted). |
| `weights` | Non-negative unit weights multiplying the kernel weights. |
| `covs_eff` | Efficiency-improving covariates, included additively and through W interactions; they do not define CATE cells. |
| `vce`, `cluster` | `hc0`-`hc3` without clusters; `cr1`-`cr3` with `cluster=`. Legacy `hc*` plus `cluster=` is remapped to the matching `cr*` with a warning. |
| `level` | Confidence level in percent, e.g. `95`, not `0.95`. |
| `bwjoint` | For categorical `covs_hte`, force one shared bandwidth across groups. Continuous/formula paths already use a shared bandwidth. |
| `data`, `subset` | Resolve column-name strings against a DataFrame and optionally keep a subset of rows. |

`bwselect` accepts `mserd`, `msetwo`, `msesum`, `msecomb1`, `msecomb2`,
`cerrd`, `certwo`, `cersum`, `cercomb1`, and `cercomb2`. MSE denotes
mean square error; CER denotes coverage error rate.

## Return Objects

`rdhte()` returns an `RdhteResult` with fields for estimates, inference
results, bandwidths, labels, and metadata:
`Estimate`, `Estimate_bc`, `coef`, `coef_bc`, `se_rb`, `ci_rb`, `t_rb`,
`pv_rb`, `vcov`, `W_lev`, `W_names`, `kernel`, `bwselect`, `vce`,
`vce_select`, `c`, `h`, `p`, `q`, `N`, `Nh`, `level`, and `rdmodel`.
Convenience methods include `summary()`, `tidy()`, and `glance()`.

`rdbwhte()` returns an `RdbwhteResult` with bandwidths `h`, effective
sample sizes `Nh`, group labels, and the same metadata. `rdhte_lincom()`
returns a dictionary with `individual` and `joint` test-result DataFrames.

### Multiple heterogeneity variables

For factor-by-factor interactions and formula expressions, supply
`covs_hte` as a list of arrays / DataFrame / formula string:

```python
import pandas as pd
df = pd.DataFrame({"y": y, "x": x, "W1": W1, "W2": W2})

# Cartesian-product factor cells.
m_2d = rdhte(y="y", x="x", c=0.0, covs_hte=[df["W1"], df["W2"]], data=df)

# Formula syntax (linear-in-Wk CATE per term).
m_form = rdhte(y="y", x="x", c=0.0, covs_hte="W1 * W2", data=df)
# Coefs: T, T:W1, T:W2, T:W1.W2.

# factor() is accepted as an alias for patsy C().
m_mix = rdhte(y="y", x="x", c=0.0, covs_hte="factor(W1) * W2", data=df)
```

## Supported Features

| Feature                                          | Status |
|--------------------------------------------------|--------|
| Categorical `covs_hte` (one variable)            | yes |
| Continuous `covs_hte` (single numeric column)    | yes |
| No `covs_hte` (ATE-only path)                    | yes |
| Manual bandwidth (`h=`) or per-side `h_l`/`h_r`  | yes |
| `rdrobust.rdbwselect` automatic bandwidth        | yes |
| Cluster-robust SE (`cluster=`, `vce="cr1"`/`cr2`/`cr3`) | yes |
| Heteroskedastic SE (`vce="hc0"`/`hc1`/`hc2`/`hc3`) | yes |
| Robust bias correction (q-order parallel fit)    | yes |
| Linear-combination tests (`rdhte_lincom`)        | yes |
| Forest plot for categorical W (`rdhte.plot.plot`) | yes |
| `summary()` / `tidy()` / `glance()` methods      | yes |
| Factor-by-factor interaction `covs_hte`          | yes (pass `[W1, W2]` / `pd.DataFrame` / 2-D ndarray) |
| Formula-string `covs_hte` (e.g. `"w_left*w_strength"`) | yes (requires `data=`; `factor(X)` is accepted as an alias for patsy `C(X)`) |
| `covs_eff` (efficiency-improving covariates)     | yes (1-D / 2-D / DataFrame; enters additively + W interactions; never with T or Xp) |

`vce` accepts: `nn`, `hc0`, `hc1`, `hc2`, `hc3`, `cr1`, `cr2`, `cr3`.
When `cluster=` is supplied the default switches to `cr1`; passing an
`hc*` value with a cluster issues a warning and upgrades to the
matching `cr*` form.

## References

- Calonico, Cattaneo, Farrell, Palomba, and Titiunik (2025a):
  [Treatment Effect Heterogeneity in Regression Discontinuity Designs](https://arxiv.org/abs/2503.13696).
- Calonico, Cattaneo, Farrell, Palomba, and Titiunik (2025b):
  [rdhte: Conditional Average Treatment Effects in RD Designs](https://arxiv.org/abs/2507.01128).
- Cattaneo and Titiunik (2022):
  [Regression Discontinuity Designs](https://rdpackages.github.io/references/Cattaneo-Titiunik_2022_ARE.pdf).
- Calonico, Cattaneo, and Farrell (2020):
  [Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ.pdf).
- Calonico, Cattaneo, Farrell, and Titiunik (2019):
  [Regression Discontinuity Designs using Covariates](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT.pdf).
- Calonico, Cattaneo, and Titiunik (2014):
  [Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf).
- Granzier, Pons, and Tricaud (2023):
  [Coordination and Bandwagon Effects: How Past Rankings Shape the Behavior of Voters and Candidates](https://www.aeaweb.org/articles?id=10.1257/app.20210840).

Related RD software is available at <https://rdpackages.github.io/>.

## Authors

Sebastian Calonico, Matias D. Cattaneo, Max H. Farrell, Filippo Palomba,
Rocio Titiunik.

## License

GPL-3; see `LICENSE` in this package and `LICENSE.md` at the repository root.

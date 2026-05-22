"""Phase-8 tests: ``covs_eff`` (efficiency-improving covariates).

Mirrors R rdhte ``covs.eff = ...``. The covariates enter the regression
additively and (in W-aware paths) also interact with W, but never with
T or the running-variable polynomial. They tighten the SE of the CATE
without changing its identification.

Numerical parity vs R on the bundled `rdhte_dataset.csv` fixture.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd
import pytest

from rdhte import rdhte


HERE = os.path.dirname(__file__)


def _load_rdhte_dataset():
    fixture = os.path.join(HERE, "..", "..", "..", "R", "rdhte_dataset.csv")
    if not os.path.exists(fixture):
        pytest.skip(f"missing fixture: {fixture}")
    return pd.read_csv(fixture)


# ----------------------------- API tests ----------------------------------

@pytest.fixture(scope="module")
def synth_w_eff():
    rng = np.random.default_rng(0)
    n = 1500
    x = rng.uniform(-1, 1, n)
    W = np.where(rng.random(n) < 0.5, "A", "B")
    z1 = rng.normal(0, 1, n)
    z2 = rng.normal(0, 1, n)
    y = (
        (x >= 0).astype(float) * (1.0 + 0.5 * (W == "B"))
        + 0.3 * x + 0.4 * z1 - 0.2 * z2 + rng.normal(0, 0.3, n)
    )
    return dict(y=y, x=x, W=W, z1=z1, z2=z2)


def test_covs_eff_categorical_tightens_se(synth_w_eff):
    """SE should drop when adding genuinely-relevant covs_eff."""
    d = synth_w_eff
    m0 = rdhte(y=d["y"], x=d["x"], c=0.0, covs_hte=d["W"], h=0.3)
    m1 = rdhte(y=d["y"], x=d["x"], c=0.0, covs_hte=d["W"],
               covs_eff=np.column_stack([d["z1"], d["z2"]]), h=0.3)
    assert (m1.se_rb < m0.se_rb).all(), \
        f"covs_eff should shrink SEs; got {m0.se_rb} -> {m1.se_rb}"


def test_covs_eff_continuous_path(synth_w_eff):
    """Continuous-W + covs_eff runs and returns 2 coefs (T, slope)."""
    d = synth_w_eff
    m = rdhte(y=d["y"], x=d["x"], c=0.0,
              covs_hte=d["z1"],  # continuous W
              covs_eff=d["z2"], h=0.3)
    assert list(m.coef.index) == ["T", "T:covs.hte"]


def test_covs_eff_ate_only_path(synth_w_eff):
    """ATE-only + covs_eff: T coef present, SE finite."""
    d = synth_w_eff
    m = rdhte(y=d["y"], x=d["x"], c=0.0,
              covs_eff=np.column_stack([d["z1"], d["z2"]]), h=0.3)
    assert "T" in m.coef.index
    assert np.isfinite(m.se_rb[0])


def test_covs_eff_formula_path(synth_w_eff):
    """Formula-string covs_hte + covs_eff."""
    d = synth_w_eff
    df = pd.DataFrame(d)
    m = rdhte(y="y", x="x", c=0.0,
              covs_hte="W", covs_eff=df[["z1", "z2"]], h=0.3, data=df)
    # Formula path treats W (object dtype) as a categorical via patsy → C(W).
    assert m.coef is not None


def test_covs_eff_dataframe_input(synth_w_eff):
    """DataFrame input for covs_eff should match a 2-D ndarray of the same data."""
    d = synth_w_eff
    Z_df = pd.DataFrame({"z1": d["z1"], "z2": d["z2"]})
    m_df = rdhte(y=d["y"], x=d["x"], c=0.0, covs_hte=d["W"],
                 covs_eff=Z_df, h=0.3)
    m_arr = rdhte(y=d["y"], x=d["x"], c=0.0, covs_hte=d["W"],
                  covs_eff=Z_df.values, h=0.3)
    np.testing.assert_allclose(m_df.coef.values, m_arr.coef.values, atol=1e-12)


def test_covs_eff_length_mismatch_errors(synth_w_eff):
    d = synth_w_eff
    bad = d["z1"][:100]  # wrong length
    with pytest.raises(ValueError, match="covs_eff"):
        rdhte(y=d["y"], x=d["x"], c=0.0, covs_hte=d["W"],
              covs_eff=bad, h=0.3)


# ----------------------------- R-parity tests -----------------------------

# Frozen R reference values, harvested with:
#   library(rdhte)
#   d <- read.csv("rdhte_dataset.csv", stringsAsFactors = TRUE)
#   m  <- rdhte(y = y, x = x, covs.hte = factor(w_left), covs.eff = w_strength,
#               cluster = cluster_var, h = 0.1, data = d)
#   m2 <- rdhte(y = y, x = x, covs.hte = w_strength, covs.eff = w_ideology,
#               cluster = cluster_var, h = 0.1, data = d)

R_REF_CATEGORICAL = {
    "T":         0.02274,      # CATE for w_left = 0  (R reference category)
    "tau_1":     0.08780,      # CATE for w_left = 1
    "se_T":      0.00685,
    "se_1":      0.01367,
}

R_REF_CONTINUOUS = {
    "T":          -0.06121,
    "T:covs.hte":  0.27658,
    "se_T":        0.03024,
    "se_slope":    0.08017,
}

TOL = 1e-3


def test_r_parity_categorical_with_covs_eff():
    d = _load_rdhte_dataset()
    m = rdhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
              covs_hte=d["w_left"].to_numpy(),
              covs_eff=d["w_strength"].to_numpy(),
              cluster=d["cluster_var"].to_numpy(), h=0.1)
    # Python's saturated-categorical path returns one row per W level by name.
    assert sorted(m.W_lev) == ["0", "1"]
    tau_0 = float(m.coef["0"])
    tau_1 = float(m.coef["1"])
    assert abs(tau_0 - R_REF_CATEGORICAL["T"])     < TOL, f"tau(0): {tau_0} vs R {R_REF_CATEGORICAL['T']}"
    assert abs(tau_1 - R_REF_CATEGORICAL["tau_1"]) < TOL, f"tau(1): {tau_1} vs R {R_REF_CATEGORICAL['tau_1']}"


def test_r_parity_continuous_with_covs_eff():
    d = _load_rdhte_dataset()
    m = rdhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
              covs_hte=d["w_strength"].to_numpy(),
              covs_eff=d["w_ideology"].to_numpy(),
              cluster=d["cluster_var"].to_numpy(), h=0.1)
    assert abs(float(m.coef["T"])          - R_REF_CONTINUOUS["T"])         < TOL
    assert abs(float(m.coef["T:covs.hte"]) - R_REF_CONTINUOUS["T:covs.hte"]) < TOL
    # SEs should also agree within tolerance.
    se_T     = float(m.se_rb[0])
    se_slope = float(m.se_rb[1])
    assert abs(se_T     - R_REF_CONTINUOUS["se_T"])     < TOL
    assert abs(se_slope - R_REF_CONTINUOUS["se_slope"]) < TOL

"""Phase-7 tests: formula-string `covs_hte`.

Mirrors R rdhte cases 4 (additive) and 5 (full interaction):

    rdhte(..., covs.hte = "W2 + W3")
    rdhte(..., covs.hte = "w_left*w_strength")
    rdhte(..., covs.hte = "factor(W1)*W3")

R-style `factor(X)` is auto-translated to patsy's `C(X)` so users can
paste R formulas verbatim.

Numerical parity against R (rdhte 0.2.0 + rdrobust 4.0.0) is verified
on the bundled `rdhte_dataset.csv` fixture.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd
import pytest

from rdhte import rdhte


HERE = os.path.dirname(__file__)


@pytest.fixture(scope="module")
def synth_two_var():
    rng = np.random.default_rng(0)
    n = 1500
    df = pd.DataFrame({
        "x":          rng.uniform(-1, 1, n),
        "w_left":     rng.binomial(1, 0.5, n).astype(float),
        "w_strength": rng.uniform(0, 1, n),
    })
    df["y"] = (
        (df["x"] >= 0).astype(float)
        * (0.2 + 0.5 * df["w_left"] + 0.3 * df["w_strength"]
           + 0.4 * df["w_left"] * df["w_strength"])
        + 0.3 * df["x"] + rng.normal(0, 0.3, n)
    )
    return df


# -------------------------- shape / API tests ----------------------------

def test_formula_additive_shape(synth_two_var):
    df = synth_two_var
    m = rdhte(y="y", x="x", c=0.0, covs_hte="w_left + w_strength", h=0.3, data=df)
    assert list(m.coef.index) == ["T", "T:w_left", "T:w_strength"]
    assert m.coef.shape == (3,)
    assert (m.se_rb > 0).all()


def test_formula_full_interaction_shape(synth_two_var):
    df = synth_two_var
    m = rdhte(y="y", x="x", c=0.0, covs_hte="w_left*w_strength", h=0.3, data=df)
    # R names: T, T:w_left, T:w_strength, T:w_left:w_strength
    # Python uses '.' in the interaction term so rdhte_lincom backticks work.
    assert m.coef.shape == (4,)
    assert "T" in m.coef.index
    assert "T:w_left" in m.coef.index
    assert "T:w_strength" in m.coef.index
    assert "T:w_left.w_strength" in m.coef.index


def test_formula_factor_translation(synth_two_var):
    """R's `factor(X)` auto-translates to patsy's `C(X)`."""
    df = synth_two_var
    m = rdhte(y="y", x="x", c=0.0, covs_hte="factor(w_left)*w_strength",
              h=0.3, data=df)
    # Just verify it doesn't crash; column naming differs from R for factors
    # because patsy uses `C(X)[T.lvl]` style.
    assert m.coef is not None
    assert (m.se_rb > 0).all()


def test_formula_requires_data():
    with pytest.raises(ValueError, match="data="):
        rdhte(y=np.zeros(10), x=np.zeros(10), c=0.0,
              covs_hte="w_left + w_strength", h=0.3, data=None)


# -------------------------- R-parity tests --------------------------------

# Frozen R reference values, harvested by:
#   library(rdhte)
#   d <- read.csv("rdhte_dataset.csv", stringsAsFactors = TRUE)
#   m4 <- rdhte(y = y, x = x, covs.hte = "w_left + w_strength",
#               cluster = cluster_var, data = d)
#   m5 <- rdhte(y = y, x = x, covs.hte = "w_left*w_strength",
#               cluster = cluster_var, data = d)

R_REF_CASE4 = {
    "T":            -0.05151712,
    "T:w_left":      0.05731526,
    "T:w_strength":  0.18956900,
    "h":             0.09654631,
}

R_REF_CASE5 = {
    "T":                     -0.03251453,
    "T:w_left":              -0.05903098,
    "T:w_strength":           0.14085000,
    "T:w_left.w_strength":    0.27448080,
    "h":                      0.09654631,
}

TOL = 1e-3   # rdrobust internal h-selection drift gives 4-decimal agreement here


def _load_rdhte_dataset():
    fixture = os.path.join(HERE, "..", "..", "..", "R", "rdhte_dataset.csv")
    if not os.path.exists(fixture):
        pytest.skip(f"missing fixture: {fixture}")
    return pd.read_csv(fixture)


def test_r_parity_additive_formula():
    df = _load_rdhte_dataset()
    m = rdhte(y="y", x="x", covs_hte="w_left + w_strength",
              cluster="cluster_var", data=df)
    for name, ref in R_REF_CASE4.items():
        if name == "h":
            assert abs(m.h[0, 0] - ref) < TOL, f"h: {m.h[0,0]} vs R {ref}"
            continue
        got = float(m.coef[name])
        assert abs(got - ref) < TOL, f"{name}: Py {got} vs R {ref}"


def test_r_parity_full_interaction_formula():
    df = _load_rdhte_dataset()
    m = rdhte(y="y", x="x", covs_hte="w_left*w_strength",
              cluster="cluster_var", data=df)
    for name, ref in R_REF_CASE5.items():
        if name == "h":
            assert abs(m.h[0, 0] - ref) < TOL, f"h: {m.h[0,0]} vs R {ref}"
            continue
        got = float(m.coef[name])
        assert abs(got - ref) < TOL, f"{name}: Py {got} vs R {ref}"

"""Phase-5 tests: __repr__, summary(), tidy(), glance() on RdhteResult / RdbwhteResult.

The R-side reference behavior is documented in
`packages/R/rdhte/R/rdhte.R::print.summary.rdhte` (and the matching
``tidy.rdhte`` / ``glance.rdhte``). We check column names and a small
set of value invariants rather than diffing the full ASCII table — the
table contents are already verified by the parity tests, this file
locks in the *shape* of the methods.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd
import pytest

from rdhte import rdhte, rdbwhte


HERE = os.path.dirname(__file__)


@pytest.fixture(scope="module")
def fit_categorical():
    rng = np.random.default_rng(0)
    n = 800
    x = rng.uniform(-1, 1, n)
    W = np.where(rng.random(n) < 0.5, "A", "B")
    y = ((x >= 0) * (1 + 0.5 * (W == "B")) + 0.3 * x + rng.normal(0, 0.3, n))
    return rdhte(y=y, x=x, c=0.0, covs_hte=W, h=0.3)


@pytest.fixture(scope="module")
def fit_continuous():
    rng = np.random.default_rng(1)
    n = 800
    x = rng.uniform(-1, 1, n)
    W = rng.uniform(0, 1, n)
    y = ((x >= 0) * (1 + 0.5 * W) + 0.3 * x + rng.normal(0, 0.3, n))
    return rdhte(y=y, x=x, c=0.0, covs_hte=W, h=0.3)


@pytest.fixture(scope="module")
def fit_bw_categorical():
    rng = np.random.default_rng(2)
    n = 800
    x = rng.uniform(-1, 1, n)
    W = np.where(rng.random(n) < 0.5, "A", "B")
    y = ((x >= 0) * (1 + 0.5 * (W == "B")) + 0.3 * x + rng.normal(0, 0.3, n))
    return rdbwhte(y=y, x=x, c=0.0, covs_hte=W)


# ----------------------------- RdhteResult ---------------------------------

def test_rdhte_repr_short_categorical(fit_categorical):
    out = repr(fit_categorical)
    # Header line is the rdhte tagline; should NOT contain the full table.
    assert "RD Heterogeneous Treatment Effects" in out
    assert "Number of Obs." in out
    assert "Bandwidth (h)" in out
    # No table separator.
    assert "===" not in out


def test_rdhte_summary_has_table_categorical(fit_categorical):
    s = fit_categorical.summary()
    # Top metadata.
    assert "Sharp RD" in s
    assert "BW type" in s
    # Estimate table column labels.
    for label in ("Point", "Robust Inference", "Estimate", "z", "Pr(>|z|)", "C.I.", "Nh-", "Nh+", "h-", "h+"):
        assert label in s, f"summary missing {label!r}"
    # Row labels for each cell.
    assert "A" in s and "B" in s
    # Table separators present.
    assert "===" in s and "---" in s


def test_rdhte_summary_continuous_skips_Nh_h_columns(fit_continuous):
    s = fit_continuous.summary()
    # Continuous case: per-row h/Nh columns are not in the table.
    assert "Nh-" not in s
    assert "h-" not in s
    # But BW est. (h) line IS present in the header block.
    assert "BW est. (h)" in s


def test_rdhte_tidy_columns_and_rows(fit_categorical):
    df = fit_categorical.tidy()
    expected = ["term", "estimate", "std.error", "statistic", "p.value",
                "conf.low", "conf.high", "estimate.bc",
                "h.left", "h.right", "n.eff.left", "n.eff.right"]
    assert list(df.columns) == expected
    assert list(df["term"]) == ["A", "B"]
    # estimate equals coef Series values.
    np.testing.assert_allclose(df["estimate"].values, fit_categorical.coef.values)
    # estimate.bc equals coef_bc.
    np.testing.assert_allclose(df["estimate.bc"].values, fit_categorical.coef_bc.values)
    # std.error positive.
    assert (df["std.error"] > 0).all()


def test_rdhte_tidy_continuous(fit_continuous):
    df = fit_continuous.tidy()
    # Continuous case: rows are "T" and "T:covs.hte".
    assert list(df["term"]) == ["T", "T:covs.hte"]


def test_rdhte_glance_one_row(fit_categorical):
    g = fit_categorical.glance()
    assert isinstance(g, pd.DataFrame)
    assert len(g) == 1
    for col in ("n", "cutoff", "p", "q", "kernel", "vce",
                "bwselect", "level", "n.terms", "covs.continuous"):
        assert col in g.columns
    assert g["n.terms"].iloc[0] == 2
    assert g["covs.continuous"].iloc[0] is False or g["covs.continuous"].iloc[0] == False  # noqa: E712


# ----------------------------- RdbwhteResult -------------------------------

def test_rdbwhte_repr_short(fit_bw_categorical):
    out = repr(fit_bw_categorical)
    assert "MSE-Optimal" in out
    assert "Number of Obs." in out
    assert "===" not in out


def test_rdbwhte_summary_table(fit_bw_categorical):
    s = fit_bw_categorical.summary()
    assert "BW type" in s
    assert "Group" in s
    assert "h-" in s and "h+" in s
    assert "===" in s


def test_rdbwhte_tidy(fit_bw_categorical):
    df = fit_bw_categorical.tidy()
    assert list(df.columns) == ["term", "h.left", "h.right"]
    # Per-cell rows for the two levels.
    assert sorted(df["term"]) == ["A", "B"]
    assert (df["h.left"] > 0).all()
    assert (df["h.right"] > 0).all()


def test_rdbwhte_glance(fit_bw_categorical):
    g = fit_bw_categorical.glance()
    assert len(g) == 1
    for col in ("n", "cutoff", "p", "q", "kernel", "vce", "bwselect", "n.terms", "covs.continuous"):
        assert col in g.columns

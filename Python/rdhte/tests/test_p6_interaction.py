"""Phase-6 tests: factor-by-factor interaction `covs_hte`.

Verifies the three accepted multi-factor input shapes (list of arrays,
pandas DataFrame, 2-D ndarray) produce identical results, and that the
cell labels follow R's ``interaction(..., sep=".")`` convention.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from rdhte import rdhte, rdbwhte


@pytest.fixture(scope="module")
def two_factor_sample():
    rng = np.random.default_rng(0)
    n = 1200
    x = rng.uniform(-1, 1, n)
    W1 = np.where(rng.random(n) < 0.5, "A", "B")
    W2 = np.where(rng.random(n) < 0.5, "lo", "hi")
    y = (
        (x >= 0) * (1.0 + 0.5 * (W1 == "B") + 0.3 * (W2 == "hi"))
        + 0.3 * x + rng.normal(0, 0.3, n)
    )
    return dict(y=y, x=x, W1=W1, W2=W2)


def test_interaction_list_of_arrays(two_factor_sample):
    d = two_factor_sample
    m = rdhte(y=d["y"], x=d["x"], c=0.0, covs_hte=[d["W1"], d["W2"]], h=0.3)
    assert sorted(m.W_lev) == ["A.hi", "A.lo", "B.hi", "B.lo"]
    assert m.coef.shape == (4,)
    assert (m.se_rb > 0).all()


def test_interaction_three_input_forms_match(two_factor_sample):
    d = two_factor_sample
    m_list = rdhte(y=d["y"], x=d["x"], c=0.0,
                   covs_hte=[d["W1"], d["W2"]], h=0.3)
    m_df = rdhte(y=d["y"], x=d["x"], c=0.0,
                 covs_hte=pd.DataFrame({"W1": d["W1"], "W2": d["W2"]}), h=0.3)
    m_arr = rdhte(y=d["y"], x=d["x"], c=0.0,
                  covs_hte=np.column_stack([d["W1"], d["W2"]]), h=0.3)
    # Reorder by W_lev so column-order differences (if any) don't matter.
    def sort(mm):
        idx = np.argsort(mm.W_lev)
        return np.asarray(mm.coef.values)[idx], np.asarray(mm.se_rb)[idx]
    c1, s1 = sort(m_list); c2, s2 = sort(m_df); c3, s3 = sort(m_arr)
    np.testing.assert_allclose(c1, c2, atol=1e-12)
    np.testing.assert_allclose(c1, c3, atol=1e-12)
    np.testing.assert_allclose(s1, s2, atol=1e-12)
    np.testing.assert_allclose(s1, s3, atol=1e-12)


def test_interaction_rejects_continuous_component(two_factor_sample):
    d = two_factor_sample
    W3 = np.random.default_rng(1).uniform(0, 1, len(d["y"]))
    with pytest.raises(ValueError, match="categorical"):
        rdhte(y=d["y"], x=d["x"], c=0.0, covs_hte=[d["W1"], W3], h=0.3)


def test_rdbwhte_interaction_per_cell_h(two_factor_sample):
    d = two_factor_sample
    bw = rdbwhte(y=d["y"], x=d["x"], c=0.0, covs_hte=[d["W1"], d["W2"]])
    assert sorted(bw.W_lev) == ["A.hi", "A.lo", "B.hi", "B.lo"]
    assert bw.h.shape == (4, 2)
    assert (bw.h > 0).all()


def test_interaction_drops_rows_with_any_factor_na(two_factor_sample):
    d = two_factor_sample
    W1 = pd.Series(d["W1"]).astype("object").copy()
    W1.iloc[:50] = None
    m_with_na = rdhte(y=d["y"], x=d["x"], c=0.0,
                       covs_hte=[W1.values, d["W2"]], h=0.3)
    # Drop the same rows manually to get the gold reference.
    keep = ~W1.isna().values
    m_clean = rdhte(y=d["y"][keep], x=d["x"][keep], c=0.0,
                    covs_hte=[d["W1"][keep], d["W2"][keep]], h=0.3)
    # NA-dropped fit should match the clean fit's coef and SEs exactly.
    np.testing.assert_allclose(m_with_na.coef.values, m_clean.coef.values, atol=1e-12)
    np.testing.assert_allclose(m_with_na.se_rb,        m_clean.se_rb,        atol=1e-12)


def test_interaction_tidy_has_one_row_per_cell(two_factor_sample):
    d = two_factor_sample
    m = rdhte(y=d["y"], x=d["x"], c=0.0, covs_hte=[d["W1"], d["W2"]], h=0.3)
    df = m.tidy()
    assert sorted(df["term"]) == ["A.hi", "A.lo", "B.hi", "B.lo"]
    assert len(df) == 4

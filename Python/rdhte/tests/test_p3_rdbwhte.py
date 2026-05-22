"""Phase-3 parity: rdbwhte vs R rdbwhte."""

from __future__ import annotations

import os
import numpy as np
import pandas as pd
import pytest


HERE = os.path.dirname(__file__)


def _load(filename):
    path = os.path.join(HERE, filename)
    if not os.path.exists(path):
        pytest.skip(f"missing fixture: {filename}")
    return pd.read_csv(path)


# R reference (from R rdbwhte on parity_fixture.csv):
R_BWHTE_PER_CELL = {
    "A": (0.3931628, 0.2701059),
    "B": (0.2415309, 0.3474842),
}
R_BWHTE_JOINT = (0.2399831, 0.3749588)

TOL = 1e-5


def test_rdbwhte_per_cell():
    from rdhte import rdbwhte
    d = _load("parity_fixture.csv")
    res = rdbwhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
                  c=0.0, covs_hte=d["W"].to_numpy(),
                  bwselect="msetwo")
    # Levels A, B should be in that order.
    assert res.W_lev == ["A", "B"]
    for j, lev in enumerate(res.W_lev):
        ref_l, ref_r = R_BWHTE_PER_CELL[lev]
        assert abs(res.h[j, 0] - ref_l) < TOL, f"{lev} h_l"
        assert abs(res.h[j, 1] - ref_r) < TOL, f"{lev} h_r"
    # Per-cell mode also: covs_cont False, returned shape correct.
    assert res.covs_cont is False
    assert res.h.shape == (2, 2)
    assert res.Nh.shape == (2, 2)


def test_rdbwhte_joint():
    from rdhte import rdbwhte
    d = _load("parity_fixture.csv")
    res = rdbwhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
                  c=0.0, covs_hte=d["W"].to_numpy(),
                  bwselect="msetwo", bwjoint=True)
    for j in range(len(res.W_lev)):
        assert abs(res.h[j, 0] - R_BWHTE_JOINT[0]) < TOL
        assert abs(res.h[j, 1] - R_BWHTE_JOINT[1]) < TOL


def test_rdbwhte_no_covs_hte():
    from rdhte import rdbwhte
    d = _load("parity_fixture.csv")
    res = rdbwhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
                  c=0.0, bwselect="msetwo")
    # ATE-only path: one bw pair, single-row matrix.
    assert res.h.shape == (1, 2)
    assert res.covs_cont is False
    assert (res.h > 0).all()


def test_rdbwhte_continuous_covs_hte():
    from rdhte import rdbwhte
    d = _load("parity_fixture_continuous.csv")
    res = rdbwhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
                  c=0.0, covs_hte=d["W"].to_numpy(),
                  bwselect="msetwo")
    # Continuous W -> one shared bw pair, covs_cont=True.
    assert res.h.shape == (1, 2)
    assert res.covs_cont is True
    assert (res.h > 0).all()

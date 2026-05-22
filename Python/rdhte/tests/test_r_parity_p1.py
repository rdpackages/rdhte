"""Phase-1 R parity tests.

Compares the Python rdhte() output against R rdhte::rdhte on a shared
CSV fixture. R reference values are frozen below (harvested by the
companion R script in ``_r_parity_fixture.R``). Tests pass with a
tolerance of 1e-5 on point estimates and SEs.
"""

from __future__ import annotations

import os
import numpy as np
import pandas as pd
import pytest


HERE = os.path.dirname(__file__)


def _load(filename):
    path = os.path.join(HERE, filename)
    if not os.path.exists(path):
        pytest.skip(f"missing fixture: {filename}; run tests/_make_parity_fixture.py")
    return pd.read_csv(path)


# --- Frozen R reference values ----------------------------------------------
# R script that produced these (see also `_r_parity_fixture.R`):
#
#   library(rdhte)
#   d <- read.csv("parity_fixture.csv")
#   m <- rdhte(y=d$y, x=d$x, c=0, covs.hte=factor(d$W), h=0.3)
#   m_auto <- rdhte(y=d$y, x=d$x, c=0, covs.hte=factor(d$W), bwselect="msetwo", bw.joint=TRUE)
#   d2 <- read.csv("parity_fixture_cluster.csv")
#   m_cl <- rdhte(y=d2$y, x=d2$x, c=0, covs.hte=factor(d2$W), cluster=d2$cl, h=0.3)
#
R_REF = {
    "manual_h": {
        "Estimate":     {"A": 1.082484, "B": 1.474088},
        "Estimate.bc":  {"A": 1.096741, "B": 1.561586},
        "se.rb":        {"A": 0.1300707, "B": 0.1270590},
    },
    "auto_bw_joint": {
        "Estimate":     {"A": 1.081221, "B": 1.477122},
        "Estimate.bc":  {"A": 1.101711, "B": 1.559012},
        "se.rb":        {"A": 0.1317766, "B": 0.1284379},
        "h":            [0.2399831, 0.3749588],
    },
    "cluster": {
        "Estimate":     {"A": 0.5706814, "B": 1.4086951},
        "Estimate.bc":  {"A": 0.3714099, "B": 1.2733193},
        "se.rb":        {"A": 0.1763464, "B": 0.1612806},
    },
}

TOL = 1e-5


def _check(m, ref, *, has_h=False):
    est    = m.coef
    est_bc = m.coef_bc
    se     = pd.Series(m.se_rb, index=m.W_lev)
    assert abs(est["A"]    - ref["Estimate"]["A"])    < TOL
    assert abs(est["B"]    - ref["Estimate"]["B"])    < TOL
    assert abs(est_bc["A"] - ref["Estimate.bc"]["A"]) < TOL
    assert abs(est_bc["B"] - ref["Estimate.bc"]["B"]) < TOL
    assert abs(se["A"]     - ref["se.rb"]["A"])       < TOL
    assert abs(se["B"]     - ref["se.rb"]["B"])       < TOL
    if has_h:
        h_l, h_r = m.h[0, 0], m.h[0, 1]
        assert abs(h_l - ref["h"][0]) < TOL
        assert abs(h_r - ref["h"][1]) < TOL


def test_parity_manual_h():
    from rdhte import rdhte
    d = _load("parity_fixture.csv")
    m = rdhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
              c=0.0, covs_hte=d["W"].to_numpy(), h=0.3)
    _check(m, R_REF["manual_h"])


def test_parity_auto_bw_joint():
    from rdhte import rdhte
    d = _load("parity_fixture.csv")
    # `bwjoint=True` forces a single shared bandwidth across cells, which
    # is what the R reference fixture was generated against. (Default is
    # now per-cell, matching R rdhte's `bw.joint = FALSE` default.)
    m = rdhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
              c=0.0, covs_hte=d["W"].to_numpy(),
              bwselect="msetwo", bwjoint=True)
    _check(m, R_REF["auto_bw_joint"], has_h=True)


def test_parity_cluster_manual_h():
    from rdhte import rdhte
    d = _load("parity_fixture_cluster.csv")
    m = rdhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
              c=0.0, covs_hte=d["W"].to_numpy(),
              cluster=d["cl"].to_numpy(), h=0.3)
    _check(m, R_REF["cluster"])


def test_no_hte_average_effect_runs():
    """ATE-only path (no covs_hte) should run cleanly."""
    from rdhte import rdhte
    d = _load("parity_fixture.csv")
    m = rdhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(), c=0.0, h=0.3)
    assert m.coef is not None
    assert "T" in m.coef.index
    assert m.se_rb[0] > 0

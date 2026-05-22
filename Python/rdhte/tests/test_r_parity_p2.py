"""Phase-2 R parity tests.

Covers the continuous-W path and rdhte_lincom. Frozen R reference values
are embedded below; the script that generated them is at the bottom.
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
        pytest.skip(f"missing fixture: {filename}")
    return pd.read_csv(path)


# Frozen R reference values.
R_REF_CONTINUOUS = {
    "T":          {"Estimate": 0.5485485, "Estimate.bc": 0.7601077, "se.rb": 0.1517786},
    "T:covs.hte": {"Estimate": 1.5047208, "Estimate.bc": 1.1774729, "se.rb": 0.3320932},
}

R_REF_LINCOM_CAT = {
    # Hypothesis: factor(W)A - factor(W)B = 0 (same as Python "A - B = 0")
    "A_minus_B": {"estimate": -0.392, "z": -2.556, "p": 0.011,
                  "ci": (-0.821, -0.108)},
    "A_eq_0":    {"estimate":  1.082, "z":  8.432, "p": 0.000,
                  "ci": ( 0.842,  1.352)},
    "joint":     {"chi2": 222.147, "df": 2, "p": 0.0},
}

R_REF_LINCOM_CLUSTER = {
    "A_minus_B": {"estimate": -0.838, "z": -4.354, "p": 0.0,
                  "ci": (-1.308, -0.496)},
}

R_REF_LINCOM_CONT_SLOPE = {
    "slope_eq_0": {"estimate": 1.505, "z": 3.546, "p": 0.0,
                   "ci": (0.527, 1.828)},
}

TOL_ROUND_3 = 1.5e-3   # R output is rounded to 3 digits
TOL_NUMERIC = 1e-5     # raw numeric outputs


# ----------------------------- Continuous W ---------------------------------

def test_parity_continuous_W_manual_h():
    from rdhte import rdhte
    d = _load("parity_fixture_continuous.csv")
    m = rdhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
              c=0.0, covs_hte=d["W"].to_numpy(), h=0.3)
    assert m.covs_cont is True
    for name, ref in R_REF_CONTINUOUS.items():
        assert abs(m.coef.loc[name]    - ref["Estimate"])    < TOL_NUMERIC, name
        assert abs(m.coef_bc.loc[name] - ref["Estimate.bc"]) < TOL_NUMERIC, name
        idx = list(m.coef.index).index(name)
        assert abs(m.se_rb[idx]        - ref["se.rb"])       < TOL_NUMERIC, name


# ----------------------------- lincom: categorical --------------------------

def test_parity_lincom_categorical_manual_h():
    from rdhte import rdhte, rdhte_lincom
    d = _load("parity_fixture.csv")
    m = rdhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
              c=0.0, covs_hte=d["W"].to_numpy(), h=0.3)
    res = rdhte_lincom(m, linfct=["A - B = 0", "A = 0"])
    ind = res["individual"]
    row_AB = ind.iloc[0]
    row_A  = ind.iloc[1]
    ref_AB = R_REF_LINCOM_CAT["A_minus_B"]
    ref_A  = R_REF_LINCOM_CAT["A_eq_0"]
    assert abs(row_AB["estimate"]  - ref_AB["estimate"])  < TOL_ROUND_3
    assert abs(row_AB["z_stat"]    - ref_AB["z"])         < TOL_ROUND_3
    assert abs(row_AB["p_value"]   - ref_AB["p"])         < TOL_ROUND_3
    assert abs(row_AB["conf.low"]  - ref_AB["ci"][0])     < TOL_ROUND_3
    assert abs(row_AB["conf.high"] - ref_AB["ci"][1])     < TOL_ROUND_3
    assert abs(row_A["estimate"]   - ref_A["estimate"])   < TOL_ROUND_3
    assert abs(row_A["z_stat"]     - ref_A["z"])          < TOL_ROUND_3
    joint = res["joint"].iloc[0]
    assert abs(joint["statistic"]  - R_REF_LINCOM_CAT["joint"]["chi2"]) < TOL_ROUND_3
    assert joint["df"] == R_REF_LINCOM_CAT["joint"]["df"]
    assert abs(joint["p_value"]    - R_REF_LINCOM_CAT["joint"]["p"])    < TOL_ROUND_3


# ----------------------------- lincom: cluster ------------------------------

def test_parity_lincom_cluster_manual_h():
    from rdhte import rdhte, rdhte_lincom
    d = _load("parity_fixture_cluster.csv")
    m = rdhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
              c=0.0, covs_hte=d["W"].to_numpy(),
              cluster=d["cl"].to_numpy(), h=0.3)
    res = rdhte_lincom(m, linfct="A - B = 0")
    row = res["individual"].iloc[0]
    ref = R_REF_LINCOM_CLUSTER["A_minus_B"]
    assert abs(row["estimate"]  - ref["estimate"])  < TOL_ROUND_3
    assert abs(row["z_stat"]    - ref["z"])         < TOL_ROUND_3
    assert abs(row["conf.low"]  - ref["ci"][0])     < TOL_ROUND_3
    assert abs(row["conf.high"] - ref["ci"][1])     < TOL_ROUND_3


# ----------------------------- lincom: continuous slope ---------------------

def test_parity_lincom_continuous_slope():
    from rdhte import rdhte, rdhte_lincom
    d = _load("parity_fixture_continuous.csv")
    m = rdhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
              c=0.0, covs_hte=d["W"].to_numpy(), h=0.3)
    res = rdhte_lincom(m, linfct="`T:covs.hte` = 0")
    row = res["individual"].iloc[0]
    ref = R_REF_LINCOM_CONT_SLOPE["slope_eq_0"]
    assert abs(row["estimate"]  - ref["estimate"])  < TOL_ROUND_3
    assert abs(row["z_stat"]    - ref["z"])         < TOL_ROUND_3
    assert abs(row["conf.low"]  - ref["ci"][0])     < TOL_ROUND_3
    assert abs(row["conf.high"] - ref["ci"][1])     < TOL_ROUND_3


def test_lincom_rejects_fraction_level():
    from rdhte import rdhte, rdhte_lincom
    d = _load("parity_fixture.csv")
    m = rdhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
              c=0.0, covs_hte=d["W"].to_numpy(), h=0.3)
    with pytest.raises(ValueError, match="percentage form"):
        rdhte_lincom(m, linfct="A - B = 0", level=0.95)


# R script that produced the reference values:
#
#   library(rdhte)
#   d  <- read.csv("parity_fixture.csv")
#   d2 <- read.csv("parity_fixture_cluster.csv")
#   d3 <- read.csv("parity_fixture_continuous.csv")
#   m  <- rdhte(y=d$y, x=d$x, c=0, covs.hte=factor(d$W), h=0.3)
#   r  <- rdhte_lincom(m, linfct=c("`factor(d$W)A` - `factor(d$W)B` = 0",
#                                  "`factor(d$W)A` = 0"))
#   mc <- rdhte(y=d2$y, x=d2$x, c=0, covs.hte=factor(d2$W), cluster=d2$cl, h=0.3)
#   rc <- rdhte_lincom(mc, linfct=c("`factor(d2$W)A` - `factor(d2$W)B` = 0"))
#   mw <- rdhte(y=d3$y, x=d3$x, c=0, covs.hte=d3$W, h=0.3)
#   rw <- rdhte_lincom(mw, linfct=c("`T:d3$W` = 0"))

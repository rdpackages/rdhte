"""Phase-0 smoke test.

Verify the package imports cleanly, exposes the documented public API,
and reports its version. Functional tests will land in Phase 1+.
"""

from __future__ import annotations

import pytest


def test_package_imports():
    import rdhte
    assert rdhte.__version__ == "0.1.0"


def test_public_api_is_exposed():
    import rdhte

    expected = {"rdhte", "rdbwhte", "rdhte_lincom", "RdhteResult", "RdbwhteResult"}
    missing = expected - set(rdhte.__all__)
    assert not missing, f"Missing from __all__: {missing}"
    for name in expected:
        assert hasattr(rdhte, name), f"rdhte module does not expose {name!r}"


def test_phase1_rdhte_runs_on_minimal_data():
    """Phase-1 contract: rdhte() runs and returns an RdhteResult.

    Phase-1 implementation covers categorical W with manual or auto bandwidth.
    rdbwhte() and rdhte_lincom() are still scaffolds and will raise.
    """
    import numpy as np
    from rdhte import rdhte, rdbwhte, rdhte_lincom, RdhteResult

    rng = np.random.default_rng(0)
    n = 500
    x = rng.uniform(-1.0, 1.0, n)
    W = np.where(rng.random(n) < 0.5, "A", "B")
    y = ((x >= 0).astype(float) * (1.0 + 0.5 * (W == "B"))
         + 0.3 * x + rng.normal(0, 0.3, n))

    m = rdhte(y, x, c=0.0, covs_hte=W, h=0.3)
    assert isinstance(m, RdhteResult)
    assert set(m.W_lev) == {"A", "B"}
    assert m.se_rb.shape == (2,)
    assert (m.se_rb > 0).all()

    # All four public functions are implemented as of Phase 4.
    bw = rdbwhte(y, x, c=0.0, covs_hte=W)
    assert bw.h.shape[0] == 2          # one row per W level

    res = rdhte_lincom(m, linfct="A - B = 0")
    assert "individual" in res
    assert "joint" in res
    assert res["joint"]["df"].iloc[0] == 1


def test_plot_module_imports():
    """Phase-4: plot.plot is implemented; just verify the import path works."""
    from rdhte.plot import plot
    # Calling plot() on a real fitted model is exercised in test_p4_plot.py;
    # here we just verify the function is importable and callable.
    assert callable(plot)


def test_dataclasses_round_trip():
    """Dataclass return objects are constructible with no args (all optional)."""
    from rdhte import RdhteResult, RdbwhteResult

    r = RdhteResult(Estimate=None)
    assert r.N == 0
    assert r.covs_cont is False
    assert isinstance(r.extras, dict)

    b = RdbwhteResult()
    assert b.N == 0
    assert isinstance(b.extras, dict)

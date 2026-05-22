"""Phase-4 tests for plot.rdhte (plotnine forest plot)."""

from __future__ import annotations

import os
import numpy as np
import pandas as pd
import pytest


HERE = os.path.dirname(__file__)


pytest.importorskip("plotnine")


def _load(filename):
    path = os.path.join(HERE, filename)
    if not os.path.exists(path):
        pytest.skip(f"missing fixture: {filename}")
    return pd.read_csv(path)


def test_plot_returns_ggplot_categorical():
    from rdhte import rdhte
    from rdhte.plot import plot
    import plotnine

    d = _load("parity_fixture.csv")
    m = rdhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
              c=0.0, covs_hte=d["W"].to_numpy(), h=0.3)
    p = plot(m)
    assert isinstance(p, plotnine.ggplot)
    # Three layers: hline (zero), errorbar, point.
    assert len(p.layers) == 3
    # Data has one row per group.
    assert len(p.data) == 2
    assert set(p.data["group"].astype(str)) == {"A", "B"}


def test_plot_sort_reorders():
    from rdhte import rdhte
    from rdhte.plot import plot

    d = _load("parity_fixture.csv")
    m = rdhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
              c=0.0, covs_hte=d["W"].to_numpy(), h=0.3)
    p = plot(m, sort=True)
    # After sort, estimates should be in ascending order.
    est = p.data["estimate"].to_numpy()
    assert np.all(np.diff(est) >= 0)


def test_plot_no_zero_line():
    from rdhte import rdhte
    from rdhte.plot import plot

    d = _load("parity_fixture.csv")
    m = rdhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
              c=0.0, covs_hte=d["W"].to_numpy(), h=0.3)
    p = plot(m, zero_line=False)
    # Without the hline layer: only errorbar + point.
    assert len(p.layers) == 2


def test_plot_continuous_W_errors():
    from rdhte import rdhte
    from rdhte.plot import plot

    d = _load("parity_fixture_continuous.csv")
    m = rdhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
              c=0.0, covs_hte=d["W"].to_numpy(), h=0.3)
    with pytest.raises(NotImplementedError, match="categorical"):
        plot(m)


def test_plot_data_matches_model():
    from rdhte import rdhte
    from rdhte.plot import plot

    d = _load("parity_fixture.csv")
    m = rdhte(y=d["y"].to_numpy(), x=d["x"].to_numpy(),
              c=0.0, covs_hte=d["W"].to_numpy(), h=0.3)
    p = plot(m, sort=False)
    # Estimate column should equal m.coef.
    pd.testing.assert_series_equal(
        p.data["estimate"].reset_index(drop=True),
        m.coef.reset_index(drop=True),
        check_names=False,
    )
    # CI columns equal m.ci_rb.
    assert np.allclose(p.data["ci_low"].to_numpy(),  m.ci_rb[:, 0])
    assert np.allclose(p.data["ci_high"].to_numpy(), m.ci_rb[:, 1])

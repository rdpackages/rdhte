"""Generate `parity_fixture.csv` and the R reference fit. Run once.

After this runs, `parity_fixture.csv` is the shared dataset both R and
Python read. `parity_fixture_R.csv` holds the R fit results.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(__file__)
OUT  = os.path.join(HERE, "parity_fixture.csv")


def main():
    rng = np.random.default_rng(seed=2026)
    n = 1500
    x = rng.uniform(-1.0, 1.0, n)
    W = np.where(rng.random(n) < 0.4, "A", "B")
    y = ((x >= 0).astype(float)
         * (1.0 * (W == "A").astype(float) + 1.5 * (W == "B").astype(float))
         + 0.3 * x
         + rng.normal(0, 0.3, n))
    pd.DataFrame({"y": y, "x": x, "W": W}).to_csv(OUT, index=False)
    print(f"Wrote {OUT}  (n = {n})")


if __name__ == "__main__":
    main()

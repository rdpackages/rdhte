"""Generate a fixture with a cluster variable for cluster-vcov parity test."""

from __future__ import annotations
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(__file__)
OUT  = os.path.join(HERE, "parity_fixture_cluster.csv")


def main():
    rng = np.random.default_rng(seed=4242)
    n = 2000
    x = rng.uniform(-1.0, 1.0, n)
    W = np.where(rng.random(n) < 0.4, "A", "B")
    cl = rng.integers(1, 51, size=n)               # 50 clusters
    eps_cl = rng.normal(0, 0.4, size=51)            # within-cluster shocks
    y = ((x >= 0).astype(float)
         * (1.0 * (W == "A").astype(float) + 1.5 * (W == "B").astype(float))
         + 0.3 * x
         + rng.normal(0, 0.3, n)
         + eps_cl[cl])                              # cluster shocks
    pd.DataFrame({"y": y, "x": x, "W": W, "cl": cl}).to_csv(OUT, index=False)
    print(f"Wrote {OUT}  (n = {n}, 50 clusters)")


if __name__ == "__main__":
    main()

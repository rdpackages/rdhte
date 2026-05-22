#!/usr/bin/env python
from __future__ import annotations

import json
import math
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


SCRIPT = Path(__file__).resolve()
REPO_ROOT = SCRIPT.parents[2]
SRC = REPO_ROOT / "Python" / "rdhte" / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from rdhte import rdbwhte, rdhte, rdhte_lincom  # noqa: E402


def clean_name(value: Any) -> str:
    name = "unnamed" if value is None or str(value) == "" else str(value)
    name = re.sub(r"[^0-9A-Za-z_]+", "_", name).strip("_")
    if not name:
        name = "unnamed"
    if name[0].isdigit():
        name = f"X{name}"
    return name


def num(value: Any) -> float | None:
    if value is None:
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def named_numbers(values: Iterable[Any], names: Iterable[Any] | None = None) -> dict[str, float | None]:
    vals = list(values)
    nms = list(names) if names is not None else [f"v{i + 1}" for i in range(len(vals))]
    return {clean_name(name): num(value) for name, value in zip(nms, vals)}


def matrix_numbers(
    values: Any,
    row_names: Iterable[Any] | None = None,
    col_names: Iterable[Any] | None = None,
) -> dict[str, float | None]:
    mat = np.asarray(values, dtype=float)
    if mat.ndim == 1:
        mat = mat.reshape(-1, 1)
    rows = list(row_names) if row_names is not None else [f"r{i + 1}" for i in range(mat.shape[0])]
    cols = list(col_names) if col_names is not None else (
        ["left", "right"] if mat.shape[1] == 2 else [f"c{j + 1}" for j in range(mat.shape[1])]
    )
    out: dict[str, float | None] = {}
    for i, row in enumerate(rows):
        for j, col in enumerate(cols):
            out[f"{clean_name(row)}_{clean_name(col)}"] = num(mat[i, j])
    return out


def coef_names(obj: Any, override: Iterable[str] | None = None) -> list[str]:
    if override is not None:
        return list(override)
    if getattr(obj, "coef", None) is not None:
        return [str(x) for x in obj.coef.index]
    return [str(x) for x in (obj.W_lev or [])]


def h_rows(obj: Any) -> list[str] | None:
    h = np.asarray(obj.h)
    if obj.W_lev is not None and h.ndim == 2 and len(obj.W_lev) == h.shape[0]:
        return [str(x) for x in obj.W_lev]
    return None


def summarize_rdhte(
    obj: Any,
    names: Iterable[str] | None = None,
    n_lr: Iterable[int] | None = None,
) -> dict[str, Any]:
    nms = coef_names(obj, names)
    ci = np.asarray(obj.ci_rb)
    return {
        "estimate": named_numbers(obj.coef.to_numpy(), nms),
        "estimate_bc": named_numbers(obj.coef_bc.to_numpy(), nms),
        "se_rb": named_numbers(obj.se_rb, nms),
        "t_rb": named_numbers(obj.t_rb, nms),
        "pv_rb": named_numbers(obj.pv_rb, nms),
        "ci_rb": matrix_numbers(ci, nms, ["lower", "upper"]),
        "h": matrix_numbers(obj.h, h_rows(obj), ["left", "right"]),
        "n": named_numbers(n_lr if n_lr is not None else [obj.N], None if n_lr is not None else ["N"]),
        "nh": matrix_numbers(obj.Nh, h_rows(obj), ["left", "right"]),
        "vcov": matrix_numbers(obj.vcov, nms, nms),
        "kernel": obj.kernel,
        "vce": obj.vce,
        "bwselect": obj.bwselect,
        "rdmodel": obj.rdmodel,
    }


def summarize_rdbwhte(obj: Any, n_lr: Iterable[int] | None = None) -> dict[str, Any]:
    return {
        "h": matrix_numbers(obj.h, h_rows(obj), ["left", "right"]),
        "n": named_numbers(n_lr if n_lr is not None else [obj.N], None if n_lr is not None else ["N"]),
        "nh": matrix_numbers(obj.Nh, h_rows(obj), ["left", "right"]),
        "kernel": obj.kernel,
        "vce": obj.vce,
        "bwselect": obj.bwselect,
        "rdmodel": obj.rdmodel,
    }


def summarize_lincom(obj: dict[str, pd.DataFrame]) -> dict[str, Any]:
    individual = obj["individual"]
    joint = obj["joint"]
    stat_col = "statistic_chi2" if "statistic_chi2" in joint.columns else "statistic"
    return {
        "individual": [
            {
                "hypothesis": row["hypothesis"],
                "estimate": num(row["estimate"]),
                "z_stat": num(row["z_stat"]),
                "p_value": num(row["p_value"]),
                "conf_low": num(row["conf.low"]),
                "conf_high": num(row["conf.high"]),
            }
            for _, row in individual.iterrows()
        ],
        "joint": {
            "statistic_chi2": num(joint.loc[0, stat_col]),
            "df": num(joint.loc[0, "df"]),
            "p_value": num(joint.loc[0, "p_value"]),
        },
    }


def main() -> int:
    output = Path(sys.argv[1]) if len(sys.argv) >= 2 else REPO_ROOT / "docs" / "audit" / "baselines" / "python-current.json"
    output.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(REPO_ROOT / "R" / "rdhte_dataset.csv")
    n_lr = [int((df["x"] < 0).sum()), int((df["x"] >= 0).sum())]

    rd_left = rdhte(y="y", x="x", covs_hte=pd.Categorical(df["w_left"]), cluster="cluster_var", data=df)
    rd_left_joint = rdhte(y="y", x="x", covs_hte=pd.Categorical(df["w_left"]), cluster="cluster_var", bwjoint=True, data=df)
    bw_left = rdbwhte(y="y", x="x", covs_hte=pd.Categorical(df["w_left"]), cluster="cluster_var", data=df)
    rd_ideology = rdhte(y="y", x="x", covs_hte=pd.Categorical(df["w_ideology"]), cluster="cluster_var", data=df)
    rd_strength = rdhte(y="y", x="x", covs_hte=df["w_strength"].to_numpy(), kernel="uni", cluster="cluster_var", data=df)
    rd_interaction = rdhte(y="y", x="x", covs_hte="w_left*w_strength", h=0.1, cluster="cluster_var", data=df)
    rd_average = rdhte(y="y", x="x", h=0.1, vce="hc3", data=df)

    cases = {
        "binary_left": {
            "rdhte": summarize_rdhte(rd_left, ["factor(w_left)0", "factor(w_left)1"], n_lr),
            "lincom": summarize_lincom(rdhte_lincom(rd_left, linfct=["`1` - `0` = 0"], digits=12)),
        },
        "binary_left_joint": {
            "rdhte": summarize_rdhte(rd_left_joint, ["w_left0", "w_left1"], n_lr),
        },
        "binary_left_bw": {
            "rdbwhte": summarize_rdbwhte(bw_left, n_lr),
        },
        "categorical_ideology": {
            "rdhte": summarize_rdhte(
                rd_ideology,
                ["factor(w_ideology)1", "factor(w_ideology)2", "factor(w_ideology)3", "factor(w_ideology)4"],
                n_lr,
            ),
        },
        "continuous_strength": {
            "rdhte": summarize_rdhte(rd_strength, ["T", "T:w_strength"], n_lr),
        },
        "interaction_strength": {
            "rdhte": summarize_rdhte(
                rd_interaction,
                ["T", "T:w_left", "T:w_strength", "T:w_left:w_strength"],
                n_lr,
            ),
        },
        "average_manual": {
            "rdhte": summarize_rdhte(rd_average, ["T"], n_lr),
        },
    }

    result = {
        "schema_version": 1,
        "package": "rdhte",
        "language": "python",
        "source": "working-tree",
        "timestamp_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "environment": {
            "python": sys.version.split()[0],
            "platform": sys.platform,
        },
        "cases": cases,
    }

    output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any


RDHTE_CASES = [
    "binary_left",
    "binary_left_joint",
    "categorical_ideology",
    "continuous_strength",
    "interaction_strength",
    "average_manual",
]


def values(obj: dict[str, Any]) -> list[Any]:
    return list(obj.values())


def matrix_values(obj: Any) -> list[Any]:
    out: list[Any] = []
    for row in obj:
        if isinstance(row, list):
            out.extend(row)
        else:
            out.append(row)
    return out


def lower_values(obj: dict[str, Any]) -> list[Any]:
    return [value for key, value in obj.items() if key.endswith("_lower")]


def upper_values(obj: dict[str, Any]) -> list[Any]:
    return [value for key, value in obj.items() if key.endswith("_upper")]


def is_num(value: Any) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def compare_vec(
    label: str,
    expected: list[Any],
    actual: list[Any],
    failures: list[str],
    *,
    atol: float,
    rtol: float,
    prefix_ok: bool = False,
) -> None:
    exp = [x for x in expected if x is not None]
    act = [x for x in actual if x is not None]
    if prefix_ok:
        act = act[: len(exp)]
    if len(exp) != len(act):
        failures.append(f"{label}: length mismatch expected {len(exp)}, actual {len(act)}")
        return
    for i, (a, b) in enumerate(zip(exp, act)):
        if not (is_num(a) and is_num(b)):
            continue
        if not math.isclose(float(a), float(b), rel_tol=rtol, abs_tol=atol):
            failures.append(
                f"{label}.{i}: expected {float(a):.17g}, actual {float(b):.17g}, "
                f"diff {float(b) - float(a):.3g}"
            )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare an R/Python named baseline against the Stata baseline."
    )
    parser.add_argument("reference")
    parser.add_argument("stata")
    parser.add_argument("--atol", type=float, default=1e-5)
    parser.add_argument("--rtol", type=float, default=1e-5)
    args = parser.parse_args()

    reference = json.loads(Path(args.reference).read_text(encoding="utf-8"))
    stata = json.loads(Path(args.stata).read_text(encoding="utf-8"))

    failures: list[str] = []
    compared = 0

    for case in RDHTE_CASES:
        ref = reference["cases"][case]["rdhte"]
        st = stata["cases"][case]["rdhte"]
        pairs = [
            ("estimate", values(ref["estimate"]), matrix_values(st["tau_hat"]), False),
            ("estimate_bc", values(ref["estimate_bc"]), matrix_values(st["tau_bc"]), False),
            ("se_rb", values(ref["se_rb"]), matrix_values(st["tau_se"]), False),
            ("t_rb", values(ref["t_rb"]), matrix_values(st["tau_t"]), False),
            ("pv_rb", values(ref["pv_rb"]), matrix_values(st["tau_pv"]), False),
            ("ci_low", lower_values(ref["ci_rb"]), matrix_values(st["tau_ci_lb"]), False),
            ("ci_high", upper_values(ref["ci_rb"]), matrix_values(st["tau_ci_ub"]), False),
            ("h", values(ref["h"]), matrix_values(st["h"]), True),
            ("nh", values(ref["nh"]), matrix_values(st["tau_n"]), True),
        ]
        if (
            ref.get("vcov") is not None
            and "tau_v" in st
            and any(x is not None for x in matrix_values(st["tau_v"]))
        ):
            pairs.append(("vcov", values(ref["vcov"]), matrix_values(st["tau_v"]), False))
        for metric, left, right, prefix_ok in pairs:
            compared += min(len(left), len(right))
            compare_vec(
                f"{case}.{metric}",
                left,
                right,
                failures,
                atol=args.atol,
                rtol=args.rtol,
                prefix_ok=prefix_ok,
            )

    ref_bw = reference["cases"]["binary_left_bw"]["rdbwhte"]
    st_bw = stata["cases"]["binary_left_bw"]["rdbwhte"]
    for metric, left, right in [
        ("binary_left_bw.h", values(ref_bw["h"]), matrix_values(st_bw["h"])),
        ("binary_left_bw.nh", values(ref_bw["nh"]), matrix_values(st_bw["tau_n"])),
    ]:
        compared += min(len(left), len(right))
        compare_vec(metric, left, right, failures, atol=args.atol, rtol=args.rtol)

    ref_lc = reference["cases"]["binary_left"]["lincom"]["individual"][0]
    st_lc = stata["cases"]["binary_left_lincom"]["lincom"]
    for metric in ["estimate", "conf_low", "conf_high", "p_value"]:
        compared += 1
        compare_vec(
            f"binary_left_lincom.{metric}",
            [ref_lc[metric]],
            [st_lc[metric]],
            failures,
            atol=args.atol,
            rtol=args.rtol,
        )

    if failures:
        for failure in failures:
            print(failure)
        print(f"Compared {compared} numeric entries; {len(failures)} differences outside tolerance.")
        return 1

    print(f"Compared {compared} numeric entries; all matched within tolerance.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

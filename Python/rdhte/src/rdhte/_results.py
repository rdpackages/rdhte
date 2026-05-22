"""Return-value dataclasses for rdhte() and rdbwhte().

The dataclasses keep estimates, inference results, bandwidths, labels, and
metadata together while allowing downstream code to extend the object.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional

import numpy as np
import pandas as pd


def _fmt_int(v, w: int = 10) -> str:
    return f"{int(v):>{w}d}"


def _fmt_str(v, w: int = 10, just: str = "right") -> str:
    s = str(v)
    return f"{s:>{w}s}" if just == "right" else f"{s:<{w}s}"


def _fmt_f(v, w: int = 10, k: int = 3) -> str:
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return f"{'NA':>{w}s}"
    return f"{v:>{w}.{k}f}"


@dataclass
class RdhteResult:
    """Return value of :func:`rdhte.rdhte`.

    Field names expose the main estimates, robust bias-corrected inference,
    bandwidths, labels, and fit metadata.
    """

    Estimate: pd.DataFrame                 # one row per group (categorical W) or one row per coef (continuous W)
    Estimate_bc: Optional[pd.DataFrame] = None
    se_rb: Optional[np.ndarray] = None
    coef: Optional[pd.Series] = None        # per-group CATEs (categorical) or [T, T:covs.hte] (continuous)
    coef_bc: Optional[pd.Series] = None
    vcov: Optional[np.ndarray] = None
    ci_rb: Optional[np.ndarray] = None      # shape (n_groups, 2)
    t_rb: Optional[np.ndarray] = None
    pv_rb: Optional[np.ndarray] = None
    coef_full: Optional[pd.Series] = None
    vcov_full: Optional[np.ndarray] = None
    W_lev: Optional[list] = None
    W_names: Optional[list] = None
    covs_hte_chr: Optional[str] = None
    kernel: Optional[str] = None
    bwselect: Optional[str] = None
    vce: Optional[str] = None
    vce_select: Optional[str] = None        # canonical lowercase form ("cr1", "hc3", ...)
    c: float = 0.0
    h: Optional[np.ndarray] = None
    p: int = 1
    q: int = 2
    N: int = 0
    Nh: Optional[np.ndarray] = None
    covs_cont: bool = False                 # True for continuous W
    level: float = 95.0
    rdmodel: Optional[str] = None
    call: Optional[Any] = None

    extras: dict = field(default_factory=dict)

    def __repr__(self) -> str:
        """Short header."""
        h0 = self.h[0] if self.h is not None and self.h.size else (np.nan, np.nan)
        lines = [
            "RD Heterogeneous Treatment Effects Estimation",
            f"Number of Obs.           {self.N:>10d}",
            f"Kernel                   {self.kernel!s:>10s}",
            f"VCE method               {self.vce!s:>10s}",
            f"Poly. Order (p)          {int(self.p):>10d}",
            f"Bandwidth (h)            {h0[0]:>10.3f} {h0[1]:>10.3f}",
        ]
        return "\n".join(lines)

    def summary(self) -> str:
        """Full table."""
        Nh = self.Nh if self.Nh is not None else np.zeros((1, 2), dtype=int)
        h  = self.h  if self.h  is not None else np.zeros((1, 2))
        nrows = len(self.W_lev or [])

        # Header block.
        out = []
        out.append(str(self.rdmodel or "rdhte"))
        out.append("")
        out.append(f"Number of Obs.           {self.N:>10d}")
        out.append(f"BW type                  {self.bwselect!s:>10s}")
        out.append(f"Kernel                   {self.kernel!s:>10s}")
        out.append(f"VCE method               {self.vce!s:>10s}")
        out.append("")
        # Per-side N and Nh (use the shared bandwidth row for the categorical case).
        Nh0 = Nh[0] if Nh.shape[0] else (0, 0)
        out.append(f"Eff. Number of Obs.      {int(Nh0[0]):>10d}    {int(Nh0[1]):>10d}")
        out.append(f"Order est. (p)           {self.p:>10d}    {self.p:>10d}")
        out.append(f"Order bias  (q)          {self.q:>10d}    {self.q:>10d}")
        if self.covs_cont:
            h0 = h[0] if h.shape[0] else (np.nan, np.nan)
            out.append(f"BW est. (h)              {h0[0]:>10.3f}    {h0[1]:>10.3f}")
        out.append("")

        # Estimate table.
        name_col = (self.covs_hte_chr or "Group")
        name_w   = max(12, len(name_col))
        cont = self.covs_cont
        if cont:
            tot_w = name_w + 12 + 12 + 12 + 24
        else:
            tot_w = name_w + 12 + 12 + 12 + 24 + 10 + 10 + 10 + 10
        out.append("=" * tot_w)
        # Top header row.
        hdr = (
            _fmt_str("", name_w)
            + _fmt_str("Point", 12)
            + _fmt_str("Robust Inference", 12 + 12 + 24)
        )
        out.append(hdr)
        # Sub-header.
        sub  = _fmt_str(name_col if not cont else "", name_w)
        sub += _fmt_str("Estimate", 12)
        sub += _fmt_str("z", 12)
        sub += _fmt_str("Pr(>|z|)", 12)
        sub += f"[ {int(self.level)}% C.I. ]".center(24)
        if not cont:
            sub += _fmt_str("Nh-", 10) + _fmt_str("Nh+", 10)
            sub += _fmt_str("h-", 10)  + _fmt_str("h+", 10)
        out.append(sub)
        out.append("-" * tot_w)

        est    = np.asarray(self.Estimate["Estimate"], dtype=float)
        t_rb   = np.asarray(self.t_rb,  dtype=float) if self.t_rb  is not None else np.zeros(nrows)
        pv_rb  = np.asarray(self.pv_rb, dtype=float) if self.pv_rb is not None else np.zeros(nrows)
        ci     = np.asarray(self.ci_rb, dtype=float) if self.ci_rb is not None else np.zeros((nrows, 2))
        labels = (self.W_names if cont else self.W_lev) or []

        # In the categorical case, all groups share the same (h, Nh)
        # (set in rdhte()'s build); broadcast h[0] / Nh[0] for every row.
        for i in range(nrows):
            row  = _fmt_str(labels[i], name_w)
            row += _fmt_f(est[i],   12, 3)
            row += _fmt_f(t_rb[i],  12, 3)
            row += _fmt_f(pv_rb[i], 12, 3)
            left  = f"[{ci[i, 0]:>.3f} , "
            right = f"{ci[i, 1]:>.3f}]"
            row += f"{left:>12s}" + f"{right:<12s}"
            if not cont:
                hi = h[i]  if h.shape[0]  > i else h[0]
                ni = Nh[i] if Nh.shape[0] > i else Nh[0]
                row += _fmt_int(ni[0], 10) + _fmt_int(ni[1], 10)
                row += _fmt_f(hi[0], 10, 3) + _fmt_f(hi[1], 10, 3)
            out.append(row)
        out.append("=" * tot_w)
        return "\n".join(out)

    def tidy(self) -> pd.DataFrame:
        """One row per heterogeneity term, broom-style.

        Columns: ``term``, ``estimate``, ``std.error``, ``statistic``,
        ``p.value``, ``conf.low``, ``conf.high``, ``estimate.bc``,
        ``h.left``, ``h.right``, ``n.eff.left``, ``n.eff.right``.
        """
        cont   = self.covs_cont
        terms  = (self.W_names if cont else self.W_lev) or []
        nrows  = len(terms)
        est    = np.asarray(self.Estimate["Estimate"], dtype=float)
        est_bc = np.asarray(self.Estimate["Estimate.bc"], dtype=float)
        se_rb  = np.asarray(self.se_rb, dtype=float) if self.se_rb is not None else np.full(nrows, np.nan)
        t_rb   = np.asarray(self.t_rb,  dtype=float) if self.t_rb  is not None else np.full(nrows, np.nan)
        pv_rb  = np.asarray(self.pv_rb, dtype=float) if self.pv_rb is not None else np.full(nrows, np.nan)
        ci     = np.asarray(self.ci_rb, dtype=float) if self.ci_rb is not None else np.full((nrows, 2), np.nan)

        h  = self.h  if self.h  is not None else np.zeros((1, 2))
        Nh = self.Nh if self.Nh is not None else np.zeros((1, 2), dtype=int)
        h_l  = np.array([h[i if h.shape[0]  > i else 0, 0] for i in range(nrows)], dtype=float)
        h_r  = np.array([h[i if h.shape[0]  > i else 0, 1] for i in range(nrows)], dtype=float)
        nh_l = np.array([Nh[i if Nh.shape[0] > i else 0, 0] for i in range(nrows)], dtype=float)
        nh_r = np.array([Nh[i if Nh.shape[0] > i else 0, 1] for i in range(nrows)], dtype=float)

        return pd.DataFrame({
            "term":         terms,
            "estimate":     est,
            "std.error":    se_rb,
            "statistic":    t_rb,
            "p.value":      pv_rb,
            "conf.low":     ci[:, 0],
            "conf.high":    ci[:, 1],
            "estimate.bc":  est_bc,
            "h.left":       h_l,
            "h.right":      h_r,
            "n.eff.left":   nh_l,
            "n.eff.right":  nh_r,
        })

    def glance(self) -> pd.DataFrame:
        """One-row summary of the fit."""
        return pd.DataFrame([{
            "n":               int(self.N),
            "cutoff":          float(self.c),
            "p":               int(self.p),
            "q":               int(self.q),
            "kernel":          self.kernel,
            "vce":             self.vce,
            "vce_select":      self.vce_select,
            "bwselect":        self.bwselect,
            "level":           float(self.level),
            "n.terms":         len(self.W_lev or []),
            "covs.continuous": bool(self.covs_cont),
        }])


@dataclass
class RdbwhteResult:
    """Return value of :func:`rdhte.rdbwhte`. Mirror of the R object."""

    h: Optional[np.ndarray] = None
    W_lev: Optional[list] = None
    W_names: Optional[list] = None
    covs_hte_chr: Optional[str] = None
    kernel: Optional[str] = None
    bwselect: Optional[str] = None
    vce: Optional[str] = None
    vce_select: Optional[str] = None        # canonical lowercase form
    c: float = 0.0
    p: int = 1
    q: int = 2
    N: int = 0
    Nh: Optional[np.ndarray] = None
    covs_cont: bool = False
    rdmodel: Optional[str] = None
    call: Optional[Any] = None

    extras: dict = field(default_factory=dict)

    def __repr__(self) -> str:
        lines = [
            "MSE-Optimal Bandwidth Selection for RD Heterogeneous Treatment Effects",
            f"Number of Obs.           {self.N:>10d}",
            f"Kernel                   {self.kernel!s:>10s}",
            f"VCE method               {self.vce!s:>10s}",
            f"Poly. Order (p)          {int(self.p):>10d}",
        ]
        return "\n".join(lines)

    def summary(self) -> str:
        out = []
        out.append(str(self.rdmodel or "rdbwhte"))
        out.append("")
        out.append(f"Number of Obs.           {self.N:>10d}")
        out.append(f"BW type                  {self.bwselect!s:>10s}")
        out.append(f"Kernel                   {self.kernel!s:>10s}")
        out.append(f"VCE method               {self.vce!s:>10s}")
        out.append("")
        out.append(f"Order est. (p)           {self.p:>10d}")
        out.append(f"Order bias  (q)          {self.q:>10d}")
        out.append("")

        tot_w = 15 + 15 + 15
        out.append("=" * tot_w)
        out.append(_fmt_str("Group" if not self.covs_cont else "", 15)
                   + _fmt_str("h-", 15) + _fmt_str("h+", 15))
        out.append("-" * tot_w)

        h = self.h if self.h is not None else np.zeros((1, 2))
        if not self.covs_cont:
            for i, lev in enumerate(self.W_lev or []):
                out.append(_fmt_str(lev, 15) + _fmt_f(h[i, 0], 15, 3) + _fmt_f(h[i, 1], 15, 3))
        else:
            row0 = h[0] if h.size else (np.nan, np.nan)
            out.append(_fmt_str("Overall", 15) + _fmt_f(row0[0], 15, 3) + _fmt_f(row0[1], 15, 3))
        out.append("=" * tot_w)
        return "\n".join(out)

    def tidy(self) -> pd.DataFrame:
        h = self.h if self.h is not None else np.zeros((1, 2))
        if self.covs_cont:
            return pd.DataFrame({"term": ["Overall"],
                                 "h.left": [float(h[0, 0])],
                                 "h.right": [float(h[0, 1])]})
        terms = list(self.W_lev or [])
        return pd.DataFrame({
            "term":    terms,
            "h.left":  [float(h[i, 0]) for i in range(len(terms))],
            "h.right": [float(h[i, 1]) for i in range(len(terms))],
        })

    def glance(self) -> pd.DataFrame:
        return pd.DataFrame([{
            "n":               int(self.N),
            "cutoff":          float(self.c),
            "p":               int(self.p),
            "q":               int(self.q),
            "kernel":          self.kernel,
            "vce":             self.vce,
            "vce_select":      self.vce_select,
            "bwselect":        self.bwselect,
            "n.terms":         len(self.W_lev or []),
            "covs.continuous": bool(self.covs_cont),
        }])

"""rdbwhte: standalone data-driven bandwidth selection for rdhte.

For each level of a categorical ``covs_hte``, calls
``rdrobust.rdbwselect`` on the cell subsample. For continuous
``covs_hte``, formula-string ``covs_hte``, or no ``covs_hte``, runs
``rdbwselect`` once on the full sample. Returns a
:class:`rdhte.RdbwhteResult` with per-cell ``h`` and ``Nh`` matrices.
``covs_eff`` is forwarded to ``rdrobust.rdbwselect`` via its ``covs=``
argument so the bandwidth reflects the efficiency covariates.
"""

from __future__ import annotations

import warnings
from typing import Optional, Union

import numpy as np
import pandas as pd

from ._helpers import (
    resolve_arg, resolve_covs_hte, combine_factors,
    resolve_formula_covs_hte, resolve_covs_eff,
)
from ._results import RdbwhteResult
from ._vce import _VALID_VCE, normalize_vce
from .rdhte import _arr, _looks_categorical  # reuse private helpers

ArrayLike = Union[np.ndarray, pd.Series, list]


def rdbwhte(
    y: ArrayLike,
    x: ArrayLike,
    c: float = 0.0,
    covs_hte: Optional[ArrayLike] = None,
    covs_eff: Optional[ArrayLike] = None,
    p: int = 1,
    q: Optional[int] = None,
    kernel: str = "tri",
    weights: Optional[ArrayLike] = None,
    bwselect: str = "mserd",
    vce: str = "hc3",
    cluster: Optional[ArrayLike] = None,
    subset: Optional[ArrayLike] = None,
    bwjoint: bool = False,
    data: Optional[pd.DataFrame] = None,
) -> RdbwhteResult:
    """Data-driven bandwidth selection for heterogeneous RD treatment effects.

    Parameters match :func:`rdhte.rdhte` minus the inference-only options
    (``vce``, ``level``, ``h_l``, ``h_r``, ``h``). Accepted ``covs_hte``
    forms:

    - ``None`` -> pooled single-bandwidth for the ATE.
    - Categorical -> one ``(h_l, h_r)`` row per cell (or one shared row
      when ``bwjoint=True``).
    - Continuous -> one pooled ``(h_l, h_r)`` row.
    - Multi-factor interaction (list / DataFrame / 2-D ndarray) ->
      Cartesian-product cells, one row per non-empty cell.
    - Formula string (requires ``data=``) -> one pooled ``(h_l, h_r)``
      row labelled ``"Overall"``.

    ``covs_eff`` is forwarded to the underlying ``rdrobust.rdbwselect``
    via its ``covs=`` argument so the chosen bandwidth reflects the
    efficiency covariates.

    The accepted bandwidth selectors are ``"mserd"``, ``"msetwo"``,
    ``"msesum"``, ``"msecomb1"``, ``"msecomb2"``, ``"cerrd"``,
    ``"certwo"``, ``"cersum"``, ``"cercomb1"``, and ``"cercomb2"``.
    MSE denotes mean square error; CER denotes coverage error rate.

    Returns
    -------
    RdbwhteResult
        Dataclass with ``h`` (``n_lev x 2`` ndarray of left/right
        bandwidths), ``Nh`` (effective sample sizes per cell), ``W_lev``,
        plus the usual metadata. Convenience methods: ``summary()``,
        ``tidy()``, ``glance()``.
    """
    from rdrobust import rdbwselect

    # Resolve `data=` lookups + coerce.
    y        = np.asarray(_arr(resolve_arg(y,        data, "y")),        dtype=float)
    x        = np.asarray(_arr(resolve_arg(x,        data, "x")),        dtype=float)
    cluster  = _arr(resolve_arg(cluster,  data, "cluster"))
    weights  = _arr(resolve_arg(weights,  data, "weights"))
    subset   = _arr(resolve_arg(subset,   data, "subset"))
    covs_eff = resolve_arg(covs_eff, data, "covs_eff")

    n_orig = len(y)
    covs_eff_mat = resolve_covs_eff(covs_eff, n_orig)

    # Formula-string covs_hte: materialize the patsy matrix on full-length data.
    # BARE-COLUMN-NAME shortcut (mirrors rdhte.py): if `data` is supplied and
    # `covs_hte` is a bare identifier matching a column of `data`, look it up
    # and dispatch through the array path so 0/1 binaries auto-promote to
    # factor (matching R + Stata).
    formula_W_mat = None
    formula_W_chr = None
    _bare_lookup = (
        isinstance(covs_hte, str)
        and data is not None
        and covs_hte.isidentifier()
        and covs_hte in data.columns
    )
    if isinstance(covs_hte, str) and not _bare_lookup:
        W_mat_df = resolve_formula_covs_hte(covs_hte, data)
        if len(W_mat_df) != n_orig:
            raise ValueError(
                "Formula-built covs_hte matrix length does not match `y` length."
            )
        formula_W_mat = W_mat_df.values.astype(float, copy=False)
        formula_W_chr = covs_hte
        covs_components = None
    else:
        covs_hte = resolve_arg(covs_hte, data, "covs_hte")
        covs_components = resolve_covs_hte(covs_hte, n_orig)

    if subset is not None:
        if subset.dtype == bool:
            mask = subset
        else:
            mask = np.zeros(n_orig, dtype=bool); mask[np.asarray(subset, dtype=int)] = True
        y = y[mask]; x = x[mask]
        if covs_components is not None:
            covs_components = [c[mask] for c in covs_components]
        if formula_W_mat is not None:
            formula_W_mat = formula_W_mat[mask]
        if covs_eff_mat is not None:
            covs_eff_mat = covs_eff_mat[mask]
        if cluster  is not None: cluster  = cluster[mask]
        if weights  is not None: weights  = weights[mask]

    # Drop NA on core columns.
    keep = np.isfinite(y) & np.isfinite(x)
    if covs_components is not None:
        for ch in covs_components:
            ch = np.asarray(ch)
            if np.issubdtype(ch.dtype, np.number):
                keep &= np.isfinite(ch.astype(float))
            else:
                keep &= pd.notna(pd.Series(ch)).to_numpy()
    if formula_W_mat is not None:
        keep &= np.isfinite(formula_W_mat).all(axis=1)
    if covs_eff_mat is not None:
        keep &= np.isfinite(covs_eff_mat).all(axis=1)
    if cluster is not None:
        cl = np.asarray(cluster)
        if np.issubdtype(cl.dtype, np.number):
            keep &= np.isfinite(cl.astype(float))
        else:
            keep &= pd.notna(pd.Series(cl)).to_numpy()
    if weights is not None:
        wgt = np.asarray(weights, dtype=float)
        keep &= np.isfinite(wgt) & (wgt > 0)
    y, x = y[keep], x[keep]
    if covs_components is not None:
        covs_components = [c[keep] for c in covs_components]
    if formula_W_mat is not None:
        formula_W_mat = formula_W_mat[keep]
    if covs_eff_mat is not None:
        covs_eff_mat = covs_eff_mat[keep]
    if cluster  is not None: cluster  = cluster[keep]
    if weights  is not None: weights  = weights[keep]

    # Resolve covs_hte to either None, a 1-D array, or a combined Categorical.
    if covs_components is None:
        covs_hte = None
    elif len(covs_components) == 1:
        covs_hte = np.asarray(covs_components[0])
    else:
        for j, comp in enumerate(covs_components):
            s = pd.Series(np.asarray(comp))
            if not _looks_categorical(s):
                raise ValueError(
                    f"Multi-component covs_hte requires all components to be "
                    f"categorical (component {j} looks continuous)."
                )
        covs_hte = combine_factors(covs_components)

    if q is None:
        q = p + 1
    kernel = kernel.lower()
    kernel = {"triangular": "tri", "epanechnikov": "epa", "uniform": "uni"}.get(kernel, kernel)
    Kernel = {"tri": "Triangular", "epa": "Epanechnikov", "uni": "Uniform"}[kernel]

    Xc = x - c
    T  = (Xc >= 0).astype(float)

    # vce normalization via shared helper (mirrors rdhte). The
    # canonical name maps to nn/hc0-hc3 below for the rdbwselect call,
    # which only accepts those forms.
    user_vce = vce is not None and vce != "hc3"
    vce, vce_label = normalize_vce(vce, cluster, user_vce)
    bw_vce_internal = {"cr1": "hc1", "cr2": "hc2", "cr3": "hc3"}.get(vce, vce)
    if bw_vce_internal not in {"nn", "hc0", "hc1", "hc2", "hc3"}:
        bw_vce_internal = "hc3"

    def _bw_once(sub_idx=None):
        """Call rdrobust.rdbwselect on a (possibly subset) sample."""
        if sub_idx is None:
            ys, xs = y, Xc
            cl = cluster
            ce = covs_eff_mat
        else:
            ys, xs = y[sub_idx], Xc[sub_idx]
            cl = cluster[sub_idx] if cluster is not None else None
            ce = covs_eff_mat[sub_idx, :] if covs_eff_mat is not None else None
        bw = rdbwselect(y=ys, x=xs, c=0.0, p=p, q=q, kernel=kernel,
                        bwselect=bwselect, vce=bw_vce_internal, cluster=cl,
                        covs=ce)
        bws_row = bw.bws.iloc[0]
        return float(bws_row.iloc[0]), float(bws_row.iloc[1])

    # Formula-string covs_hte: treat as pooled single-bandwidth (matches R rdbwhte).
    if formula_W_mat is not None:
        h_l, h_r = _bw_once()
        h_mat  = np.array([[h_l, h_r]])
        nh_mat = np.array([[int(((Xc < 0) & (np.abs(Xc) <= h_l)).sum()),
                            int(((Xc >= 0) & (np.abs(Xc) <= h_r)).sum())]])
        return RdbwhteResult(
            h=h_mat,
            W_lev=["Overall"], W_names=["Overall"],
            covs_hte_chr=formula_W_chr,
            kernel=Kernel, bwselect=bwselect, vce=vce_label, vce_select=vce,
            c=c, p=p, q=q, N=int(len(y)), Nh=nh_mat,
            covs_cont=True,
            rdmodel="MSE-Optimal BW for RD HTE (formula covs.hte).",
        )

    # Decide categorical vs continuous mode.
    if covs_hte is None:
        # ATE-only.
        h_l, h_r = _bw_once()
        h_mat  = np.array([[h_l, h_r]])
        nh_mat = np.array([[int(((Xc < 0) & (np.abs(Xc) <= h_l)).sum()),
                            int(((Xc >= 0) & (np.abs(Xc) <= h_r)).sum())]])
        return RdbwhteResult(
            h=h_mat,
            W_lev=["T"], W_names=["T"], covs_hte_chr=None,
            kernel=Kernel, bwselect=bwselect, vce=vce_label, vce_select=vce,
            c=c, p=p, q=q, N=int(len(y)), Nh=nh_mat,
            covs_cont=False,
            rdmodel="MSE-Optimal BW for RD Average Treatment Effect.",
        )

    if isinstance(covs_hte, pd.Categorical):
        W_cat_pre = covs_hte
        is_factor = True
    else:
        W_arr = pd.Series(np.asarray(covs_hte))
        is_factor = _looks_categorical(W_arr)
        W_cat_pre = pd.Categorical(W_arr) if is_factor else None
    if not is_factor:
        # Continuous covs_hte -> one shared bandwidth pair (R behavior).
        h_l, h_r = _bw_once()
        h_mat  = np.array([[h_l, h_r]])
        nh_mat = np.array([[int(((Xc < 0) & (np.abs(Xc) <= h_l)).sum()),
                            int(((Xc >= 0) & (np.abs(Xc) <= h_r)).sum())]])
        return RdbwhteResult(
            h=h_mat,
            W_lev=["T", "T:covs.hte"], W_names=["T", "T:covs.hte"],
            covs_hte_chr="covs_hte",
            kernel=Kernel, bwselect=bwselect, vce=vce_label, vce_select=vce,
            c=c, p=p, q=q, N=int(len(y)), Nh=nh_mat,
            covs_cont=True,
            rdmodel="MSE-Optimal BW for RD HTE (continuous covs.hte).",
        )

    # Categorical: one row per level. bw.joint=True -> repeat the global bw
    # across all levels (R bwjoint=TRUE).
    W_cat = W_cat_pre
    W_levels = list(W_cat.categories)
    n_lev = len(W_levels)
    h_mat  = np.zeros((n_lev, 2))
    nh_mat = np.zeros((n_lev, 2), dtype=int)

    if bwjoint:
        h_l_g, h_r_g = _bw_once()
        for j, lev in enumerate(W_levels):
            ind = np.asarray(W_cat) == lev
            h_mat[j, :] = [h_l_g, h_r_g]
            nh_mat[j, 0] = int(((Xc[ind] < 0)  & (np.abs(Xc[ind]) <= h_l_g)).sum())
            nh_mat[j, 1] = int(((Xc[ind] >= 0) & (np.abs(Xc[ind]) <= h_r_g)).sum())
    else:
        for j, lev in enumerate(W_levels):
            ind = np.asarray(W_cat) == lev
            try:
                h_l_j, h_r_j = _bw_once(sub_idx=ind)
            except Exception as e:
                warnings.warn(f"rdbwselect failed for level {lev!r}: {e}")
                h_l_j = h_r_j = np.nan
            h_mat[j, :] = [h_l_j, h_r_j]
            nh_mat[j, 0] = int(((Xc[ind] < 0)  & (np.abs(Xc[ind]) <= h_l_j)).sum())
            nh_mat[j, 1] = int(((Xc[ind] >= 0) & (np.abs(Xc[ind]) <= h_r_j)).sum())

    return RdbwhteResult(
        h=h_mat,
        W_lev=[str(lv) for lv in W_levels],
        W_names=[f"covs_hte{lv}" for lv in W_levels],
        covs_hte_chr="covs_hte",
        kernel=Kernel, bwselect=bwselect, vce=vce_label, vce_select=vce,
        c=c, p=p, q=q, N=int(len(y)), Nh=nh_mat,
        covs_cont=False,
        rdmodel="MSE-Optimal BW for RD HTE (categorical covs.hte).",
    )

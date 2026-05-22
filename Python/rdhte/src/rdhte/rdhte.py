"""rdhte: main estimation function.

Main estimation routine for RD heterogeneous treatment effects. Algorithm:

  1. Compute ``Xc = x - c`` and ``T = 1{Xc >= 0}``.
  2. Resolve ``covs_hte`` into one of four design paths:
       - none (ATE-only),
       - categorical (one CATE per cell; auto-detected from dtype, or
         multi-factor interaction via ``[W1, W2]`` / DataFrame / 2-D ndarray),
       - continuous (single numeric column -> linear-in-W CATE),
       - formula string (e.g. ``"w_left * w_strength"`` with ``data=``;
         ``factor(X)`` is accepted as an alias for patsy ``C(X)``).
  3. Determine bandwidth(s) (manual via ``h`` / ``h_l`` / ``h_r``, or
     via ``rdrobust.rdbwselect`` with ``bwselect=`` default ``"mserd"``).
  4. Build the saturated kernel-weighted design:
        y ~ T * Xp * W                          (categorical / continuous)
        y ~ T * Xp * (col1 + col2 + ...)        (formula path)
     plus ``+ covs_eff + covs_eff:W`` when ``covs_eff`` is supplied
     (no T or Xp interaction; efficiency-only adjustment).
  5. Fit WLS twice (p-order for point estimates, q-order for bias correction).
  6. Per-cell extraction:
       - reference cell      : tau = coef(T)
       - non-reference cell k: tau = coef(T) + coef(T:W_k)
     Likewise for ``tau.bc``.
  7. SE via heteroskedastic (HC0/1/2/3) or cluster (CR1/2/3) sandwich,
     with small-sample correction.
"""

from __future__ import annotations

import warnings
from typing import Any, Optional, Union

import numpy as np
import pandas as pd
import statsmodels.api as sm

from ._helpers import (
    resolve_arg, resolve_covs_hte, combine_factors,
    resolve_formula_covs_hte, resolve_covs_eff,
)
from ._results import RdhteResult
from ._vce import _VALID_VCE, normalize_vce

ArrayLike = Union[np.ndarray, pd.Series, list]

_VALID_KERNEL = {"tri", "epa", "uni", "triangular", "epanechnikov", "uniform"}


def rdhte(
    y: ArrayLike,
    x: ArrayLike,
    c: float = 0.0,
    covs_hte: Optional[ArrayLike] = None,
    covs_eff: Optional[ArrayLike] = None,
    p: int = 1,
    q: Optional[int] = None,
    kernel: str = "tri",
    weights: Optional[ArrayLike] = None,
    h: Optional[Union[float, ArrayLike]] = None,
    h_l: Optional[Union[float, ArrayLike]] = None,
    h_r: Optional[Union[float, ArrayLike]] = None,
    bwselect: str = "mserd",
    vce: str = "hc3",
    cluster: Optional[ArrayLike] = None,
    level: float = 95.0,
    subset: Optional[ArrayLike] = None,
    bwjoint: bool = False,
    data: Optional[pd.DataFrame] = None,
) -> RdhteResult:
    """Heterogeneous treatment effects in a sharp RD design.

    Main estimator for conditional treatment effects in a sharp RD design.
    See the module docstring for the algorithm. The theoretical framework is
    Calonico, Cattaneo, Farrell, Palomba, and Titiunik (2025a), with a
    software introduction in Calonico, Cattaneo, Farrell, Palomba, and
    Titiunik (2025b). Accepted ``covs_hte`` forms:

    - **None** -> ATE-only (no heterogeneity).
    - **Categorical**: 1-D array / pandas Series / ``pd.Categorical``.
      Auto-detected from dtype (object/string/bool/category, or numeric
      with at most 10 unique values).
    - **Multi-factor interaction**: list of arrays, ``pd.DataFrame`` with
      multiple columns, or 2-D ``np.ndarray``. Cells are the Cartesian
      product with R's ``interaction(..., sep=".")`` labels.
    - **Continuous**: 1-D numeric array (linear-in-W CATE).
    - **Formula string**: e.g. ``"w_left * w_strength"`` (requires
      ``data=``). ``factor(X)`` is accepted as an alias for patsy ``C(X)``.

    Parameters
    ----------
    y, x : array-like
        Outcome and running variable (length-n).
    c : float
        RD cutoff (default 0).
    covs_hte : optional, see above for accepted forms.
    covs_eff : optional 1-D / 2-D array / DataFrame
        Efficiency-improving covariates. Enter additively (and as
        ``covs_eff:W`` interactions in W-aware paths); never interact
        with T or the running-variable polynomial.
    p, q : int
        Local-polynomial orders for point estimation (p) and bias
        correction (q). Defaults: ``p=1``; ``q=None`` then resolved to
        ``p + 1``.
    kernel : {"tri", "epa", "uni"} or {"triangular", "epanechnikov", "uniform"}
        Default ``"tri"``.
    weights : optional 1-D array
        Unit-specific weights that multiply the kernel.
    h, h_l, h_r : float or None
        Bandwidth(s). Supply ``h`` (same on both sides), or both ``h_l``
        and ``h_r``. If none of these is given, ``rdrobust.rdbwselect``
        chooses one using ``bwselect``.
    bwselect : str
        ``rdrobust.rdbwselect`` method. Default ``"mserd"`` (one common
        MSE-optimal bandwidth). Other accepted values: ``"msetwo"``,
        ``"msesum"``, ``"msecomb1"``, ``"msecomb2"``, ``"cerrd"``,
        ``"certwo"``, ``"cersum"``, ``"cercomb1"``, ``"cercomb2"``.
        MSE denotes mean square error; CER denotes coverage error rate.
    vce : str
        Variance estimator. Without ``cluster=``: ``"hc0"``-``"hc3"`` or
        ``"nn"``; default ``"hc3"``. With ``cluster=``: ``"cr1"`` (default
        when cluster supplied), ``"cr2"`` (Bell-McCaffrey), ``"cr3"``
        (block jackknife). Legacy ``hc*`` + cluster combinations are
        warned and remapped to the matching ``cr*``.
        The fitted object exposes both ``.vce`` (display label, e.g.
        ``"CR1"``) and ``.vce_select`` (canonical lowercase form,
        e.g. ``"cr1"``).
    cluster : optional array-like
        Cluster identifier.
    level : float
        Confidence level (percent). Default 95.
    subset : optional bool or int-index array.
    bwjoint : bool
        For categorical / 0-1 ``covs_hte``: when ``False`` (the default),
        run an independent
        ``rdbwselect`` per cell so the bandwidth is cell-optimal. Set
        ``True`` to force a single shared bandwidth across cells (~n
        times faster). Continuous and multi-column W always use a
        shared bandwidth (this argument is a no-op there).
    data : optional pandas DataFrame
        When supplied, ``y``, ``x``, ``covs_hte``, ``covs_eff``,
        ``cluster``, ``weights``, ``subset`` may be bare-name strings
        referring to columns. Formula-string ``covs_hte`` requires
        ``data=`` (variables are resolved against it via patsy).

    Returns
    -------
    RdhteResult
        Dataclass with estimates, inference results, bandwidths, labels,
        and metadata. Convenience methods:
        ``summary()`` (formatted text), ``tidy()`` (broom-style
        DataFrame), ``glance()`` (one-row metadata DataFrame).

    Examples
    --------
    >>> from rdhte import rdhte                    # doctest: +SKIP
    >>> m = rdhte(y=y, x=x, c=0.0, covs_hte=W, cluster=cl)   # doctest: +SKIP
    >>> print(m.summary())                          # doctest: +SKIP
    """
    # ---- 1. Resolve `data=` lookups and coerce to ndarray ------------------
    y        = np.asarray(_arr(resolve_arg(y,        data, "y")),        dtype=float)
    x        = np.asarray(_arr(resolve_arg(x,        data, "x")),        dtype=float)
    cluster  = _arr(resolve_arg(cluster,  data, "cluster"))
    weights  = _arr(resolve_arg(weights,  data, "weights"))
    subset   = _arr(resolve_arg(subset,   data, "subset"))
    covs_eff = resolve_arg(covs_eff, data, "covs_eff")

    n_orig = len(y)
    covs_eff_mat = resolve_covs_eff(covs_eff, n_orig)
    if len(x) != n_orig:
        raise ValueError(f"'y' and 'x' must have equal length (got y={len(y)}, x={len(x)}).")

    # Detect formula-string covs_hte and materialize the model matrix early so
    # the resulting columns flow through the same subset/NA pipeline as the
    # other inputs. Mirrors R rdhte's `is.character(covs.hte)` branch.
    #
    # BARE-COLUMN-NAME shortcut: if `data` is supplied and `covs_hte` is a
    # bare Python identifier matching a column of `data` (no formula
    # operators), look up the column and dispatch through the array path
    # so 0/1 binaries auto-promote to factor (matching R + Stata). Without
    # this, ``rdhte(..., data=df, covs_hte="w_left")`` would treat w_left
    # as a continuous patsy formula even when it's a 0/1 column.
    formula_W_mat = None       # set when covs_hte is a formula string
    formula_W_cols = None      # column names of the patsy-built matrix
    formula_W_chr  = None      # user-facing label
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
        formula_W_mat  = W_mat_df.values.astype(float, copy=False)
        formula_W_cols = list(W_mat_df.columns)
        formula_W_chr  = covs_hte
        covs_components = None
    else:
        covs_hte = resolve_arg(covs_hte, data, "covs_hte")
        # Split covs_hte into 1+ components (supports factor-by-factor interaction).
        covs_components = resolve_covs_hte(covs_hte, n_orig)

    # ---- 2. Subset (if requested) ------------------------------------------
    if subset is not None:
        if subset.dtype == bool:
            if len(subset) != n_orig:
                raise ValueError("Logical subset must have length equal to length(x).")
            mask = subset
        else:
            mask = np.zeros(n_orig, dtype=bool)
            mask[np.asarray(subset, dtype=int)] = True
        y = y[mask]; x = x[mask]
        if covs_components is not None:
            covs_components = [c[mask] for c in covs_components]
        if formula_W_mat is not None:
            formula_W_mat = formula_W_mat[mask]
        if covs_eff_mat is not None:
            covs_eff_mat = covs_eff_mat[mask]
        if cluster  is not None: cluster  = cluster[mask]
        if weights  is not None: weights  = weights[mask]

    # ---- 3. Drop NA on core columns ---------------------------------------
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
    if covs_eff_mat is not None:
        covs_eff_mat = covs_eff_mat[keep]
    if formula_W_mat is not None:
        formula_W_mat = formula_W_mat[keep]
    if cluster  is not None: cluster  = cluster[keep]
    if weights  is not None: weights  = weights[keep]
    # Build the resolved covs_hte (either a 1-D array passed through or a
    # combined Cartesian-product Categorical for multi-factor interaction).
    if covs_components is None:
        covs_hte = None
    elif len(covs_components) == 1:
        covs_hte = np.asarray(covs_components[0])
    else:
        # Multi-component interaction: require all components to look factor-like.
        for j, comp in enumerate(covs_components):
            s = pd.Series(np.asarray(comp))
            if not _looks_categorical(s):
                raise ValueError(
                    f"Multi-component covs_hte requires all components to be "
                    f"categorical / factor-like (component {j} looks continuous). "
                    f"Mixed continuous-and-factor interactions require the formula-string "
                    f"covs_hte path, which is not yet supported in this release."
                )
        covs_hte = combine_factors(covs_components)

    # ---- 4. Polynomial orders + kernel + cutoff ---------------------------
    if q is None:
        q = p + 1
    if p < 0 or q <= p:
        raise ValueError(f"Need q > p; got p={p}, q={q}.")
    kernel = kernel.lower()
    if kernel not in _VALID_KERNEL:
        raise ValueError(f"Unsupported kernel {kernel!r}.")
    kernel = {"triangular": "tri", "epanechnikov": "epa", "uniform": "uni"}.get(kernel, kernel)

    Xc = x - c
    T  = (Xc >= 0).astype(float)

    # ---- 5. vce + cluster normalization. Shared with rdbwhte() via
    # _vce.normalize_vce; see that module for the full remapping rules.
    # The "user_vce" proxy treats anything other than the default "hc3"
    # as user-supplied, which means cluster+hc3 takes the silent-default
    # path to cr1 (rather than warning about hc3+cluster).
    user_vce = vce is not None and vce != "hc3"
    vce, vce_label = normalize_vce(vce, cluster, user_vce)
    # Internal vce code passed to statsmodels' get_robustcov_results
    vce_internal = {"cr1": "HC1", "cr2": "HC2", "cr3": "HC3"}.get(vce, vce.upper())

    # ---- 5b. Detect categorical W EARLY so the bw block can decide
    # whether per-cell rdbwselect (R-style bw.joint=FALSE) applies. We
    # only need to know `is_factor` and the cell levels here; the actual
    # downstream dispatch happens later.
    if formula_W_mat is not None or covs_hte is None:
        _W_is_factor = False
        _W_cat = None
        _W_levels = None
    elif isinstance(covs_hte, pd.Categorical):
        _W_is_factor = True
        _W_cat = covs_hte
        _W_levels = list(_W_cat.categories)
    else:
        _W_arr = pd.Series(np.asarray(covs_hte))
        _W_is_factor = _looks_categorical(_W_arr)
        _W_cat = pd.Categorical(_W_arr) if _W_is_factor else None
        _W_levels = list(_W_cat.categories) if _W_is_factor else None

    # ---- 6. Bandwidth --------------------------------------------------------
    # Four modes (mirror R rdhte):
    #   (a) h supplied scalar              -> h_l = h_r = h
    #   (b) h_l and h_r supplied (scalar)  -> use directly
    #   (c) auto bw + categorical W + bwjoint=False (DEFAULT, matches R
    #       bw.joint=FALSE): per-cell rdbwselect, h_lev is I x 2
    #   (d) auto bw + (continuous W OR bwjoint=True OR no covs_hte):
    #       single pooled rdbwselect, h_lev broadcast across cells
    h_lev = None  # I x 2 of per-cell bws when in mode (c); else None
    if h is not None and h_l is None and h_r is None:
        h_l = float(h); h_r = float(h)
        bwselect_label = "Manual"
    elif h_l is not None and h_r is not None:
        h_l = float(h_l); h_r = float(h_r)
        bwselect_label = "Manual"
    elif h is None and h_l is None and h_r is None:
        from rdrobust import rdbwselect
        # R rdhte forwards its own `vce` to rdbwselect (so the bw is chosen
        # consistently with the inference variance). rdrobust accepts
        # "nn"/"hc0"-"hc3" only (not the cr1/cr2/cr3 aliases); when our user
        # asks for cluster vce we pass the matching hc* form.
        bw_vce = {"cr1": "hc1", "cr2": "hc2", "cr3": "hc3"}.get(vce, vce)
        if bw_vce not in {"nn", "hc0", "hc1", "hc2", "hc3"}:
            bw_vce = "hc3"
        do_per_cell = (_W_is_factor and not bwjoint and _W_levels is not None
                       and len(_W_levels) >= 2)
        if do_per_cell:
            # Per-cell rdbwselect (mode c). Mirrors R rdhte's bw.joint=FALSE
            # branch (rdhte.R:514-525).
            W_codes = np.asarray(_W_cat.codes)
            h_lev = np.zeros((len(_W_levels), 2))
            for l_idx in range(len(_W_levels)):
                ind = (W_codes == l_idx)
                bw_l = rdbwselect(
                    y=y, x=Xc, c=0.0, p=p, q=q, kernel=kernel,
                    bwselect=bwselect, vce=bw_vce,
                    cluster=cluster if cluster is not None else None,
                    covs=covs_eff_mat if covs_eff_mat is not None else None,
                    subset=ind,
                )
                h_lev[l_idx, 0] = float(bw_l.bws.iloc[0, 0])
                h_lev[l_idx, 1] = float(bw_l.bws.iloc[0, 1])
            # Scalar h_l / h_r are the max so the (Xc, kw) downstream tools
            # have valid sentinel values; the per-obs h_vec built below uses
            # the cell-specific entries instead.
            h_l = float(h_lev[:, 0].max())
            h_r = float(h_lev[:, 1].max())
        else:
            # Pooled rdbwselect (mode d).
            bw = rdbwselect(
                y=y, x=Xc, c=0.0, p=p, q=q, kernel=kernel,
                bwselect=bwselect, vce=bw_vce,
                cluster=cluster if cluster is not None else None,
                covs=covs_eff_mat if covs_eff_mat is not None else None,
            )
            bws_df = bw.bws
            h_l = float(bws_df.iloc[0, 0])
            h_r = float(bws_df.iloc[0, 1])
        bwselect_label = bwselect
    else:
        raise ValueError("Provide either `h`, or both `h_l` and `h_r`, or neither (for auto bw).")

    # When per-cell h_lev is set, build h_vec from per-cell entries so
    # both r_bw and kw are correct cell-by-cell. Otherwise broadcast the
    # scalar h_l / h_r as before.
    if h_lev is not None:
        W_codes = np.asarray(_W_cat.codes)
        h_left_vec  = h_lev[W_codes, 0]
        h_right_vec = h_lev[W_codes, 1]
        in_bw_left  = (Xc < 0)  & (np.abs(Xc) <= h_left_vec)
        in_bw_right = (Xc >= 0) & (np.abs(Xc) <= h_right_vec)
    else:
        in_bw_left  = (Xc < 0)  & (np.abs(Xc) <= h_l)
        in_bw_right = (Xc >= 0) & (np.abs(Xc) <= h_r)
    r_bw = in_bw_left | in_bw_right

    # ---- 7. Kernel weights ---------------------------------------------------
    # When h_lev is set (per-cell bw), the per-obs h_vec already came from
    # h_left_vec / h_right_vec above; reuse those instead of broadcasting.
    if h_lev is not None:
        h_vec = np.where(Xc < 0, h_left_vec, h_right_vec)
    else:
        h_vec = np.where(Xc < 0, h_l, h_r)
    if kernel == "tri":
        kw = np.maximum(1 - np.abs(Xc) / h_vec, 0.0)
        Kernel = "Triangular"
    elif kernel == "epa":
        kw = np.maximum(1 - (Xc / h_vec) ** 2, 0.0)
        Kernel = "Epanechnikov"
    else:
        kw = (np.abs(Xc) <= h_vec).astype(float)
        Kernel = "Uniform"
    if weights is not None:
        kw = kw * np.asarray(weights, dtype=float)

    # ---- 8. Build the W factor (categorical heterogeneity) ------------------
    if formula_W_mat is not None:
        return _fit_W_matrix(
            y=y, Xc=Xc, T=T, kw=kw, r_bw=r_bw,
            W_mat=formula_W_mat, W_cols=formula_W_cols,
            p=p, q=q, cluster=cluster, vce_internal=vce_internal,
            level=level, h_l=h_l, h_r=h_r,
            kernel=kernel, Kernel=Kernel, vce_label=vce_label,
            bwselect_label=bwselect_label, c=c, keep_mask=keep,
            in_bw_left=in_bw_left, in_bw_right=in_bw_right,
            covs_hte_chr=formula_W_chr,
            covs_eff_mat=covs_eff_mat,
        )
    if covs_hte is None:
        # No heterogeneity: just the average ATE.
        return _fit_no_hte(y, Xc, T, kw, r_bw, p, q, cluster, vce_internal,
                           level, h_l, h_r, kernel, Kernel, vce_label,
                           bwselect_label, c, covs_eff_mat=covs_eff_mat)

    # Decide categorical vs continuous mode for covs_hte.
    if isinstance(covs_hte, pd.Categorical):
        W_cat = covs_hte
        is_factor = True
    else:
        W_arr = pd.Series(np.asarray(covs_hte))
        is_factor = _looks_categorical(W_arr)
        W_cat = pd.Categorical(W_arr) if is_factor else None
    if not is_factor:
        return _fit_continuous_W(
            y=y, Xc=Xc, T=T, kw=kw, r_bw=r_bw, W_raw=np.asarray(covs_hte, dtype=float),
            p=p, q=q, cluster=cluster, vce_internal=vce_internal,
            level=level, h_l=h_l, h_r=h_r,
            kernel=kernel, Kernel=Kernel, vce_label=vce_label,
            bwselect_label=bwselect_label, c=c, keep_mask=keep, in_bw_left=in_bw_left,
            in_bw_right=in_bw_right, covs_eff_mat=covs_eff_mat,
        )
    W_levels = list(W_cat.categories)
    n_lev = len(W_levels)

    # ---- 9. Saturated joint design matrix on the within-bandwidth subsample
    yb = y[r_bw]
    Xc_b = Xc[r_bw]
    T_b  = T[r_bw]
    W_b  = W_cat[r_bw]
    kw_b = kw[r_bw]
    cl_b = cluster[r_bw] if cluster is not None else None
    Z_b  = covs_eff_mat[r_bw, :] if covs_eff_mat is not None else None
    n_eff = len(yb)

    def _design(poly_order: int) -> pd.DataFrame:
        """Saturated design: 1 + T + Xp + W + T:Xp + T:W + Xp:W + T:Xp:W.

        With ``covs_eff`` supplied, also adds ``covs_eff_<k>`` main effects
        and their interactions with each non-reference W indicator
        (mirrors R rdhte: ``y ~ T*Xp*W.covs + covs*W.covs``).

        Columns named to mirror R's lm() output so we can extract by name.
        """
        d: dict[str, np.ndarray] = {"Intercept": np.ones(n_eff)}
        d["T"] = T_b
        for k in range(1, poly_order + 1):
            d[f"Xc{k}"] = Xc_b ** k
            d[f"T:Xc{k}"] = T_b * (Xc_b ** k)
        # W indicators for non-reference levels.
        W_indices = np.asarray(W_b.codes)
        for j, lev in enumerate(W_levels[1:], start=1):
            ind_j = (W_indices == j).astype(float)
            d[f"W.covs{lev}"]      = ind_j
            d[f"T:W.covs{lev}"]    = T_b * ind_j
            for k in range(1, poly_order + 1):
                d[f"Xc{k}:W.covs{lev}"]   = (Xc_b ** k) * ind_j
                d[f"T:Xc{k}:W.covs{lev}"] = T_b * (Xc_b ** k) * ind_j
        if Z_b is not None:
            for zk in range(Z_b.shape[1]):
                col = Z_b[:, zk]
                d[f"covs_eff_{zk}"] = col
                # Z * W interactions (non-reference levels only), no T or Xp.
                for j, lev in enumerate(W_levels[1:], start=1):
                    ind_j = (W_indices == j).astype(float)
                    d[f"covs_eff_{zk}:W.covs{lev}"] = col * ind_j
        return pd.DataFrame(d)

    Dp = _design(p)
    Dq = _design(q)

    # ---- 10. Weighted OLS for p (point) and q (bias-corrected) --------------
    # statsmodels.OLS supports weights via WLS; equivalently we can scale.
    Wsq = np.sqrt(kw_b)
    drop = Wsq <= 0
    if drop.any():
        yb_w   = yb[~drop]
        Dp_w   = Dp.loc[~drop, :]
        Dq_w   = Dq.loc[~drop, :]
        Wsq_w  = Wsq[~drop]
        cl_b_w = cl_b[~drop] if cl_b is not None else None
    else:
        yb_w, Dp_w, Dq_w, Wsq_w, cl_b_w = yb, Dp, Dq, Wsq, cl_b

    # Use WLS rather than scaling so that get_robustcov_results works directly.
    wls_p = sm.WLS(yb_w, Dp_w.values, weights=Wsq_w ** 2).fit()
    wls_q = sm.WLS(yb_w, Dq_w.values, weights=Wsq_w ** 2).fit()

    # Apply cluster / HC cov.
    if cluster is not None:
        cov_p = wls_p.get_robustcov_results(cov_type="cluster",
                                            groups=cl_b_w, use_correction=True,
                                            use_t=False)
        cov_q = wls_q.get_robustcov_results(cov_type="cluster",
                                            groups=cl_b_w, use_correction=True,
                                            use_t=False)
    else:
        cov_p = wls_p.get_robustcov_results(cov_type=vce_internal)
        cov_q = wls_q.get_robustcov_results(cov_type=vce_internal)

    coef_p = pd.Series(np.asarray(cov_p.params), index=Dp_w.columns)
    coef_q = pd.Series(np.asarray(cov_q.params), index=Dq_w.columns)
    V_q    = np.asarray(cov_q.cov_params())

    # ---- 11. Extract per-cell CATEs + SEs (mirrors R lines 666-740) ---------
    tau_hat    = np.full(n_lev, np.nan)
    tau_hat_bc = np.full(n_lev, np.nan)
    tau_hat_se = np.full(n_lev, np.nan)
    cell_names = [_label_lev(lev) for lev in W_levels]

    # Reference level CATE: T coefficient.
    idx_T_q = list(Dq_w.columns).index("T")
    tau_hat[0]    = coef_p["T"]
    tau_hat_bc[0] = coef_q["T"]
    tau_hat_se[0] = np.sqrt(V_q[idx_T_q, idx_T_q])

    # Non-reference levels: T + T:W_k.
    vcov = np.full((n_lev, n_lev), np.nan)
    vcov[0, 0] = tau_hat_se[0] ** 2
    for j, lev in enumerate(W_levels[1:], start=1):
        coef_name = f"T:W.covs{lev}"
        if coef_name not in coef_q.index:
            # Should not happen if design built correctly.
            continue
        idx_jk_q = list(Dq_w.columns).index(coef_name)
        tau_hat[j]    = coef_p["T"] + coef_p[coef_name]
        tau_hat_bc[j] = coef_q["T"] + coef_q[coef_name]
        var_j = (V_q[idx_T_q, idx_T_q]
                 + V_q[idx_jk_q, idx_jk_q]
                 + 2 * V_q[idx_T_q, idx_jk_q])
        tau_hat_se[j] = np.sqrt(max(var_j, 0.0))
        vcov[j, j] = var_j
        # Off-diagonal (cell 0 vs cell j): cov(T, T + T:W_j) = var(T) + cov(T, T:W_j)
        cov_0j = V_q[idx_T_q, idx_T_q] + V_q[idx_T_q, idx_jk_q]
        vcov[0, j] = vcov[j, 0] = cov_0j

    # Off-diagonal (cell i vs cell j), i, j >= 1:
    for i in range(1, n_lev):
        for j in range(i + 1, n_lev):
            name_i = f"T:W.covs{W_levels[i]}"
            name_j = f"T:W.covs{W_levels[j]}"
            idx_i = list(Dq_w.columns).index(name_i)
            idx_j = list(Dq_w.columns).index(name_j)
            cov_ij = (V_q[idx_T_q, idx_T_q]
                      + V_q[idx_T_q, idx_i]
                      + V_q[idx_T_q, idx_j]
                      + V_q[idx_i, idx_j])
            vcov[i, j] = vcov[j, i] = cov_ij

    # ---- 12. Robust BC z-stat, p-value, CIs --------------------------------
    qz   = _norm_ppf(1 - (1 - level/100) / 2)  # +1.96 at 95% (mirrors R: -qnorm(alpha/2))
    t_rb = tau_hat_bc / tau_hat_se
    pv_rb = 2 * _norm_cdf(-np.abs(t_rb))
    ci_rb = np.column_stack([tau_hat_bc - qz * tau_hat_se,
                             tau_hat_bc + qz * tau_hat_se])

    estimate_df = pd.DataFrame({
        "Group":       cell_names,
        "Estimate":    tau_hat,
        "Estimate.bc": tau_hat_bc,
        "se.rb":       tau_hat_se,
        "z.rb":        t_rb,
        "p.rb":        pv_rb,
        "ci.lo":       ci_rb[:, 0],
        "ci.hi":       ci_rb[:, 1],
    })

    rdmodel = "Sharp RD Heterogeneous Treatment Effects: Subgroups."
    if cluster is not None:
        rdmodel += " Cluster-Adjusted."

    # Build per-cell h and Nh matrices. When per-cell bw is in effect
    # (h_lev is not None), use the cell-specific entries; otherwise
    # broadcast the scalar (h_l, h_r) pair across all cells (joint bw).
    if h_lev is not None:
        h_out  = h_lev.copy()
        W_codes = np.asarray(W_cat.codes)
        Nh_out = np.zeros((n_lev, 2), dtype=int)
        for l_idx in range(n_lev):
            cell_mask = (W_codes == l_idx)
            Nh_out[l_idx, 0] = int((cell_mask & in_bw_left).sum())
            Nh_out[l_idx, 1] = int((cell_mask & in_bw_right).sum())
    else:
        h_out  = np.tile([[h_l, h_r]], (n_lev, 1))
        W_codes = np.asarray(W_cat.codes)
        Nh_out = np.zeros((n_lev, 2), dtype=int)
        for l_idx in range(n_lev):
            cell_mask = (W_codes == l_idx)
            Nh_out[l_idx, 0] = int((cell_mask & in_bw_left).sum())
            Nh_out[l_idx, 1] = int((cell_mask & in_bw_right).sum())

    return RdhteResult(
        Estimate=estimate_df,
        coef=pd.Series(tau_hat, index=cell_names),
        coef_bc=pd.Series(tau_hat_bc, index=cell_names),
        vcov=vcov,
        se_rb=tau_hat_se,
        ci_rb=ci_rb,
        t_rb=t_rb,
        pv_rb=pv_rb,
        coef_full=coef_p,
        vcov_full=V_q,
        W_lev=cell_names,
        W_names=[f"covs_hte{lev}" for lev in cell_names],
        covs_hte_chr="covs_hte",
        kernel=Kernel,
        bwselect=bwselect_label,
        vce=vce_label, vce_select=vce,
        c=c,
        h=h_out,
        p=p, q=q,
        N=int(keep.sum()),
        Nh=Nh_out,
        covs_cont=False,
        level=level,
        rdmodel=rdmodel,
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _arr(x):
    """Coerce to numpy; pass through None."""
    if x is None:
        return None
    if isinstance(x, pd.Series):
        return x.to_numpy()
    return np.asarray(x)


def _label_lev(lev) -> str:
    """Stringify a categorical level for use in coef names.

    Floats that happen to be whole numbers (e.g. ``0.0``, ``1.0``) are
    rendered as ``"0"``, ``"1"`` so the resulting coef names match the
    R ``factor(0/1)`` convention. Anything else falls back to ``str()``.
    """
    if isinstance(lev, float) and np.isfinite(lev) and lev.is_integer():
        return str(int(lev))
    # numpy float scalar (e.g. np.float64(0.0))
    try:
        f = float(lev)
        if np.isfinite(f) and f.is_integer() and not isinstance(lev, (str, bool)):
            # Only collapse to int when the *original* dtype is numeric float.
            # Booleans should stringify as "True"/"False" (caught above).
            if hasattr(lev, "dtype") and "float" in str(getattr(lev, "dtype", "")):
                return str(int(f))
    except (TypeError, ValueError):
        pass
    return str(lev)


def _looks_categorical(s: pd.Series) -> bool:
    """Heuristic: True if a Series should be treated as factor-like."""
    if isinstance(s.dtype, pd.CategoricalDtype):
        return True
    if s.dtype == object:
        return True
    if pd.api.types.is_string_dtype(s):
        return True
    if pd.api.types.is_bool_dtype(s):
        return True
    # Numeric with few unique values - treat as factor by default.
    if pd.api.types.is_numeric_dtype(s):
        return s.dropna().nunique() <= 10
    return False


def _norm_ppf(p):
    from scipy.stats import norm
    return norm.ppf(p)


def _norm_cdf(z):
    from scipy.stats import norm
    return norm.cdf(z)


def _fit_continuous_W(y, Xc, T, kw, r_bw, W_raw, p, q, cluster, vce_internal,
                      level, h_l, h_r, kernel, Kernel, vce_label,
                      bwselect_label, c, keep_mask, in_bw_left, in_bw_right,
                      covs_eff_mat=None) -> RdhteResult:
    """Continuous-W path. Mirrors R lines 673-697.

    Design: y ~ T * Xp * W  (W treated as continuous numeric).
    CATE function: tau(w) = coef(T) + coef(T:covs.hte) * w.
    Output: 2-row Estimate -- intercept (T) and slope (T:covs.hte).
    With ``covs_eff`` supplied: also adds main-effect ``covs_eff_<k>`` and
    ``covs_eff_<k>:W.covs`` columns (no T or Xp interaction).
    """
    yb   = y[r_bw]
    Xc_b = Xc[r_bw]
    T_b  = T[r_bw]
    W_b  = W_raw[r_bw]
    kw_b = kw[r_bw]
    cl_b = cluster[r_bw] if cluster is not None else None
    Z_b  = covs_eff_mat[r_bw, :] if covs_eff_mat is not None else None
    n_eff = len(yb)

    def _design(poly_order: int) -> pd.DataFrame:
        """Saturated y ~ T*Xp*W with W as a single continuous column."""
        d: dict[str, np.ndarray] = {"Intercept": np.ones(n_eff)}
        d["T"] = T_b
        for k in range(1, poly_order + 1):
            d[f"Xc{k}"] = Xc_b ** k
            d[f"T:Xc{k}"] = T_b * (Xc_b ** k)
        d["W.covs"]   = W_b
        d["T:W.covs"] = T_b * W_b
        for k in range(1, poly_order + 1):
            d[f"Xc{k}:W.covs"]   = (Xc_b ** k) * W_b
            d[f"T:Xc{k}:W.covs"] = T_b * (Xc_b ** k) * W_b
        if Z_b is not None:
            for zk in range(Z_b.shape[1]):
                col = Z_b[:, zk]
                d[f"covs_eff_{zk}"]        = col
                d[f"covs_eff_{zk}:W.covs"] = col * W_b
        return pd.DataFrame(d)

    Dp = _design(p)
    Dq = _design(q)
    Wsq = np.sqrt(kw_b)
    drop = Wsq <= 0
    if drop.any():
        yb_w   = yb[~drop]
        Dp_w   = Dp.loc[~drop, :]
        Dq_w   = Dq.loc[~drop, :]
        Wsq_w  = Wsq[~drop]
        cl_b_w = cl_b[~drop] if cl_b is not None else None
    else:
        yb_w, Dp_w, Dq_w, Wsq_w, cl_b_w = yb, Dp, Dq, Wsq, cl_b

    wls_p = sm.WLS(yb_w, Dp_w.values, weights=Wsq_w ** 2).fit()
    wls_q = sm.WLS(yb_w, Dq_w.values, weights=Wsq_w ** 2).fit()
    if cluster is not None:
        cov_p = wls_p.get_robustcov_results(cov_type="cluster",
                                            groups=cl_b_w, use_correction=True,
                                            use_t=False)
        cov_q = wls_q.get_robustcov_results(cov_type="cluster",
                                            groups=cl_b_w, use_correction=True,
                                            use_t=False)
    else:
        cov_p = wls_p.get_robustcov_results(cov_type=vce_internal)
        cov_q = wls_q.get_robustcov_results(cov_type=vce_internal)

    coef_p = pd.Series(np.asarray(cov_p.params), index=Dp_w.columns)
    coef_q = pd.Series(np.asarray(cov_q.params), index=Dq_w.columns)
    V_q    = np.asarray(cov_q.cov_params())

    # Extract T (intercept-at-W=0 CATE) and T:W.covs (slope).
    cell_names = ["T", "T:covs.hte"]   # matches R rdhte output convention
    coef_names = ["T", "T:W.covs"]
    tau_hat    = np.array([coef_p[c0] for c0 in coef_names])
    tau_hat_bc = np.array([coef_q[c0] for c0 in coef_names])
    idx        = [list(Dq_w.columns).index(c0) for c0 in coef_names]
    tau_hat_se = np.array([np.sqrt(V_q[i, i]) for i in idx])
    vcov       = np.array([[V_q[i, j] for j in idx] for i in idx])

    qz   = _norm_ppf(1 - (1 - level / 100) / 2)
    t_rb = tau_hat_bc / tau_hat_se
    pv_rb = 2 * _norm_cdf(-np.abs(t_rb))
    ci_rb = np.column_stack([tau_hat_bc - qz * tau_hat_se,
                             tau_hat_bc + qz * tau_hat_se])

    est = pd.DataFrame({
        "Group":       cell_names,
        "Estimate":    tau_hat,
        "Estimate.bc": tau_hat_bc,
        "se.rb":       tau_hat_se,
        "z.rb":        t_rb,
        "p.rb":        pv_rb,
        "ci.lo":       ci_rb[:, 0],
        "ci.hi":       ci_rb[:, 1],
    })
    rdmodel = "Sharp RD Heterogeneous Treatment Effects: Continuous."
    if cluster is not None: rdmodel += " Cluster-Adjusted."

    return RdhteResult(
        Estimate=est,
        coef=pd.Series(tau_hat, index=cell_names),
        coef_bc=pd.Series(tau_hat_bc, index=cell_names),
        vcov=vcov, se_rb=tau_hat_se,
        ci_rb=ci_rb, t_rb=t_rb, pv_rb=pv_rb,
        coef_full=coef_p, vcov_full=V_q,
        W_lev=cell_names, W_names=cell_names,
        covs_hte_chr="covs_hte", kernel=Kernel,
        bwselect=bwselect_label, vce=vce_label, vce_select=vce_internal,
        c=c, h=np.array([[h_l, h_r]]), p=p, q=q,
        N=int(keep_mask.sum()),
        Nh=np.array([[int(in_bw_left.sum()), int(in_bw_right.sum())]]),
        covs_cont=True, level=level, rdmodel=rdmodel,
    )


def _fit_no_hte(y, Xc, T, kw, r_bw, p, q, cluster, vce_internal,
                level, h_l, h_r, kernel, Kernel, vce_label,
                bwselect_label, c, covs_eff_mat=None) -> RdhteResult:
    """ATE-only path (no covs_hte). Mirrors R lines ~330-350.

    With ``covs_eff`` supplied: adds ``covs_eff_<k>`` main-effect columns
    (additive, no T or Xp interaction). Mirrors R ``y ~ T*Xp + covs``.
    """
    yb = y[r_bw]
    Xc_b = Xc[r_bw]
    T_b  = T[r_bw]
    kw_b = kw[r_bw]
    cl_b = cluster[r_bw] if cluster is not None else None
    Z_b  = covs_eff_mat[r_bw, :] if covs_eff_mat is not None else None

    def _design_simple(poly_order):
        d = {"Intercept": np.ones_like(yb), "T": T_b}
        for k in range(1, poly_order + 1):
            d[f"Xc{k}"] = Xc_b ** k
            d[f"T:Xc{k}"] = T_b * (Xc_b ** k)
        if Z_b is not None:
            for zk in range(Z_b.shape[1]):
                d[f"covs_eff_{zk}"] = Z_b[:, zk]
        return pd.DataFrame(d)

    Dp = _design_simple(p)
    Dq = _design_simple(q)
    fit_p = sm.WLS(yb, Dp.values, weights=kw_b).fit()
    fit_q = sm.WLS(yb, Dq.values, weights=kw_b).fit()
    if cluster is not None:
        cov_q = fit_q.get_robustcov_results(cov_type="cluster", groups=cl_b,
                                            use_correction=True, use_t=False)
    else:
        cov_q = fit_q.get_robustcov_results(cov_type=vce_internal)
    coef_p = pd.Series(fit_p.params,             index=Dp.columns)
    coef_q = pd.Series(np.asarray(cov_q.params), index=Dq.columns)
    V_q    = np.asarray(cov_q.cov_params())
    idx_T  = list(Dq.columns).index("T")
    tau    = float(coef_p["T"])
    tau_bc = float(coef_q["T"])
    se_rb  = float(np.sqrt(V_q[idx_T, idx_T]))
    qz     = _norm_ppf(1 - (1 - level/100) / 2)
    z_rb   = tau_bc / se_rb
    pv     = float(2 * _norm_cdf(-abs(z_rb)))
    ci     = np.array([tau_bc - qz*se_rb, tau_bc + qz*se_rb])
    est = pd.DataFrame({"Group": ["ATE"], "Estimate": [tau], "Estimate.bc": [tau_bc],
                        "se.rb": [se_rb], "z.rb": [z_rb], "p.rb": [pv],
                        "ci.lo": [ci[0]], "ci.hi": [ci[1]]})
    rdmodel = "Sharp RD Average Treatment Effect."
    if cluster is not None: rdmodel += " Cluster-Adjusted."
    return RdhteResult(
        Estimate=est, coef=pd.Series([tau], index=["T"]),
        coef_bc=pd.Series([tau_bc], index=["T"]),
        vcov=np.array([[se_rb**2]]), se_rb=np.array([se_rb]),
        ci_rb=ci.reshape(1, 2), t_rb=np.array([z_rb]), pv_rb=np.array([pv]),
        coef_full=coef_p, vcov_full=V_q,
        W_lev=["T"], W_names=["T"], covs_hte_chr=None,
        kernel=Kernel, bwselect=bwselect_label, vce=vce_label, vce_select=vce_internal,
        c=c, h=np.array([[h_l, h_r]]), p=p, q=q,
        N=int(len(y)),
        Nh=np.array([[int(((Xc < 0) & (np.abs(Xc) <= h_l)).sum()),
                      int(((Xc >= 0) & (np.abs(Xc) <= h_r)).sum())]]),
        covs_cont=False, level=level, rdmodel=rdmodel,
    )


def _fit_W_matrix(y, Xc, T, kw, r_bw, W_mat, W_cols, p, q, cluster, vce_internal,
                  level, h_l, h_r, kernel, Kernel, vce_label,
                  bwselect_label, c, keep_mask, in_bw_left, in_bw_right,
                  covs_hte_chr, covs_eff_mat=None) -> RdhteResult:
    """Multi-column continuous-W path (formula-string covs_hte).

    Design: ``y ~ T * Xp * (W_1 + W_2 + ... + W_K)``.
    CATE function: ``tau(W) = coef(T) + sum_k coef(T:W_k) * W_k``.
    Output: a row for ``T`` plus one row per column of ``W_mat``.
    Mirrors R rdhte's behavior when ``covs.hte`` is a formula string.
    With ``covs_eff`` supplied: adds ``covs_eff_<k>`` and
    ``covs_eff_<k>:W.<col>`` (no T or Xp interaction).
    """
    yb   = y[r_bw]
    Xc_b = Xc[r_bw]
    T_b  = T[r_bw]
    Wb   = W_mat[r_bw, :]
    kw_b = kw[r_bw]
    cl_b = cluster[r_bw] if cluster is not None else None
    Z_b  = covs_eff_mat[r_bw, :] if covs_eff_mat is not None else None
    n_eff, n_cols = Wb.shape

    # Replace ':' in formula-built names with '.' so rdhte_lincom backtick
    # syntax works without escaping (mirrors the R rdhte_lincom convention).
    safe_cols = [str(c0).replace(":", ".") for c0 in W_cols]

    def _design(poly_order: int) -> pd.DataFrame:
        d: dict[str, np.ndarray] = {"Intercept": np.ones(n_eff)}
        d["T"] = T_b
        for k in range(1, poly_order + 1):
            d[f"Xc{k}"] = Xc_b ** k
            d[f"T:Xc{k}"] = T_b * (Xc_b ** k)
        for idx, name in enumerate(safe_cols):
            wcol = Wb[:, idx]
            d[f"W.{name}"]   = wcol
            d[f"T:W.{name}"] = T_b * wcol
            for k in range(1, poly_order + 1):
                d[f"Xc{k}:W.{name}"]   = (Xc_b ** k) * wcol
                d[f"T:Xc{k}:W.{name}"] = T_b * (Xc_b ** k) * wcol
        if Z_b is not None:
            for zk in range(Z_b.shape[1]):
                col = Z_b[:, zk]
                d[f"covs_eff_{zk}"] = col
                for name in safe_cols:
                    d[f"covs_eff_{zk}:W.{name}"] = col * Wb[:, safe_cols.index(name)]
        return pd.DataFrame(d)

    Dp = _design(p)
    Dq = _design(q)
    Wsq = np.sqrt(kw_b)
    drop = Wsq <= 0
    if drop.any():
        yb_w   = yb[~drop]
        Dp_w   = Dp.loc[~drop, :]
        Dq_w   = Dq.loc[~drop, :]
        Wsq_w  = Wsq[~drop]
        cl_b_w = cl_b[~drop] if cl_b is not None else None
    else:
        yb_w, Dp_w, Dq_w, Wsq_w, cl_b_w = yb, Dp, Dq, Wsq, cl_b

    wls_p = sm.WLS(yb_w, Dp_w.values, weights=Wsq_w ** 2).fit()
    wls_q = sm.WLS(yb_w, Dq_w.values, weights=Wsq_w ** 2).fit()
    if cluster is not None:
        cov_p = wls_p.get_robustcov_results(cov_type="cluster",
                                            groups=cl_b_w, use_correction=True,
                                            use_t=False)
        cov_q = wls_q.get_robustcov_results(cov_type="cluster",
                                            groups=cl_b_w, use_correction=True,
                                            use_t=False)
    else:
        cov_p = wls_p.get_robustcov_results(cov_type=vce_internal)
        cov_q = wls_q.get_robustcov_results(cov_type=vce_internal)

    coef_p = pd.Series(np.asarray(cov_p.params), index=Dp_w.columns)
    coef_q = pd.Series(np.asarray(cov_q.params), index=Dq_w.columns)
    V_q    = np.asarray(cov_q.cov_params())

    # Extract: T (intercept), then T:W.<col> per matrix column.
    cell_names = ["T"] + [f"T:{nm}" for nm in safe_cols]
    coef_names = ["T"] + [f"T:W.{nm}" for nm in safe_cols]
    tau_hat    = np.array([coef_p[c0] for c0 in coef_names])
    tau_hat_bc = np.array([coef_q[c0] for c0 in coef_names])
    idx        = [list(Dq_w.columns).index(c0) for c0 in coef_names]
    tau_hat_se = np.array([np.sqrt(V_q[i, i]) for i in idx])
    vcov       = np.array([[V_q[i, j] for j in idx] for i in idx])

    qz   = _norm_ppf(1 - (1 - level / 100) / 2)
    t_rb = tau_hat_bc / tau_hat_se
    pv_rb = 2 * _norm_cdf(-np.abs(t_rb))
    ci_rb = np.column_stack([tau_hat_bc - qz * tau_hat_se,
                             tau_hat_bc + qz * tau_hat_se])

    est = pd.DataFrame({
        "Group":       cell_names,
        "Estimate":    tau_hat,
        "Estimate.bc": tau_hat_bc,
        "se.rb":       tau_hat_se,
        "z.rb":        t_rb,
        "p.rb":        pv_rb,
        "ci.lo":       ci_rb[:, 0],
        "ci.hi":       ci_rb[:, 1],
    })
    rdmodel = "Sharp RD Heterogeneous Treatment Effects: Formula."
    if cluster is not None:
        rdmodel += " Cluster-Adjusted."

    return RdhteResult(
        Estimate=est,
        coef=pd.Series(tau_hat, index=cell_names),
        coef_bc=pd.Series(tau_hat_bc, index=cell_names),
        vcov=vcov, se_rb=tau_hat_se,
        ci_rb=ci_rb, t_rb=t_rb, pv_rb=pv_rb,
        coef_full=coef_p, vcov_full=V_q,
        W_lev=cell_names, W_names=cell_names,
        covs_hte_chr=covs_hte_chr,
        kernel=Kernel, bwselect=bwselect_label, vce=vce_label, vce_select=vce_internal,
        c=c, h=np.array([[h_l, h_r]]), p=p, q=q,
        N=int(keep_mask.sum()),
        Nh=np.array([[int(in_bw_left.sum()), int(in_bw_right.sum())]]),
        covs_cont=True, level=level, rdmodel=rdmodel,
    )

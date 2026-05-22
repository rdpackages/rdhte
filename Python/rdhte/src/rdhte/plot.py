"""Plotting helpers for rdhte.

Categorical ``covs_hte`` only: draws one point per group at the
conventional point estimate (``Estimate``) with the robust bias-corrected
confidence interval (``ci.rb``). A dashed horizontal line at ``y = 0``
gives a visual reference for the null effect.

Default backend is `plotnine` (ggplot2 syntax). Returns the ggplot
object so users can compose further layers.
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import pandas as pd

from ._results import RdhteResult


def plot(
    obj: RdhteResult,
    *,
    sort: bool = False,
    point_size: float = 2.5,
    errorbar_width: float = 0.2,
    zero_line: bool = True,
    title: Optional[str] = None,
    xlab: Optional[str] = None,
    ylab: Optional[str] = None,
):
    """Forest plot of per-group CATEs with robust BC confidence intervals.

    Parameters
    ----------
    obj : RdhteResult
        Object returned by :func:`rdhte.rdhte`.
    sort : bool
        If ``True``, reorder groups along the x-axis by point estimate.
        Default ``False`` (preserve the original ``W_lev`` order).
    point_size : float
        Size of the point markers. Default ``2.5``.
    errorbar_width : float
        Width of the error-bar caps. Default ``0.2``.
    zero_line : bool
        If ``True`` (default), draw a dashed reference line at ``y = 0``.
    title, xlab, ylab : str or None
        Plot annotations. Defaults derived from ``obj`` metadata.

    Returns
    -------
    plotnine.ggplot
        Returned for further composition (``+ theme(...)`` etc.).

    Raises
    ------
    ImportError
        If plotnine is not installed.
    NotImplementedError
        If ``obj`` was fit with continuous ``covs_hte``.
    """
    if obj.covs_cont:
        raise NotImplementedError(
            "plot() supports only categorical covs_hte. The fitted rdhte object "
            "has continuous covs_hte; a CATE-curve plot is on the roadmap as a "
            "separate helper."
        )

    try:
        from plotnine import (
            ggplot, aes, geom_errorbar, geom_point, geom_hline,
            labs, theme_bw, theme, element_text, scale_x_discrete,
        )
    except ImportError as e:                                # pragma: no cover
        raise ImportError(
            "plot() requires the 'plotnine' package. Install it with "
            "`pip install plotnine`."
        ) from e

    if obj.coef is None or obj.ci_rb is None:
        raise ValueError("rdhte object is missing coef / ci_rb; cannot plot.")

    group_lab = list(obj.W_lev) if obj.W_lev else list(map(str, obj.coef.index))
    est = obj.coef.to_numpy()
    ci_lo = obj.ci_rb[:, 0]
    ci_hi = obj.ci_rb[:, 1]

    df = pd.DataFrame({
        "group":    group_lab,
        "estimate": est,
        "ci_low":   ci_lo,
        "ci_high":  ci_hi,
    })

    if sort:
        df = df.sort_values("estimate").reset_index(drop=True)

    # Preserve x-axis order using a categorical with explicit levels.
    df["group"] = pd.Categorical(df["group"], categories=df["group"].tolist(), ordered=True)

    if title is None:
        title = obj.rdmodel or "rdhte: heterogeneous treatment effects"
    if xlab is None:
        xlab = obj.covs_hte_chr or "Group"
    if ylab is None:
        ylab = "Treatment effect"

    max_label_len = max(len(str(g)) for g in df["group"])
    rotate_x = (len(df) > 6) or (max_label_len > 6)

    p = (ggplot(df, aes(x="group", y="estimate")))
    if zero_line:
        # Use a hex grey because plotnine accepts the matplotlib named-color
        # set, which does not include numbered greys.
        p = p + geom_hline(yintercept=0, linetype="dashed", color="#666666")
    p = (
        p
        + geom_errorbar(aes(ymin="ci_low", ymax="ci_high"), width=errorbar_width)
        + geom_point(size=point_size)
        + labs(title=title, x=xlab, y=ylab)
        + theme_bw()
    )
    if rotate_x:
        p = p + theme(axis_text_x=element_text(angle=45, hjust=1))

    return p

"""rdhte_lincom: robust bias-corrected linear-combination tests after rdhte.

The hypotheses reference per-cell CATEs by name (as in ``model.coef``), e.g.

    rdhte_lincom(m, linfct="`A` - `B` = 0")
    rdhte_lincom(m, linfct=["`A` - `B` = 0", "`A` = 0"])

The parsing of the symbolic hypotheses is delegated to
``patsy.DesignInfo.linear_constraint``.
"""

from __future__ import annotations

from typing import Sequence, Union

import numpy as np
import pandas as pd
from scipy.stats import chi2, norm

from ._results import RdhteResult


def rdhte_lincom(
    model: RdhteResult,
    linfct: Union[str, Sequence[str]],
    level: float = 95.0,
    digits: int = 3,
) -> dict:
    """Robust bias-corrected linear-combination tests after :func:`rdhte`.

    Parameters
    ----------
    model : RdhteResult
        Fitted object returned by :func:`rdhte`.
    linfct : str or sequence of str
        Symbolic linear hypotheses written against the entries of
        ``model.coef``. Use backticks for names with special characters.
        Each hypothesis is of the form ``"expr = 0"`` (an equals sign with
        a numeric right-hand side; constant defaults to 0 if omitted).
    level : float in [1, 100)
        Confidence level. Default 95.
    digits : int
        Number of decimal places for the formatted output. Default 3.

    Returns
    -------
    dict with keys ``"individual"`` (per-hypothesis test stats) and
    ``"joint"`` (joint Wald test of all restrictions simultaneously).
    More general linear restrictions can be supplied as multiple hypotheses
    in ``linfct``.
    """
    from patsy import DesignInfo

    if not (1.0 <= level < 100.0):
        raise ValueError(
            f"`level` must be in [1, 100), got {level!r}. "
            "Use the percentage form (e.g. 95), not the fraction form (0.95)."
        )

    if isinstance(linfct, str):
        linfct = [linfct]
    linfct = list(linfct)

    if model.coef is None or model.vcov is None:
        raise ValueError(
            "model.coef and model.vcov must both be populated; rdhte was not run."
        )

    coef_names = list(model.coef.index)
    coef_bc    = model.coef_bc.to_numpy() if model.coef_bc is not None else model.coef.to_numpy()
    coef_pe    = model.coef.to_numpy()
    V          = np.asarray(model.vcov)

    # Build a safe-name mapping: Patsy's linear_constraint parser doesn't
    # accept names with `:`, `.`, `(`, `)`, etc., so we substitute opaque
    # placeholders for each user-facing coef name. The R rdhte_lincom
    # convention is to backtick-quote names containing special chars; we
    # accept that convention and also accept bare names that happen to be
    # clean identifiers.
    safe_names = [f"__v{i}__" for i in range(len(coef_names))]
    name_map = dict(zip(coef_names, safe_names))

    def _rewrite(h: str) -> str:
        """Rewrite user hypothesis string into Patsy-friendly tokens.

        Strategy: longest-match-first substitution on the LHS only. Both
        backtick-quoted and bare occurrences of each coef name are
        replaced with the safe placeholder. The RHS of '=' (the constraint
        constant) is left untouched so that numeric literals like '0' or
        '0.0' aren't mistaken for digit-only coef names.
        """
        import re
        # Split at the FIRST '=' (or '==' / '!='); rewrite LHS only.
        m = re.search(r"\s*(==|!=|=)\s*", h)
        if m:
            lhs, rhs = h[:m.start()], h[m.end():]
            op = m.group(1)
        else:
            lhs, rhs, op = h, "", ""

        # Sort by length DESC so e.g. "T:covs.hte" is replaced before "T".
        for nm in sorted(coef_names, key=len, reverse=True):
            safe = name_map[nm]
            # Replace the backtick-quoted form first.
            lhs = lhs.replace(f"`{nm}`", safe)
            # And the bare form, but only if the name is a clean identifier
            # (otherwise the bare form would be ambiguous and the user must
            # backtick-quote it). Note: pure-digit names like "0"/"1" are
            # `isalnum()` True but match bare numeric literals too, so we
            # require at least one alpha character to be safe.
            stripped = nm.replace("_", "")
            if stripped.isalnum() and not stripped.isdigit():
                lhs = re.sub(rf"(?<![A-Za-z0-9_]){re.escape(nm)}(?![A-Za-z0-9_])", safe, lhs)
        return lhs + (f" {op} {rhs}" if op else "")

    di = DesignInfo(safe_names)
    L_blocks = []
    rhs_blocks = []
    for h in linfct:
        h_safe = _rewrite(h)
        lc = di.linear_constraint(h_safe)
        L_blocks.append(lc.coefs)
        rhs_blocks.append(lc.constants)
    L   = np.vstack(L_blocks)
    rhs = np.vstack(rhs_blocks).ravel()

    # Linear combinations of point and bc estimates.
    theta_pe = L @ coef_pe - rhs   # for the conventional "estimate" column
    theta_bc = L @ coef_bc - rhs   # for the bias-corrected test statistic

    LVL = L @ V @ L.T
    se  = np.sqrt(np.maximum(np.diag(LVL), 0.0))

    qz   = norm.ppf(1 - (1 - level/100) / 2)
    z    = np.where(se > 0, theta_bc / se, np.nan)
    pv   = 2.0 * norm.cdf(-np.abs(z))
    ci   = np.column_stack([theta_bc - qz * se, theta_bc + qz * se])

    individual = pd.DataFrame({
        "hypothesis": linfct,
        "estimate":   _round(theta_pe, digits),
        "z_stat":     _round(z, digits),
        "p_value":    _round(pv, digits),
        "conf.low":   _round(ci[:, 0], digits),
        "conf.high":  _round(ci[:, 1], digits),
    })

    # Joint Wald (chi-squared with df = number of restrictions).
    try:
        LVL_inv = np.linalg.inv(LVL)
    except np.linalg.LinAlgError:
        LVL_inv = np.linalg.pinv(LVL)
    wald   = float(theta_bc @ LVL_inv @ theta_bc)
    df_w   = int(L.shape[0])
    p_w    = float(chi2.sf(wald, df_w))
    joint  = pd.DataFrame({
        "statistic": [round(wald, digits)],
        "df":        [df_w],
        "p_value":   [round(p_w, digits)],
    })

    return {"individual": individual, "joint": joint}


def _round(arr, digits):
    return np.round(np.asarray(arr, dtype=float), digits)

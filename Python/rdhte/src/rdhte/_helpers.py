"""Private helpers for ``rdhte`` and ``rdbwhte``.

* ``resolve_arg``: bare-name lookups against a ``data=`` DataFrame.
* ``resolve_covs_hte``: split a ``covs_hte`` argument into one or more
  1-D component arrays (single factor, single continuous, multi-factor
  list / DataFrame / 2-D ndarray).
* ``combine_factors``: combine multiple factor components into one
  Cartesian-product ``pd.Categorical`` using dot-separated labels (first
  factor varies fastest; reference cell is the input's first level).
* ``resolve_formula_covs_hte``: patsy-build a model matrix from an
  formula string. ``factor(X)`` is accepted as an alias for ``C(X)``.
* ``resolve_covs_eff``: coerce ``covs_eff`` into a 2-D float matrix.
"""

from __future__ import annotations

from typing import Any, Optional, Union

import numpy as np
import pandas as pd


def resolve_arg(value: Any, data: Optional[pd.DataFrame], name: str) -> Any:
    """Resolve a possibly-string argument against a DataFrame.

    If ``data`` is None: return ``value`` unchanged (must already be array-like).
    If ``data`` is a DataFrame and ``value`` is a string: look up the column.
    Otherwise: return ``value`` (assumed already array-like).

    Phase 1 may extend this to also handle simple expression strings like
    ``"factor(W)"`` via ``patsy.dmatrix``.
    """
    if data is None or value is None:
        return value
    if isinstance(value, str):
        if value not in data.columns:
            raise KeyError(f"{name!r}: column {value!r} not found in `data`")
        return data[value]
    return value


def resolve_covs_eff(value: Any, n: int) -> Optional[np.ndarray]:
    """Coerce ``covs_eff`` into a 2-D float matrix of shape (n, K), or None.

    Accepted forms:
    - 1-D array/Series  -> (n, 1)
    - 2-D ndarray       -> (n, K) verbatim
    - ``pd.DataFrame``  -> (n, K) verbatim
    - list/tuple of arrays of length n -> (n, K)
    - None -> None
    """
    if value is None:
        return None
    if isinstance(value, pd.DataFrame):
        mat = value.values.astype(float, copy=False)
    elif isinstance(value, np.ndarray) and value.ndim == 2:
        mat = value.astype(float, copy=False)
    elif isinstance(value, (list, tuple)) and len(value) > 0 and not isinstance(value[0], (str, bytes, bool, int, float, np.bool_, np.number)):
        cols = [np.asarray(c, dtype=float) for c in value]
        mat = np.column_stack(cols)
    else:
        arr = np.asarray(value, dtype=float)
        mat = arr.reshape(-1, 1) if arr.ndim == 1 else arr
    if mat.shape[0] != n:
        raise ValueError(f"covs_eff has {mat.shape[0]} rows; expected {n} to match length(y).")
    return mat


def resolve_formula_covs_hte(formula: str, data: pd.DataFrame) -> pd.DataFrame:
    """Build a model matrix from a ``covs_hte`` formula string.

    Examples:

    - ``"W2 + W3"``                  -> cols: ``W2``, ``W3``
    - ``"w_left*w_strength"``        -> cols: ``w_left``, ``w_strength``, ``w_left:w_strength``
    - ``"factor(W1)*W3"``            -> cols: ``C(W1)[T.<lvl>]`` (one per non-ref level),
                                       ``W3``, and their pairwise interactions

    ``factor(X)`` is accepted as an alias for patsy's ``C(X)``. NaN values
    are preserved (not dropped)
    so the caller can union them into the usual NA mask.
    """
    import re
    import patsy

    if data is None:
        raise ValueError(
            "covs_hte was passed as a formula string but `data=` was not supplied. "
            "Provide a pandas DataFrame in `data=` so the formula can resolve variable names."
        )
    s = formula.strip()
    if s.startswith("~"):
        s = s.lstrip("~").strip()
    # factor(X) -> patsy C(X).
    s = re.sub(r"\bfactor\(", "C(", s)
    expr = f"~ ({s}) - 1"
    mat = patsy.dmatrix(
        expr, data=data, return_type="dataframe",
        NA_action=patsy.NAAction(NA_types=[]),
    )
    return mat


def _preserve_categorical(c: Any) -> Any:
    """Return a pd.Categorical if the input has categorical info; else ndarray."""
    if isinstance(c, pd.Categorical):
        return c
    if isinstance(c, pd.Series) and isinstance(c.dtype, pd.CategoricalDtype):
        return c.values  # underlying Categorical
    return np.asarray(c)


def resolve_covs_hte(value: Any, n: int) -> Optional[list]:
    """Split a ``covs_hte`` input into a list of 1-D components.

    Returns:
    - ``None`` if value is None.
    - ``[arr]`` (length 1) for a single 1-D input (the usual case).
    - ``[arr1, arr2, ...]`` for a multi-factor interaction input.

    Multi-factor input forms accepted:
    - ``pd.DataFrame`` with multiple columns
    - 2-D ``numpy.ndarray``
    - list/tuple whose first element is itself an array of length ``n``

    A list/tuple of scalars or strings is treated as a single 1-D factor
    (since that's the canonical Python way to write a categorical).

    Each returned component preserves its ``pd.Categorical`` identity so
    that downstream code (``combine_factors``) can keep the original
    category ordering. Non-categorical inputs become plain ndarrays.
    """
    if value is None:
        return None
    if isinstance(value, pd.DataFrame):
        return [_preserve_categorical(value[c]) for c in value.columns]
    if isinstance(value, np.ndarray) and value.ndim == 2:
        return [value[:, j] for j in range(value.shape[1])]
    if isinstance(value, (list, tuple)) and len(value) > 0:
        first = value[0]
        scalar_like = (
            isinstance(first, (str, bytes, bool, int, float, np.bool_, np.number))
            or np.ndim(first) == 0
        )
        if not scalar_like:
            elems = [_preserve_categorical(c) for c in value]
            if all(getattr(e, "ndim", 1) == 1 and len(e) == n for e in elems):
                return elems
            # Fall through: treat as a (badly-shaped) single 1-D factor.
    return [_preserve_categorical(value)]


def combine_factors(components: list) -> pd.Categorical:
    """Combine a list of 1-D factor arrays into one Categorical.

    Each row's label is the dot-separated join of the corresponding
    entries in each component.
    A single-element list returns the component as a Categorical directly.

    Category ordering uses the first component as the fastest-varying
    index, so the joint levels are
    ``[c1[0].c2[0], c1[1].c2[0], ..., c1[-1].c2[0], c1[0].c2[1], ...]``.
    Component category order is taken from the input ``Categorical``
    when supplied; otherwise from ``pd.unique`` (preserving first-seen
    order) on the raw values.
    """
    if len(components) == 1:
        comp = components[0]
        if isinstance(comp, pd.Categorical):
            return comp
        s = pd.Series(np.asarray(comp))
        return pd.Categorical(s, categories=list(pd.unique(s.dropna())))
    parts = []
    cats  = []
    for c in components:
        if isinstance(c, pd.Categorical):
            parts.append(pd.Series(c).astype(str))
            cats.append([str(lv) for lv in c.categories])
        else:
            s = pd.Series(np.asarray(c)).astype(str)
            parts.append(s)
            cats.append([str(lv) for lv in pd.unique(s.dropna())])
    joined = parts[0]
    for p in parts[1:]:
        joined = joined.str.cat(p, sep=".")
    # Mirror R: first factor varies fastest. Build the cross-product accordingly.
    joint_levels = list(cats[0])
    for cs in cats[1:]:
        joint_levels = [f"{a}.{b}" for b in cs for a in joint_levels]
    # Restrict to levels that actually appear (drops empty Cartesian cells).
    seen = set(pd.unique(joined.dropna()))
    joint_levels = [lv for lv in joint_levels if lv in seen]
    return pd.Categorical(joined, categories=joint_levels)

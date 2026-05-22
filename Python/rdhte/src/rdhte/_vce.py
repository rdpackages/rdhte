"""Shared vce normalization for ``rdhte`` and ``rdbwhte``.

Canonical names are ``hc0/hc1/hc2/hc3`` (no cluster) and ``cr1/cr2/cr3``
(with cluster), with documented warn-and-remap behavior for legacy
aliases. This module is the single source of truth for both estimators.

Downstream callers convert the canonical name to whatever their
particular consumer expects (``statsmodels.get_robustcov_results``
takes ``"HC0".."HC3"``; ``rdrobust.rdbwselect`` takes ``"nn"`` or
``"hc0".."hc3"`` only).
"""

from __future__ import annotations

import warnings

_VALID_VCE = {"nn", "hc0", "hc1", "hc2", "hc3", "cr1", "cr2", "cr3"}


def normalize_vce(vce: str, cluster, user_vce: bool):
    """Validate ``vce`` and apply the cluster<->no-cluster remappings.

    Parameters
    ----------
    vce : str
        The user-supplied vce string (case-insensitive).
    cluster : array-like or None
        Cluster vector; only its non-None-ness is used here.
    user_vce : bool
        ``True`` if the user supplied a non-default value of ``vce``.
        Required because the silent ``cr1`` default when cluster is
        supplied only applies if the user did not pass ``vce``
        explicitly. The caller must compute this with reference to the
        caller's own signature default.

    Returns
    -------
    (vce, vce_label) : tuple[str, str]
        Canonical vce name (``hc0``..``hc3`` or ``cr1``..``cr3``) and
        the matching uppercase display label.

    Raises
    ------
    ValueError
        If ``vce`` is not one of the accepted forms.
    """
    if not isinstance(vce, str):
        raise ValueError(f"vce must be a string; got {type(vce).__name__}.")
    vce = vce.lower()
    if vce not in _VALID_VCE:
        raise ValueError(f"vce {vce!r} not supported.")

    if cluster is not None:
        if not user_vce:
            vce = "cr1"  # silent default when cluster supplied
        elif vce in ("hc0", "hc1", "cr0"):
            warnings.warn(
                f"vce={vce!r} is not a cluster option. Switching to vce='cr1'."
            )
            vce = "cr1"
        elif vce == "hc2":
            warnings.warn(
                "vce='hc2' is not a cluster option. Switching to vce='cr2'."
            )
            vce = "cr2"
        elif vce == "hc3":
            warnings.warn(
                "vce='hc3' is not a cluster option. Switching to vce='cr3'."
            )
            vce = "cr3"
    else:
        if vce == "cr1":
            warnings.warn(
                "vce='cr1' requires a cluster variable. Falling back to vce='hc1'."
            )
            vce = "hc1"
        elif vce == "cr2":
            warnings.warn(
                "vce='cr2' requires a cluster variable. Falling back to vce='hc2'."
            )
            vce = "hc2"
        elif vce == "cr3":
            warnings.warn(
                "vce='cr3' requires a cluster variable. Falling back to vce='hc3'."
            )
            vce = "hc3"

    return vce, vce.upper()

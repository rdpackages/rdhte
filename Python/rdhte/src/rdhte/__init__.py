"""rdhte: Heterogeneous treatment effects in regression discontinuity designs.

Estimation, robust bias-corrected inference, automatic bandwidth selection,
cluster-robust standard errors, and linear-combination tests.

Public API
----------
- rdhte         : main estimation function.
- rdbwhte       : bandwidth-selection wrapper.
- rdhte_lincom  : linear-combination tests.
- RdhteResult   : dataclass returned by rdhte().
- RdbwhteResult : dataclass returned by rdbwhte().

The plotting helper :func:`rdhte.plot.plot` requires ``plotnine`` and is
imported lazily; the rest of the package has no plotting dependency.
"""

from __future__ import annotations

__version__ = "0.1.0"

from ._results import RdhteResult, RdbwhteResult
from .rdhte import rdhte
from .rdbwhte import rdbwhte
from .rdhte_lincom import rdhte_lincom

__all__ = [
    "__version__",
    "rdhte",
    "rdbwhte",
    "rdhte_lincom",
    "RdhteResult",
    "RdbwhteResult",
]

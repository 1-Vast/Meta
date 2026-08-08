"""Verified mathematical primitives and the P1B geometry bridge.

There is intentionally no assembled production DTA model. The former biological
frontend did not identify incremental affinity information and is recorded only
in ``history.md``. Unadmitted research algebra is not retained in this package.
"""

from .config import DEFAULT, MetaSieveConfig, profile
from .encoders import LigandEncoder, ProteinEncoder
from .mechanism import MechanisticInteractionBridge
from .meta_operator import BandOperator, CSMO, build_band_operator
from .runtime import configure_cuda, require_cuda

__all__ = [
    "BandOperator",
    "CSMO",
    "DEFAULT",
    "LigandEncoder",
    "MechanisticInteractionBridge",
    "MetaSieveConfig",
    "ProteinEncoder",
    "configure_cuda",
    "build_band_operator",
    "profile",
    "require_cuda",
]

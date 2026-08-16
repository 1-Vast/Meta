"""Explicit permission boundary for affinity training."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


STATE_MANIFEST = Path("manifests/state.v1.json")


def assertauthorized(manifest_path: Path = STATE_MANIFEST) -> dict[str, Any]:
    state = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not state["permissions"]["affinity_training_authorized"]:
        raise PermissionError(
            "HTL-DTA affinity training is not authorized; natural-tail Gate D0 "
            "is DATA_NOT_READY"
        )
    return state


def trainmeta(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Validate explicit training authorization before a runner consumes labels."""

    return assertauthorized()

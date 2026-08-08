"""Compatibility CLI for sealing canonical staging rows for model runtime."""
from __future__ import annotations

import csv
from pathlib import Path

from scripts.seal_compiled_dataset import main


def cluster_targets(_sequences, identity: float, verbose: bool):
    """Deprecated injection point retained for external split-map auditing."""
    raise RuntimeError("provide a homology cluster implementation before using this audit")


def audit_external_split(raw_path, split_map):
    """Reject a supplied map that separates one homology component."""
    with Path(raw_path).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    sequences = {str(row["target_id"]): str(row["sequence"]) for row in rows}
    missing = set(sequences) - set(split_map)
    extra = set(split_map) - set(sequences)
    if missing or extra:
        raise ValueError("external split-map target mismatch")
    components, _ = cluster_targets(sequences, identity=0.4, verbose=False)
    grouped = {}
    for target, component in components.items():
        grouped.setdefault(component, set()).add(split_map[target])
    straddling = sum(len(parts) > 1 for parts in grouped.values())
    if straddling:
        raise ValueError(f"external split map straddles {straddling} homology components")
    return {"valid": True, "homology_components_straddling": 0}


if __name__ == "__main__":
    raise SystemExit(main())

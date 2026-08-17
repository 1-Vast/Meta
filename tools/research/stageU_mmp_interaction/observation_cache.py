"""Deterministic disk cache for the U0 observation bank.

Building the MMP bank re-parses every ligand with RDKit MMPA (~60 s). The bank
is a pure function of structure and split metadata, so it is cached once as a
gzip JSON lines artifact. The cache is rebuilt only on `--force`; downstream
modules load it. No label chooses any cached field; the label is stored with
the row it describes and is never read during construction of keys, splits or
coverage.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
import sys

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from tools.research.stageU_mmp_interaction.observations import MMPObservation

HERE = Path(__file__).resolve().parent
CACHE = HERE / "U0_OBSERVATIONS.jsonl.gz"

FIELDS = (
    "target", "component", "core", "exact_key", "coarse_key", "cell_a",
    "cell_b", "ligand_a", "ligand_b", "delta_y", "same_panel", "stratum",
    "tanimoto", "activity_cliff", "stereo_edit", "charge_change", "r_a",
    "r_b", "context",
)


def save_observations(observations: list[MMPObservation], path: Path | None = None):
    destination = path or CACHE
    with gzip.open(destination, "wt", encoding="utf-8") as handle:
        for item in observations:
            payload = {field: getattr(item, field) for field in FIELDS}
            payload["context"] = list(item.context)
            handle.write(json.dumps(payload, sort_keys=True) + "\n")
    return destination


def load_observations(path: Path | None = None) -> list[MMPObservation]:
    source = path or CACHE
    out: list[MMPObservation] = []
    with gzip.open(source, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            row = json.loads(line)
            out.append(MMPObservation(
                target=row["target"], component=row["component"],
                core=row["core"], exact_key=row["exact_key"],
                coarse_key=row["coarse_key"], cell_a=row["cell_a"],
                cell_b=row["cell_b"], ligand_a=row["ligand_a"],
                ligand_b=row["ligand_b"], delta_y=float(row["delta_y"]),
                same_panel=bool(row["same_panel"]), stratum=row["stratum"],
                tanimoto=float(row["tanimoto"]),
                activity_cliff=bool(row["activity_cliff"]),
                stereo_edit=bool(row["stereo_edit"]),
                charge_change=int(row["charge_change"]),
                r_a=row["r_a"], r_b=row["r_b"], context=tuple(row["context"])))
    return out


def cache_sha256(path: Path | None = None) -> str:
    digest = hashlib.sha256()
    with (path or CACHE).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()

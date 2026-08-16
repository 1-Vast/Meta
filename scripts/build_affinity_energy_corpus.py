"""Build deterministic EnergyPilot.v1 rows from a verified ChEMBL37 SQLite."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.source_affinity.canonicalize import build_corpus


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument(
        "--sql",
        type=Path,
        default=ROOT / "contracts" / "source_affinity" / "chembl37_e0_core.sql",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "dataset" / "processed" / "source_affinity" / "energy_pilot_v1",
    )
    args = parser.parse_args()
    print(json.dumps(build_corpus(args.database, args.sql, args.output), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

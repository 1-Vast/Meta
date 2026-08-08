"""Acquire and verify an immutable source-affinity database release."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.source_affinity.chembl_static import acquire_release


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", choices=("chembl37",), default="chembl37")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "dataset" / "raw" / "source_affinity" / "chembl37_sqlite_v1",
    )
    parser.add_argument("--expected-bytes", type=int)
    args = parser.parse_args()
    result = acquire_release(args.output_dir, args.expected_bytes)
    print(json.dumps(result, sort_keys=True, indent=2))


if __name__ == "__main__":
    main()

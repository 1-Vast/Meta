"""Read-only verification of a compiled canonical DTA dataset."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from scripts.preprocess_dataset import audit_dataset


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    result = audit_dataset(args.directory)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())

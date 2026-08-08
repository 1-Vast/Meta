"""Command-line wrapper for :func:`test.smoke_check.run_smoke`."""
from __future__ import annotations

import argparse
import json

from test.smoke_check import run_smoke


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset", nargs="?", default="DAVIS")
    parser.add_argument("--profile", default="base")
    parser.add_argument("--sealed-dir", required=True)
    parser.add_argument("--protein-cache", required=True)
    args = parser.parse_args()
    print(json.dumps(run_smoke(args.dataset, profile_name=args.profile,
                               sealed_dir=args.sealed_dir,
                               protein_cache_path=args.protein_cache).to_dict(),
                     indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

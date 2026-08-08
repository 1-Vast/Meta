"""Runtime-bound leakage and model-invariant audit for a sealed source view."""
from __future__ import annotations

import argparse
import json

from scripts.data import load_episodes
from test.smoke_check import run_smoke


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset")
    parser.add_argument("--sealed-dir", required=True)
    parser.add_argument("--protein-cache", required=True)
    parser.add_argument("--profile", default="base")
    parser.add_argument("--seed", type=int, default=20260804)
    args = parser.parse_args()
    episodes = load_episodes(
        args.dataset, sealed_dir=args.sealed_dir, protein_cache_path=args.protein_cache,
    )
    if episodes.visible_splits != ("source",):
        raise AssertionError("audit must mount source labels only")
    result = {
        "seal_audit": episodes.seal_audit(),
        "smoke": run_smoke(
            args.dataset, profile_name=args.profile, seed=args.seed,
            sealed_dir=args.sealed_dir, protein_cache_path=args.protein_cache,
        ).to_dict(),
    }
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Remove undeclared formats from the three ESM-2 cache snapshots.

Dry-run is the default. Use --apply only when no bootstrap process is running.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from bootstrap_remote import MODEL_FILES, MODELS, ROOT


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()

    cache = (ROOT / "weights" / "hub").resolve()
    unwanted: list[tuple[Path, Path]] = []
    for model in MODELS:
        repository = cache / f"models--{model.replace('/', '--')}"
        for snapshot in (repository / "snapshots").glob("*"):
            for path in snapshot.iterdir():
                if path.name not in MODEL_FILES:
                    target = path.resolve(strict=False)
                    if not target.is_relative_to(cache):
                        raise RuntimeError(f"refusing target outside cache: {target}")
                    unwanted.append((path, target))

    for path, _ in unwanted:
        print(f"remove snapshot entry: {path.relative_to(ROOT)}")
        if args.apply:
            path.unlink(missing_ok=True)

    referenced = {
        path.resolve(strict=False)
        for path in cache.glob("models--*/snapshots/*/*")
        if path.is_symlink()
    }
    candidates = {target for _, target in unwanted}
    candidates.update(cache.glob("models--*/blobs/*.incomplete"))
    for path in sorted(candidates):
        if path not in referenced and path.is_file():
            print(f"remove unreferenced blob: {path.relative_to(ROOT)}")
            if args.apply:
                path.unlink()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

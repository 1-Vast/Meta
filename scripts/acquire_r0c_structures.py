"""Download only the frozen R0-C predownload panel from RCSB."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import json
from pathlib import Path

from scripts.data_contract import read_jsonl
from scripts.structure_sources.rcsb import (
    RCSB_FILE_SERVICE,
    WWPDB_LICENSE,
    atomic_download,
    ccd_url,
    mmcif_url,
    sha256_file,
)


FROZEN_PANEL_SHA256 = "754c4a9980e535d447a9da541d79b7aac4b8109af8c8b3e40bdbb3fc0d5757ee"


def acquire(
    panel_path: str | Path, output_root: str | Path, *, workers: int = 8
) -> dict:
    panel_file = Path(panel_path)
    root = Path(output_root)
    manifest_path = root / "acquisition_manifest.json"
    if manifest_path.exists():
        raise FileExistsError(f"R0-C acquisition manifest already exists: {manifest_path}")
    if workers < 1:
        raise ValueError("workers must be positive")
    if sha256_file(panel_file) != FROZEN_PANEL_SHA256:
        raise ValueError("R0-C predownload panel does not match the frozen SHA256")
    rows = read_jsonl(panel_file)
    if len(rows) != 512:
        raise ValueError("R0-C predownload panel must contain exactly 512 records")
    pdb_ids = sorted({str(row["pdb_id"]).lower() for row in rows})
    comp_ids = sorted({str(row["ligand_comp_id"]).upper() for row in rows})
    if len(pdb_ids) != len(rows) or len(comp_ids) != len(rows):
        raise ValueError("R0-C predownload PDB and ligand IDs must be one-to-one")
    jobs = [(mmcif_url(value), root / "mmcif" / f"{value}.cif.gz") for value in pdb_ids]
    jobs += [(ccd_url(value), root / "ccd" / f"{value}.cif") for value in comp_ids]
    completed, failed = [], []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(atomic_download, url, path): (url, path) for url, path in jobs}
        for future in as_completed(futures):
            url, path = futures[future]
            try:
                completed.append(future.result())
            except Exception as error:
                failed.append({"url": url, "path": str(path), "error": str(error)})
    manifest = {
        "schema": "MetaSieve.R0C.StructureAcquisition.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source": "BioLiP2_metadata+RCSB_PDB_CCD",
        "coordinate_service": RCSB_FILE_SERVICE,
        "coordinate_license": WWPDB_LICENSE,
        "panel": {"path": str(panel_file.resolve()), "sha256": sha256_file(panel_file)},
        "records": len(rows),
        "unique_pdb_ids": len(pdb_ids),
        "unique_ccd_ids": len(comp_ids),
        "files_ok": len(completed),
        "files_failed": len(failed),
        "files": sorted(completed, key=lambda item: item["path"]),
        "failures": sorted(failed, key=lambda item: item["path"]),
    }
    root.mkdir(parents=True, exist_ok=True)
    temporary = manifest_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    temporary.replace(manifest_path)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("panel")
    parser.add_argument("output_root")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    result = acquire(args.panel, args.output_root, workers=args.workers)
    print(json.dumps({key: result[key] for key in (
        "records", "unique_pdb_ids", "unique_ccd_ids", "files_ok", "files_failed"
    )}, indent=2))
    return 0 if not result["files_failed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())

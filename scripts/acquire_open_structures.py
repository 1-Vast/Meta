"""Acquire a bounded BioLiP2/RCSB pilot with resumable, hash-bound downloads."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
import json
from pathlib import Path

from scripts.structure_sources.biolip import (BIOLIP_ANNOTATION_URL,
    BIOLIP_DOWNLOAD_PAGE, BIOLIP_LIGAND_SUMMARY_URL, BIOLIP_PROTEIN_FASTA_URL,
    BIOLIP_README_URL, pilot_candidates, regular_ligand_ids)
from scripts.structure_sources.rcsb import (RCSB_FILE_SERVICE, WWPDB_LICENSE,
    atomic_download, ccd_url, mmcif_url)


def _download_many(jobs: list[tuple[str, Path]], workers: int) -> tuple[list[dict], list[dict]]:
    completed, failed = [], []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(atomic_download, url, path): (url, path)
                   for url, path in jobs}
        for future in as_completed(futures):
            url, path = futures[future]
            try:
                completed.append(future.result())
            except Exception as error:
                failed.append({"url": url, "path": str(path), "error": str(error)})
    return completed, failed


def acquire_open_structures(output_root: str | Path, *, candidate_limit: int,
                            workers: int = 8) -> dict:
    if candidate_limit < 1 or workers < 1:
        raise ValueError("candidate_limit and workers must be positive")
    root = Path(output_root)
    metadata = root / "biolip2"
    annotation = metadata / "BioLiP.txt.gz"
    fasta = metadata / "protein.fasta.gz"
    ligand_summary = metadata / "ligand.tsv.gz"
    metadata_results = [
        atomic_download(BIOLIP_ANNOTATION_URL, annotation),
        atomic_download(BIOLIP_PROTEIN_FASTA_URL, fasta),
        atomic_download(BIOLIP_LIGAND_SUMMARY_URL, ligand_summary),
    ]
    allowed_ligands = regular_ligand_ids(ligand_summary)
    candidates = pilot_candidates(annotation, limit=candidate_limit,
                                  allowed_ligands=allowed_ligands)
    pdb_ids = sorted({entry.pdb_id for entry in candidates})
    comp_ids = sorted({entry.ligand_comp_id for entry in candidates})
    jobs = [(mmcif_url(value), root / "mmcif" / f"{value}.cif.gz") for value in pdb_ids]
    jobs += [(ccd_url(value), root / "ccd" / f"{value}.cif") for value in comp_ids]
    downloaded, failed = _download_many(jobs, workers)
    manifest = {
        "schema": "MetaSieve.OpenStructureAcquisition.v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source": "BioLiP2+RCSB_PDB_CCD",
        "annotation_source": BIOLIP_DOWNLOAD_PAGE,
        "annotation_readme": BIOLIP_README_URL,
        "coordinate_service": RCSB_FILE_SERVICE,
        "coordinate_license": WWPDB_LICENSE,
        "candidate_limit": candidate_limit,
        "candidate_complexes": len(candidates),
        "prefilter_ligand_ids": len(allowed_ligands),
        "unique_pdb_ids": len(pdb_ids),
        "unique_ccd_ids": len(comp_ids),
        "metadata_files": metadata_results,
        "files_ok": len(downloaded),
        "files_failed": len(failed),
        "files": sorted(downloaded, key=lambda item: item["path"]),
        "failures": sorted(failed, key=lambda item: item["path"]),
    }
    root.mkdir(parents=True, exist_ok=True)
    path = root / "acquisition_manifest.json"
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    temporary.replace(path)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_root")
    parser.add_argument("--candidate-limit", type=int, default=15000)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    result = acquire_open_structures(args.output_root,
                                     candidate_limit=args.candidate_limit,
                                     workers=args.workers)
    print(json.dumps({key: result[key] for key in (
        "candidate_complexes", "unique_pdb_ids", "unique_ccd_ids",
        "files_ok", "files_failed")}, indent=2))
    return 0 if not result["files_failed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())

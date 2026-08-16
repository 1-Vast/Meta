"""RCSB PDB/CCD acquisition and canonical mmCIF chemistry helpers."""
from __future__ import annotations

import hashlib
import os
from pathlib import Path
import time

import requests


RCSB_MMCIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif.gz"
RCSB_CCD_URL = "https://files.rcsb.org/ligands/download/{comp_id}.cif"
RCSB_FILE_SERVICE = "https://www2.rcsb.org/docs/programmatic-access/file-download-services"
WWPDB_LICENSE = "CC0-1.0"


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_download(url: str, destination: str | Path, *, retries: int = 3,
                    timeout: float = 60.0) -> dict:
    """Download once, resume a partial file, then atomically publish it."""
    target = Path(destination)
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.is_file() and target.stat().st_size:
        return {"path": str(target), "bytes": target.stat().st_size,
                "sha256": sha256_file(target), "status": "cached", "url": url}
    partial = target.with_name(target.name + ".part")
    last_error: Exception | None = None
    for attempt in range(retries + 1):
        try:
            offset = partial.stat().st_size if partial.exists() else 0
            headers = {"Range": f"bytes={offset}-"} if offset else {}
            with requests.get(url, headers=headers, stream=True, timeout=timeout) as response:
                if offset and response.status_code == 200:
                    partial.unlink()
                    offset = 0
                response.raise_for_status()
                mode = "ab" if offset and response.status_code == 206 else "wb"
                with partial.open(mode) as handle:
                    for chunk in response.iter_content(1024 * 1024):
                        if chunk:
                            handle.write(chunk)
            if not partial.stat().st_size:
                raise ValueError(f"empty response from {url}")
            os.replace(partial, target)
            return {"path": str(target), "bytes": target.stat().st_size,
                    "sha256": sha256_file(target), "status": "downloaded", "url": url}
        except Exception as error:
            last_error = error
            if attempt < retries:
                time.sleep(min(2 ** attempt, 8))
    raise RuntimeError(f"download failed after {retries + 1} attempts: {url}: {last_error}")


def mmcif_url(pdb_id: str) -> str:
    return RCSB_MMCIF_URL.format(pdb_id=pdb_id.lower())


def ccd_url(comp_id: str) -> str:
    return RCSB_CCD_URL.format(comp_id=comp_id.upper())


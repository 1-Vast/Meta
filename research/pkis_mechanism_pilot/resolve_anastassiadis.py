"""Resolve Anastassiadis compound identifiers without reading assay outcomes.

Only workbook header rows 1--3 are read.  The cached output is an auditable
identity/structure sidecar used by the preregistered v2 transfer runner.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import time
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

import pandas as pd


PUBCHEM = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/"
    "property/ConnectivitySMILES,SMILES/JSON"
)


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _lookup(identifier: str, attempts: int = 5) -> dict:
    url = PUBCHEM.format(name=quote(identifier, safe=""))
    error = None
    for attempt in range(attempts):
        try:
            request = Request(url, headers={"User-Agent": "MetaSieve-external-pilot/2.0"})
            with urlopen(request, timeout=45) as response:
                payload = json.load(response)
            row = payload["PropertyTable"]["Properties"][0]
            smiles = row.get("SMILES") or row.get("ConnectivitySMILES")
            if not smiles:
                raise ValueError("PubChem response contained no SMILES")
            return {
                "status": "resolved", "query": identifier, "url": url,
                "cid": int(row["CID"]), "smiles": str(smiles),
                "connectivity_smiles": str(row.get("ConnectivitySMILES", smiles)),
            }
        except (HTTPError, URLError, TimeoutError, ValueError, KeyError) as exc:
            error = f"{type(exc).__name__}: {exc}"
            time.sleep(0.5 * (2**attempt))
    return {"status": "unresolved", "query": identifier, "url": url, "error": error}


def workbook_compounds(path: str | Path) -> list[dict]:
    # Deliberately do not load rows 4 onward: they contain the sealed outcomes.
    header = pd.read_excel(path, sheet_name=0, header=None, nrows=3)
    records = []
    for column in range(1, header.shape[1]):
        name = header.iat[1, column]
        cas = header.iat[2, column]
        if pd.isna(name) and pd.isna(cas):
            continue
        records.append({
            "workbook_column": int(column),
            "compound_name": "" if pd.isna(name) else str(name).strip(),
            "cas": "" if pd.isna(cas) else str(cas).strip(),
        })
    return records


def run(workbook: str | Path, output: str | Path, workers: int = 4) -> dict:
    workbook = Path(workbook).resolve()
    output = Path(output).resolve()
    compounds = workbook_compounds(workbook)
    resolved: dict[int, dict] = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        future_map = {}
        for record in compounds:
            query = record["cas"] or record["compound_name"]
            future_map[pool.submit(_lookup, query)] = record["workbook_column"]
        for completed, future in enumerate(as_completed(future_map), start=1):
            resolved[future_map[future]] = future.result()
            if completed % 25 == 0 or completed == len(future_map):
                print(f"PubChem identities: {completed}/{len(future_map)}", flush=True)

    entries = []
    for record in compounds:
        lookup = resolved[record["workbook_column"]]
        # CAS lookup is primary.  A failed CAS is not silently replaced by a
        # potentially ambiguous compound-name match.
        entries.append({**record, **lookup})
    payload = {
        "schema": "MetaSieve.AnastassiadisIdentitySidecar.v1",
        "outcome_rows_read": 0,
        "source_workbook": str(workbook),
        "source_sha256": sha256_file(workbook),
        "pubchem_endpoint": PUBCHEM,
        "n_compounds": len(entries),
        "n_resolved": sum(item["status"] == "resolved" for item in entries),
        "entries": entries,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: payload[key] for key in ("n_compounds", "n_resolved")}, indent=2))
    return payload


def parse_args(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--workbook", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--workers", type=int, default=4)
    return parser.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    run(args.workbook, args.output, args.workers)

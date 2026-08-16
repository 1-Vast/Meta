"""Label-blind structural coverage gate for A2S-CFES.

This gate reads source split metadata but no affinity values. It queries public
AlphaFold DB and PDBe metadata for source ``fit`` accessions and records whether
the conformational-state branch has a deployable structural substrate. Passing
this census authorizes only the external structure semantic gate, not affinity
training or a support-to-state operator.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
import time
from typing import Callable
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

import pandas as pd

from research.a2s.a2s_information_gate import (
    canonical,
    load_metadata,
    sha256_file,
    verify_lock,
)
from research.a2s.a2s_trace_stratum import DEFAULT_LOCK


ROOT = Path(__file__).resolve().parents[2]
TARGET_UNIPROT = ROOT / "dataset" / "public" / "chembl_37" / "processed" / "target_uniprot.json"
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "a2s_cfes_coverage_gate_2026-08-02.json"
DEFAULT_RECORDS = (
    ROOT / "reports" / "active" / "a2s_cfes_coverage_gate_records_2026-08-02.parquet"
)

USER_AGENT = "A2S-CFES/1.0 (public structural coverage audit)"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"
PDBE_API = "https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{accession}"
MAX_WORKERS = 12
REQUEST_ATTEMPTS = 3
REQUEST_TIMEOUT_SECONDS = 25

MIN_MAPPING_FRACTION = 0.95
MIN_ALPHAFOLD_FRACTION = 0.90
MIN_PDB_ANY_FRACTION = 0.70
MIN_PDB_TWO_FRACTION = 0.60
MIN_COMPONENTS_WITH_TWO_PDB = 80


@dataclass(frozen=True)
class FetchResult:
    status: str
    payload: object | None
    error: str | None


def fetch_json(url: str) -> FetchResult:
    """Fetch public metadata with bounded retries and explicit status."""

    for attempt in range(REQUEST_ATTEMPTS):
        request = Request(url, headers={"User-Agent": USER_AGENT})
        try:
            with urlopen(request, timeout=REQUEST_TIMEOUT_SECONDS) as response:
                return FetchResult("ok", json.load(response), None)
        except HTTPError as exc:
            if exc.code == 404:
                return FetchResult("not_found", None, "HTTP 404")
            error = f"HTTP {exc.code}"
        except (URLError, TimeoutError, OSError, json.JSONDecodeError) as exc:
            error = f"{type(exc).__name__}: {exc}"
        if attempt + 1 < REQUEST_ATTEMPTS:
            time.sleep(0.5 * (attempt + 1))
    return FetchResult("request_failed", None, error)


def parse_alphafold(result: FetchResult) -> dict[str, object]:
    if result.status != "ok" or not isinstance(result.payload, list) or not result.payload:
        return {
            "alphafold_status": result.status,
            "alphafold_version": None,
            "alphafold_model_date": None,
            "alphafold_cif_url": None,
            "alphafold_error": result.error,
        }
    entry = result.payload[0]
    if not isinstance(entry, dict):
        return {
            "alphafold_status": "invalid_response",
            "alphafold_version": None,
            "alphafold_model_date": None,
            "alphafold_cif_url": None,
            "alphafold_error": "first response entry is not an object",
        }
    return {
        "alphafold_status": "ok",
        "alphafold_version": entry.get("latestVersion"),
        "alphafold_model_date": entry.get("modelCreatedDate"),
        "alphafold_cif_url": entry.get("cifUrl"),
        "alphafold_error": None,
    }


def parse_pdbe(result: FetchResult, accession: str) -> dict[str, object]:
    entries: list[object] = []
    if result.status == "ok" and isinstance(result.payload, dict):
        candidate = result.payload.get(accession)
        if candidate is None:
            candidate = result.payload.get(accession.lower())
        if isinstance(candidate, list):
            entries = candidate
    pdb_ids = sorted(
        {
            str(entry["pdb_id"]).lower()
            for entry in entries
            if isinstance(entry, dict) and entry.get("pdb_id")
        }
    )
    status = result.status
    if result.status == "ok" and not pdb_ids:
        status = "no_structures"
    return {
        "pdbe_status": status,
        "pdb_structures": len(pdb_ids),
        "pdb_ids": "|".join(pdb_ids),
        "pdbe_error": result.error,
    }


def query_accession(
    target: str,
    component: str,
    accession: str | None,
    fetcher: Callable[[str], FetchResult] = fetch_json,
) -> dict[str, object]:
    row: dict[str, object] = {
        "target": target,
        "component": component,
        "accession": accession,
        "mapped": accession is not None,
    }
    if accession is None:
        row.update(parse_alphafold(FetchResult("unmapped", None, None)))
        row.update(parse_pdbe(FetchResult("unmapped", None, None), ""))
        return row
    row.update(parse_alphafold(fetcher(ALPHAFOLD_API.format(accession=accession))))
    row.update(parse_pdbe(fetcher(PDBE_API.format(accession=accession)), accession))
    return row


def summarise(records: pd.DataFrame) -> tuple[dict[str, object], dict[str, bool]]:
    targets = int(len(records))
    if targets == 0:
        raise ValueError("coverage gate has no fit targets")
    mapped = records.mapped.astype(bool)
    alphafold = records.alphafold_status.eq("ok")
    pdb_any = records.pdb_structures.ge(1)
    pdb_two = records.pdb_structures.ge(2)
    components_two = int(records.loc[pdb_two, "component"].nunique())
    summary: dict[str, object] = {
        "fit_targets": targets,
        "fit_components": int(records.component.nunique()),
        "mapped_targets": int(mapped.sum()),
        "unique_accessions": int(records.loc[mapped, "accession"].nunique()),
        "alphafold_targets": int(alphafold.sum()),
        "pdb_targets": int(pdb_any.sum()),
        "pdb_two_state_targets": int(pdb_two.sum()),
        "components_with_two_pdb_states": components_two,
        "mapping_fraction": float(mapped.mean()),
        "alphafold_fraction": float(alphafold.mean()),
        "pdb_any_fraction": float(pdb_any.mean()),
        "pdb_two_fraction": float(pdb_two.mean()),
        "median_pdb_structures": float(records.pdb_structures.median()),
        "total_target_pdb_structures": int(records.pdb_structures.sum()),
        "alphafold_status": records.alphafold_status.value_counts().sort_index().to_dict(),
        "pdbe_status": records.pdbe_status.value_counts().sort_index().to_dict(),
    }
    checks = {
        "mapping_fraction": summary["mapping_fraction"] >= MIN_MAPPING_FRACTION,
        "alphafold_fallback_fraction": (
            summary["alphafold_fraction"] >= MIN_ALPHAFOLD_FRACTION
        ),
        "experimental_pdb_fraction": summary["pdb_any_fraction"] >= MIN_PDB_ANY_FRACTION,
        "two_pdb_state_fraction": summary["pdb_two_fraction"] >= MIN_PDB_TWO_FRACTION,
        "components_with_two_pdb_states": components_two >= MIN_COMPONENTS_WITH_TWO_PDB,
    }
    return summary, checks


def run(output: Path, records_path: Path) -> dict[str, object]:
    lock = verify_lock(DEFAULT_LOCK)
    metadata = load_metadata(lock)
    if "affinity" in metadata.columns:
        raise AssertionError("coverage gate metadata unexpectedly contains affinity")
    fit = metadata.loc[metadata.role == "fit", ["target", "component"]].drop_duplicates()
    if fit.target.duplicated().any():
        raise RuntimeError("a fit target maps to more than one component")
    target_uniprot = json.loads(TARGET_UNIPROT.read_text(encoding="utf-8"))
    jobs = [
        (str(row.target), str(row.component), target_uniprot.get(str(row.target)))
        for row in fit.sort_values("target").itertuples(index=False)
    ]
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
        records = pd.DataFrame.from_records(
            executor.map(lambda item: query_accession(*item), jobs)
        )
    summary, checks = summarise(records)
    passed = all(checks.values())

    output.parent.mkdir(parents=True, exist_ok=True)
    records_path.parent.mkdir(parents=True, exist_ok=True)
    records.to_parquet(records_path, index=False)
    payload: dict[str, object] = {
        "schema": "a2s-cfes-coverage-gate-v1",
        "status": "LABEL_BLIND_SOURCE_FIT_METADATA_ONLY",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "source_roles_used": ["fit"],
            "source_affinity_read": False,
            "source_probe_labels_read": False,
            "source_locked_labels_read": False,
            "recipient_labels_read": False,
            "external_affinity_read": False,
            "alphafold_api": ALPHAFOLD_API,
            "pdbe_api": PDBE_API,
            "request_attempts": REQUEST_ATTEMPTS,
            "request_timeout_seconds": REQUEST_TIMEOUT_SECONDS,
        },
        "thresholds": {
            "mapping_fraction": MIN_MAPPING_FRACTION,
            "alphafold_fraction": MIN_ALPHAFOLD_FRACTION,
            "pdb_any_fraction": MIN_PDB_ANY_FRACTION,
            "pdb_two_fraction": MIN_PDB_TWO_FRACTION,
            "components_with_two_pdb_states": MIN_COMPONENTS_WITH_TWO_PDB,
        },
        "summary": summary,
        "decision": {
            "checks": checks,
            "pass": passed,
            "verdict": (
                "CFES_C0A_COVERAGE_PASS_PROCEED_EXTERNAL_SEMANTIC_GATE"
                if passed
                else "CFES_C0A_COVERAGE_FAIL_STOP"
            ),
            "authorizes": (
                "external structure-only semantic gate C0B"
                if passed
                else "no structural model or affinity training"
            ),
        },
        "artifacts": {
            "source_lock_sha256": sha256_file(DEFAULT_LOCK),
            "target_uniprot_sha256": sha256_file(TARGET_UNIPROT),
            "records": str(records_path.relative_to(ROOT)),
            "records_sha256": sha256_file(records_path),
        },
    }
    payload["content_sha256"] = sha256(canonical(payload).encode()).hexdigest()
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    args = parser.parse_args()
    payload = run(args.output, args.records)
    print(payload["decision"]["verdict"])
    print(canonical(payload["decision"]["checks"]))


if __name__ == "__main__":
    main()


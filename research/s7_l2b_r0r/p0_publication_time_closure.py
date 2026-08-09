"""Phase 0 item 4 — publication/time closure from official RCSB records.

For every development and confirmation PDB entry, retrieve the authoritative
initial release date and the primary citation DOI / PubMed ID from the RCSB Data
API, then build the closure:

    document key = canonical DOI, else canonical PubMed ID, else QUARANTINE.

A PDB ID is never substituted for missing publication metadata. The temporal
cutoff is frozen here, before any confirmation label or score is inspected.

Raw responses, the normalized table and the missingness report are all hashed.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import sys
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(r"D:\MetaSieve")
CORPUS = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "r0r1_raw_corpus"
CACHE = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "rcsb"
OUT = ROOT / "report" / "s7_l2b_r0r"
GRAPHQL = "https://data.rcsb.org/graphql"
BATCH = 50
TIME_CUTOFF = "2019-01-01"   # frozen here, before any confirmation score

QUERY = """query($ids:[String!]!){
  entries(entry_ids:$ids){
    rcsb_id
    rcsb_accession_info { initial_release_date }
    citation { id pdbx_database_id_DOI pdbx_database_id_PubMed }
  }
}"""


def sha_bytes(b):
    return hashlib.sha256(b).hexdigest()


def load_ids():
    ids = {}
    for name, cohort in (("monn_development_edge_corpus.jsonl.gz", "development"),
                         ("monn_additional_pdb_edge_corpus.jsonl.gz", "additional_pdb")):
        with gzip.open(CORPUS / name, "rt", encoding="utf-8") as f:
            for line in f:
                r = json.loads(line)
                ids.setdefault(r["pdb_id"].upper(), set()).add(cohort)
    return ids


def fetch(ids):
    body = json.dumps({"query": QUERY, "variables": {"ids": ids}}).encode()
    req = urllib.request.Request(GRAPHQL, data=body,
                                 headers={"Content-Type": "application/json"})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except (urllib.error.URLError, TimeoutError) as exc:
            if attempt == 4:
                raise
            time.sleep(2 * (attempt + 1))
    raise RuntimeError("unreachable")


def main():
    CACHE.mkdir(parents=True, exist_ok=True)
    ids_map = load_ids()
    all_ids = sorted(ids_map)
    print(f"unique PDB entries to resolve: {len(all_ids)}", flush=True)

    raw_path = CACHE / "rcsb_raw_responses.jsonl.gz"
    records = {}
    raw_hashes = []
    t0 = time.time()
    with gzip.open(raw_path, "wt", encoding="utf-8") as raw:
        for i in range(0, len(all_ids), BATCH):
            chunk = all_ids[i:i + BATCH]
            blob = fetch(chunk)
            raw_hashes.append(sha_bytes(blob))
            raw.write(json.dumps({"ids": chunk, "sha256": raw_hashes[-1],
                                  "body": json.loads(blob.decode())}) + "\n")
            data = json.loads(blob.decode()).get("data") or {}
            for e in (data.get("entries") or []):
                if not e:
                    continue
                pid = (e.get("rcsb_id") or "").upper()
                rel = ((e.get("rcsb_accession_info") or {})
                       .get("initial_release_date"))
                doi = pmid = None
                for c in (e.get("citation") or []):
                    if (c.get("id") or "").lower() == "primary":
                        doi = c.get("pdbx_database_id_DOI")
                        pmid = c.get("pdbx_database_id_PubMed")
                        break
                records[pid] = {"release_date": rel, "doi": doi, "pubmed": pmid}
            if (i // BATCH) % 20 == 0:
                print(f"  {min(i+BATCH, len(all_ids))}/{len(all_ids)} "
                      f"({time.time()-t0:.0f}s)", flush=True)

    miss = Counter()
    doc_key = {}
    for pid in all_ids:
        r = records.get(pid)
        if r is None:
            miss["entry_not_returned"] += 1
            continue
        if not r.get("release_date"):
            miss["missing_release_date"] += 1
        doi = (r.get("doi") or "").strip().lower() or None
        pm = str(r.get("pubmed") or "").strip() or None
        if doi:
            doc_key[pid] = f"doi:{doi}"
        elif pm and pm not in ("None", "0"):
            doc_key[pid] = f"pmid:{pm}"
        else:
            miss["no_primary_publication_identifier"] += 1

    dev_docs, add_docs = set(), set()
    dev_dates, add_dates = [], []
    for pid, cohorts in ids_map.items():
        k = doc_key.get(pid)
        d = (records.get(pid) or {}).get("release_date")
        if "development" in cohorts:
            if k:
                dev_docs.add(k)
            if d:
                dev_dates.append(d[:10])
        if "additional_pdb" in cohorts:
            if k:
                add_docs.add(k)
            if d:
                add_dates.append(d[:10])

    shared_docs = dev_docs & add_docs
    conf_after_cutoff = sum(1 for d in add_dates if d >= TIME_CUTOFF)

    norm_path = CACHE / "rcsb_normalized.json"
    norm_path.write_text(json.dumps({"records": records, "document_key": doc_key},
                                    sort_keys=True), encoding="utf-8")
    out = {
        "schema": "MetaSieve.S7L2B.P0.PublicationTimeClosureAudit.v1",
        "created_utc": "2026-08-09",
        "source": "RCSB Data API GraphQL, data.rcsb.org",
        "query_sha256": sha_bytes(QUERY.encode()),
        "batch_size": BATCH,
        "entries_requested": len(all_ids),
        "entries_returned": len(records),
        "raw_responses_path": str(raw_path.relative_to(ROOT)).replace("\\", "/"),
        "raw_responses_sha256": sha_bytes(raw_path.read_bytes()),
        "raw_batch_hashes_sha256": sha_bytes(json.dumps(raw_hashes).encode()),
        "normalized_path": str(norm_path.relative_to(ROOT)).replace("\\", "/"),
        "normalized_sha256": sha_bytes(norm_path.read_bytes()),
        "document_key_rule": "canonical DOI, else canonical PubMed ID, else QUARANTINE; "
                             "a PDB ID is NEVER substituted for missing publication metadata",
        "missingness": dict(miss),
        "entries_with_document_key": len(doc_key),
        "frozen_time_cutoff": TIME_CUTOFF,
        "cutoff_frozen_before_any_confirmation_score": True,
        "development_documents": len(dev_docs),
        "additional_pdb_documents": len(add_docs),
        "documents_shared_by_both_cohorts": len(shared_docs),
        "additional_pdb_entries_released_on_or_after_cutoff": conf_after_cutoff,
        "reading": "documents shared by both cohorts must be removed from any "
                   "confirmation set; entries without a primary publication "
                   "identifier are quarantined from formal confirmation",
    }
    (OUT / "PUBLICATION_TIME_CLOSURE_AUDIT.json").write_text(json.dumps(out, indent=2),
                                                             encoding="utf-8")
    print(json.dumps({k: out[k] for k in
                      ("entries_requested", "entries_returned", "missingness",
                       "entries_with_document_key", "development_documents",
                       "additional_pdb_documents", "documents_shared_by_both_cohorts",
                       "additional_pdb_entries_released_on_or_after_cutoff")}, indent=2))
    print(f"\nwrote {OUT / 'PUBLICATION_TIME_CLOSURE_AUDIT.json'}")


if __name__ == "__main__":
    sys.exit(main())

"""S7_L2B R0R-2d — new edge corpus manifest, including a typed-channel census.

Label-blind. Records what the new corpus actually contains so a future typed
channel Gate cannot be registered against channels the supervision does not
carry.
"""
from __future__ import annotations

import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
CORPUS = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "r0r1_raw_corpus"
OUT = ROOT / "report" / "s7_l2b_r0r"


def sha_file(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def survey(path):
    kinds = Counter()
    complexes_with_kind = Counter()
    n = 0
    res_counts, atom_counts, edge_counts, seq_lens = [], [], [], []
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            r = json.loads(line)
            n += 1
            ks = {e[2] for e in r["positive_typed_edges"]}
            for k in ks:
                complexes_with_kind[k] += 1
            for e in r["positive_typed_edges"]:
                kinds[e[2]] += 1
            atom_counts.append(len(r["atom_names"]))
            seq_lens.append(len(r["uniprot_sequence"]))
            edge_counts.append(len(r["positive_binary_edges"]))
            res_counts.append(len({e[0] for e in r["positive_binary_edges"]}))
    q = lambda a, p: float(np.percentile(a, p))  # noqa: E731
    return {
        "complexes": n,
        "typed_edges_by_channel": dict(kinds.most_common()),
        "complexes_containing_channel": dict(complexes_with_kind.most_common()),
        "channel_prevalence_fraction": {k: round(v / n, 4)
                                        for k, v in complexes_with_kind.most_common()},
        "ligand_heavy_atoms": {"p5": q(atom_counts, 5), "p50": q(atom_counts, 50),
                               "p95": q(atom_counts, 95), "max": int(max(atom_counts))},
        "sequence_length": {"p5": q(seq_lens, 5), "p50": q(seq_lens, 50),
                            "p95": q(seq_lens, 95), "max": int(max(seq_lens))},
        "positive_binary_edges_per_complex": {"p5": q(edge_counts, 5),
                                              "p50": q(edge_counts, 50),
                                              "p95": q(edge_counts, 95)},
        "distinct_positive_residues_per_complex": {"p5": q(res_counts, 5),
                                                   "p50": q(res_counts, 50),
                                                   "p95": q(res_counts, 95)},
        "complete_matrix_cells": int(sum(a * s for a, s in zip(atom_counts, seq_lens))),
        "positive_rate_over_complete_matrix": float(
            sum(edge_counts) / sum(a * s for a, s in zip(atom_counts, seq_lens))),
    }


def main():
    dev_p = CORPUS / "monn_development_edge_corpus.jsonl.gz"
    add_p = CORPUS / "monn_additional_pdb_edge_corpus.jsonl.gz"
    closure = json.loads((OUT / "NEW_CLOSURE_AND_SPLIT_MANIFEST.json").read_text())
    man = {
        "schema": "MetaSieve.S7L2B.R0R2d.NewEdgeCorpusManifest.v1",
        "created_utc": "2026-08-09",
        "identity_note": "this is a NEW corpus with NEW names. It does not inherit or "
                         "assert the historical 4067/701, 8646 or 524/157 identities, "
                         "nor any historical AP value.",
        "provenance": {
            "source_repository": "https://github.com/lishuya17/MONN",
            "pinned_commit": "f2b62ccf49c18a9502aa0eb0d582c6e0735ef200",
            "licence": "NON COMMERCIAL use only (MONN README)",
            "redistribution": "forbidden; sources live under gitignored dataset/raw/monn/",
            "affinity_tables_opened": 0,
        },
        "files": {
            "development": {"path": str(dev_p.relative_to(ROOT)).replace("\\", "/"),
                            "sha256": sha_file(dev_p), "bytes": dev_p.stat().st_size},
            "additional_pdb": {"path": str(add_p.relative_to(ROOT)).replace("\\", "/"),
                               "sha256": sha_file(add_p), "bytes": add_p.stat().st_size},
        },
        "ligand_identity": closure["ligand_identity"],
        "frozen_policies": closure["frozen_policies"],
        "development": survey(dev_p),
        "additional_pdb": survey(add_p),
    }
    dev = man["development"]["channel_prevalence_fraction"]
    man["typed_channel_evaluability"] = {
        "rule": "a channel is EVALUABLE only if it appears in at least 5% of development "
                "complexes; otherwise a per-complex AP for it is undefined in most units",
        "evaluable": [k for k, v in dev.items() if v >= 0.05],
        "not_evaluable": [k for k, v in dev.items() if v < 0.05],
        "consequence": "a typed-channel Gate may only be registered over the evaluable set",
    }
    (OUT / "NEW_EDGE_CORPUS_MANIFEST.json").write_text(json.dumps(man, indent=2),
                                                       encoding="utf-8")
    print(json.dumps({"development_channels": man["development"]["channel_prevalence_fraction"],
                      "additional_channels": man["additional_pdb"]["channel_prevalence_fraction"],
                      "typed_channel_evaluability": man["typed_channel_evaluability"],
                      "dev_positive_rate": man["development"]["positive_rate_over_complete_matrix"],
                      "dev_matrix_cells": man["development"]["complete_matrix_cells"]},
                     indent=2))
    print(f"\nwrote {OUT / 'NEW_EDGE_CORPUS_MANIFEST.json'}")


if __name__ == "__main__":
    main()

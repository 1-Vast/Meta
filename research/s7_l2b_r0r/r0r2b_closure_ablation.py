"""S7_L2B R0R-2b — closure relation ablation.

Diagnostic only. Measures which frozen closure relation creates the giant
component. This does NOT relax any closure; the frozen rule set is unchanged and
the ablation is reported alongside it so the cause of the topology is on record.
Label-blind.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import pickle
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(r"D:\MetaSieve")
CORPUS = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "r0r1_raw_corpus"
MONN = ROOT / "dataset" / "raw" / "monn" / "MONN" / "data"
CACHE = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "r0r2_closure"
OUT = ROOT / "report" / "s7_l2b_r0r"


def sha(b):
    return hashlib.sha256(b).hexdigest()


def load(p):
    with gzip.open(p, "rt", encoding="utf-8") as f:
        return [json.loads(l) for l in f]


class UF:
    def __init__(self, items):
        self.p = {i: i for i in items}

    def find(self, a):
        while self.p[a] != a:
            self.p[a] = self.p[self.p[a]]
            a = self.p[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def main():
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")

    dev = load(CORPUS / "monn_development_edge_corpus.jsonl.gz")
    add = load(CORPUS / "monn_additional_pdb_edge_corpus.jsonl.gz")
    for r in dev:
        r["cohort"] = "development"
    for r in add:
        r["cohort"] = "additional_pdb"
    rows = dev + add

    md = [pickle.load((MONN / "mol_dict").open("rb"), encoding="bytes"),
          pickle.load((MONN / "independent_dataset_mol_dict").open("rb"),
                      encoding="bytes")]
    ident = {}
    for ccd in {r["ligand_ccd"] for r in rows}:
        mol = None
        for d in md:
            m = d.get(ccd.encode("ascii", "ignore"))
            if m is not None:
                mol = m
                break
        if mol is None:
            continue
        w = Chem.Mol(mol)
        try:
            Chem.SanitizeMol(w)
            smi = Chem.MolToSmiles(w, isomericSmiles=False, canonical=True)
            scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=w)
        except Exception:
            continue
        ident[ccd] = (sha(smi.encode()), scaf)
    for r in rows:
        g = ident.get(r["ligand_ccd"])
        r["graph_key"], r["scaffold_smiles"] = g if g else (None, "")
        r["seq_key"] = sha(r["uniprot_sequence"].encode())
    rows = [r for r in rows if r["graph_key"]]

    hom = json.loads((CACHE / "homology_alignment_stats.json").read_text())
    hom_edges = [(a, b) for a, b, _i, _c, ok in hom if ok]

    def build(relations):
        ids = [r["source_key"] for r in rows]
        uf = UF(ids)
        if "exact_pdb" in relations:
            g = defaultdict(list)
            for r in rows:
                g[r["pdb_id"]].append(r["source_key"])
            for v in g.values():
                for m in v[1:]:
                    uf.union(v[0], m)
        if "exact_sequence" in relations:
            g = defaultdict(list)
            for r in rows:
                g[r["seq_key"]].append(r["source_key"])
            for v in g.values():
                for m in v[1:]:
                    uf.union(v[0], m)
        if "uniprot" in relations:
            g = defaultdict(list)
            for r in rows:
                g[r["uniprot_id"]].append(r["source_key"])
            for v in g.values():
                for m in v[1:]:
                    uf.union(v[0], m)
        if "exact_ligand_graph" in relations:
            g = defaultdict(list)
            for r in rows:
                g[r["graph_key"]].append(r["source_key"])
            for v in g.values():
                for m in v[1:]:
                    uf.union(v[0], m)
        if "murcko_scaffold" in relations:
            g = defaultdict(list)
            for r in rows:
                if r["scaffold_smiles"]:
                    g[r["scaffold_smiles"]].append(r["source_key"])
            for v in g.values():
                for m in v[1:]:
                    uf.union(v[0], m)
        if "homology_40pct" in relations:
            by = defaultdict(list)
            for r in rows:
                by[r["seq_key"]].append(r["source_key"])
            for a, b in hom_edges:
                if by.get(a) and by.get(b):
                    uf.union(by[a][0], by[b][0])
        comp = {i: uf.find(i) for i in ids}
        sizes = Counter(comp.values())
        dev_c = {comp[r["source_key"]] for r in rows if r["cohort"] == "development"}
        add_c = {comp[r["source_key"]] for r in rows if r["cohort"] == "additional_pdb"}
        return {
            "relations": sorted(relations),
            "components": len(sizes),
            "largest": max(sizes.values()),
            "largest_fraction": round(max(sizes.values()) / len(ids), 4),
            "components_ge5": sum(1 for v in sizes.values() if v >= 5),
            "shared_components_dev_vs_additional": len(dev_c & add_c),
        }

    PROTEIN = {"exact_pdb", "exact_sequence", "uniprot", "homology_40pct"}
    scenarios = {
        "protein_only": PROTEIN,
        "protein_plus_exact_ligand_graph": PROTEIN | {"exact_ligand_graph"},
        "protein_plus_scaffold": PROTEIN | {"murcko_scaffold"},
        "FROZEN_FULL_RULE_SET": PROTEIN | {"exact_ligand_graph", "murcko_scaffold"},
        "ligand_only_scaffold": {"murcko_scaffold"},
        "ligand_only_exact_graph": {"exact_ligand_graph"},
        "exact_sequence_only": {"exact_sequence"},
        "uniprot_only": {"uniprot"},
        "homology_only": {"homology_40pct"},
    }
    res = {k: build(v) for k, v in scenarios.items()}
    out = {"schema": "MetaSieve.S7L2B.R0R2b.ClosureAblation.v1",
           "created_utc": "2026-08-09",
           "purpose": "localize which frozen closure relation creates the giant component; "
                      "diagnostic only, the frozen rule set is unchanged",
           "total_complexes": len(rows),
           "scenarios": res}
    (OUT / "CLOSURE_RELATION_ABLATION.json").write_text(json.dumps(out, indent=2),
                                                        encoding="utf-8")
    print(f"{'scenario':38s} {'comps':>6s} {'largest':>8s} {'frac':>7s} {'>=5':>5s} {'shared':>7s}")
    for k, v in res.items():
        print(f"{k:38s} {v['components']:>6d} {v['largest']:>8d} "
              f"{v['largest_fraction']:>7.4f} {v['components_ge5']:>5d} "
              f"{v['shared_components_dev_vs_additional']:>7d}")
    print(f"\nwrote {OUT / 'CLOSURE_RELATION_ABLATION.json'}")


if __name__ == "__main__":
    main()

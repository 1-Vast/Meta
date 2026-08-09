"""S7_L2B R0R-2e — development split under a PROTEIN partition with ligand
disjointness enforced as a FILTER.

Methodological correction established by R0R-2: union-MERGING ligand identity
into the inference partition chains transitively through promiscuous cofactors
and produces a 93.19% giant component. Enforcing the same leakage control as a
DISJOINTNESS FILTER between train and held-out sets achieves protein AND ligand
closure without destroying the topology.

Inference units  = protein closure components (exact PDB, exact sequence,
                   UniProt, 40% alignment-verified homology).
Ligand closure   = held-out complexes sharing an exact ligand graph (and,
                   in the strict stratum, a Bemis-Murcko scaffold) with ANY
                   training complex are removed from evaluation.

Label-blind: only edge counts and identities are used. No model, no AP.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import pickle
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
CORPUS = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "r0r1_raw_corpus"
MONN = ROOT / "dataset" / "raw" / "monn" / "MONN" / "data"
CACHE = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "r0r2_closure"
OUT = ROOT / "report" / "s7_l2b_r0r"
SEED = 20260810
TEST_TARGET_FRACTION = 0.20


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


def build_rows():
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
            ident[ccd] = (sha(smi.encode()), MurckoScaffold.MurckoScaffoldSmiles(mol=w))
        except Exception:
            continue
    out = []
    for r in rows:
        g = ident.get(r["ligand_ccd"])
        if not g:
            continue
        r["graph_key"], r["scaffold"] = g
        r["seq_key"] = sha(r["uniprot_sequence"].encode())
        out.append(r)
    return out


def protein_components(rows):
    hom = json.loads((CACHE / "homology_alignment_stats.json").read_text())
    hom_edges = [(a, b) for a, b, _i, _c, ok in hom if ok]
    ids = [r["source_key"] for r in rows]
    uf = UF(ids)
    for keyfn in (lambda r: r["pdb_id"], lambda r: r["seq_key"], lambda r: r["uniprot_id"]):
        g = defaultdict(list)
        for r in rows:
            g[keyfn(r)].append(r["source_key"])
        for v in g.values():
            for m in v[1:]:
                uf.union(v[0], m)
    by = defaultdict(list)
    for r in rows:
        by[r["seq_key"]].append(r["source_key"])
    for a, b in hom_edges:
        if by.get(a) and by.get(b):
            uf.union(by[a][0], by[b][0])
    return {i: uf.find(i) for i in ids}


def main():
    rows = build_rows()
    comp = protein_components(rows)
    dev = [r for r in rows if r["cohort"] == "development"]
    add = [r for r in rows if r["cohort"] == "additional_pdb"]
    dev_comps = {comp[r["source_key"]] for r in dev}

    # ---- development split: whole protein components to train / held-out ----
    by_comp = defaultdict(list)
    for r in dev:
        by_comp[comp[r["source_key"]]].append(r)
    comps = sorted(by_comp, key=lambda c: (-len(by_comp[c]), c))
    rng = np.random.default_rng(SEED)
    order = list(comps)
    rng.shuffle(order)
    target = int(TEST_TARGET_FRACTION * len(dev))
    test_comps, n = set(), 0
    for c in order:
        if n >= target:
            break
        test_comps.add(c)
        n += len(by_comp[c])
    train = [r for r in dev if comp[r["source_key"]] not in test_comps]
    heldout_raw = [r for r in dev if comp[r["source_key"]] in test_comps]

    train_graphs = {r["graph_key"] for r in train}
    train_scaffolds = {r["scaffold"] for r in train if r["scaffold"]}

    ho_graph_disjoint = [r for r in heldout_raw if r["graph_key"] not in train_graphs]
    ho_strict = [r for r in ho_graph_disjoint
                 if not r["scaffold"] or r["scaffold"] not in train_scaffolds]

    def census(rs, label):
        cs = Counter(comp[r["source_key"]] for r in rs)
        return {
            "label": label,
            "complexes": len(rs),
            "components": len(cs),
            "largest_component_fraction": round(max(cs.values()) / len(rs), 4) if rs else 0,
            "components_ge5_complexes": sum(1 for v in cs.values() if v >= 5),
            "exact_sequences": len({r["seq_key"] for r in rs}),
            "exact_ligand_graphs": len({r["graph_key"] for r in rs}),
            "scaffolds": len({r["scaffold"] for r in rs if r["scaffold"]}),
            "positive_binary_edges": sum(len(r["positive_binary_edges"]) for r in rs),
        }

    # ---- confirmation cohort: additional-PDB, protein-disjoint from ALL development ----
    conf_protein = [r for r in add if comp[r["source_key"]] not in dev_comps]
    dev_graphs = {r["graph_key"] for r in dev}
    dev_scaffolds = {r["scaffold"] for r in dev if r["scaffold"]}
    conf_graph = [r for r in conf_protein if r["graph_key"] not in dev_graphs]
    conf_strict = [r for r in conf_graph
                   if not r["scaffold"] or r["scaffold"] not in dev_scaffolds]

    report = {
        "schema": "MetaSieve.S7L2B.R0R2e.SplitConstruction.v1",
        "created_utc": "2026-08-09",
        "methodological_correction": {
            "rejected": "union-MERGING ligand identity into the inference partition; "
                        "measured 93.19% giant component (R0R-2)",
            "adopted": "PROTEIN closure components as inference units; ligand-graph and "
                       "scaffold closure enforced as DISJOINTNESS FILTERS between train "
                       "and held-out, which achieves the same leakage control without "
                       "transitive chaining",
            "claim_scope": "generalisation over PROTEINS; the ligand-disjoint strata "
                           "additionally test generalisation to unseen ligands",
        },
        "seed": SEED,
        "test_target_fraction": TEST_TARGET_FRACTION,
        "development": {
            "train": census(train, "train"),
            "heldout_protein_disjoint": census(heldout_raw, "heldout_protein_disjoint"),
            "heldout_plus_ligand_graph_disjoint": census(ho_graph_disjoint, "heldout_graph_disjoint"),
            "heldout_plus_scaffold_disjoint_STRICT": census(ho_strict, "heldout_strict"),
        },
        "confirmation": {
            "protein_disjoint": census(conf_protein, "conf_protein_disjoint"),
            "plus_ligand_graph_disjoint": census(conf_graph, "conf_graph_disjoint"),
            "plus_scaffold_disjoint_STRICT": census(conf_strict, "conf_strict"),
        },
    }
    (OUT / "SPLIT_CONSTRUCTION_CENSUS.json").write_text(json.dumps(report, indent=2),
                                                        encoding="utf-8")
    for section in ("development", "confirmation"):
        print(f"\n=== {section} ===")
        for k, v in report[section].items():
            print(f"  {v['label']:34s} n={v['complexes']:6d} comp={v['components']:5d} "
                  f"largest={v['largest_component_fraction']:.4f} "
                  f">=5={v['components_ge5_complexes']:4d} "
                  f"graphs={v['exact_ligand_graphs']:5d} edges={v['positive_binary_edges']:7d}")
    print(f"\nwrote {OUT / 'SPLIT_CONSTRUCTION_CENSUS.json'}")


if __name__ == "__main__":
    main()

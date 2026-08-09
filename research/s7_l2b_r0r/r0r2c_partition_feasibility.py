"""S7_L2B R0R-2c — feasibility projection for an alternative inference partition.

Label-blind topology arithmetic ONLY. This does NOT adopt, relax or replace the
frozen closure. It answers one decision-relevant question with numbers so the
adjudication in R0R-2 can be made on evidence rather than preference:

    IF a future preregistration were to define the inference partition on the
    PROTEIN side alone (exact PDB, exact sequence, UniProt, 40% homology), with
    ligand leakage handled by held-out-ligand evaluation and the BX wrong-ligand
    control rather than by closure, what cohort would remain?

No model is fitted. No label is read. No threshold is changed.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import pickle
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
            ident[ccd] = (sha(Chem.MolToSmiles(w, isomericSmiles=False,
                                               canonical=True).encode()),
                          MurckoScaffold.MurckoScaffoldSmiles(mol=w))
        except Exception:
            continue
    for r in rows:
        g = ident.get(r["ligand_ccd"])
        r["graph_key"], r["scaffold_smiles"] = g if g else (None, "")
        r["seq_key"] = sha(r["uniprot_sequence"].encode())
    rows = [r for r in rows if r["graph_key"]]

    hom = json.loads((CACHE / "homology_alignment_stats.json").read_text())
    hom_edges = [(a, b) for a, b, _i, _c, ok in hom if ok]

    ids = [r["source_key"] for r in rows]
    uf = UF(ids)
    for keyfn in (lambda r: r["pdb_id"], lambda r: r["seq_key"],
                  lambda r: r["uniprot_id"]):
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

    comp = {i: uf.find(i) for i in ids}
    dev_comps = {comp[r["source_key"]] for r in rows if r["cohort"] == "development"}
    conf_rows = [r for r in rows
                 if r["cohort"] == "additional_pdb" and comp[r["source_key"]] not in dev_comps]
    dev_rows = [r for r in rows if r["cohort"] == "development"]

    def census(rs, label):
        cs = Counter(comp[r["source_key"]] for r in rs)
        edges = sum(len(r["positive_binary_edges"]) for r in rs)
        return {
            "label": label,
            "complexes": len(rs),
            "components": len(cs),
            "largest_component": max(cs.values()) if cs else 0,
            "largest_component_fraction": round(max(cs.values()) / len(rs), 4) if rs else 0,
            "components_ge5_complexes": sum(1 for v in cs.values() if v >= 5),
            "exact_sequences": len({r["seq_key"] for r in rs}),
            "exact_ligand_graphs": len({r["graph_key"] for r in rs}),
            "scaffolds": len({r["scaffold_smiles"] for r in rs if r["scaffold_smiles"]}),
            "positive_binary_edges": edges,
        }

    dev_c = census(dev_rows, "development_under_protein_only_partition")
    conf_c = census(conf_rows, "confirmation_candidate_protein_disjoint")

    MIN = {"complexes": 500, "components": 100, "exact_sequences": 250,
           "exact_ligand_graphs": 250, "positive_binary_edges": 2000,
           "largest_component_fraction_max": 0.25, "components_ge5_complexes": 60}
    checks = {
        "complexes": (conf_c["complexes"], MIN["complexes"], conf_c["complexes"] >= MIN["complexes"]),
        "components": (conf_c["components"], MIN["components"], conf_c["components"] >= MIN["components"]),
        "exact_sequences": (conf_c["exact_sequences"], MIN["exact_sequences"], conf_c["exact_sequences"] >= MIN["exact_sequences"]),
        "exact_ligand_graphs": (conf_c["exact_ligand_graphs"], MIN["exact_ligand_graphs"], conf_c["exact_ligand_graphs"] >= MIN["exact_ligand_graphs"]),
        "positive_binary_edges": (conf_c["positive_binary_edges"], MIN["positive_binary_edges"], conf_c["positive_binary_edges"] >= MIN["positive_binary_edges"]),
        "largest_component_fraction": (conf_c["largest_component_fraction"], MIN["largest_component_fraction_max"], conf_c["largest_component_fraction"] <= MIN["largest_component_fraction_max"]),
        "components_ge5_complexes": (conf_c["components_ge5_complexes"], MIN["components_ge5_complexes"], conf_c["components_ge5_complexes"] >= MIN["components_ge5_complexes"]),
    }
    out = {
        "schema": "MetaSieve.S7L2B.R0R2c.PartitionFeasibilityProjection.v1",
        "created_utc": "2026-08-09",
        "status": "PROJECTION_ONLY_NOT_AN_ADOPTED_CLOSURE",
        "hypothetical_partition": "exact_pdb + exact_sequence + uniprot + 40% homology; "
                                  "ligand leakage NOT closed, to be handled by held-out-ligand "
                                  "evaluation and the BX wrong-ligand control",
        "development": dev_c,
        "confirmation_candidate": conf_c,
        "confirmation_capacity_checks": {
            k: {"observed": v[0], "required": v[1], "pass": bool(v[2])}
            for k, v in checks.items()},
        "confirmation_all_checks_pass": all(v[2] for v in checks.values()),
        "caveat": "this projection does NOT establish admissibility. Publication/time closure "
                  "(R0R-3) is not applied here, and ligand-side leakage would remain uncontrolled "
                  "by closure. Adopting it would be a new registered estimand decision.",
    }
    (OUT / "PARTITION_FEASIBILITY_PROJECTION.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()

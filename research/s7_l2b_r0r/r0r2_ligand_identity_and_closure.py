"""S7_L2B R0R-2 — new ligand identities and closure topology census.

Label-blind. Reads only the R0R-1 edge corpora and the two MONN molecule
dictionaries. Never opens either MONN TSV.

Frozen policies (recorded in the output manifest):
  P1 RDKit version is pinned and recorded.
  P2 Deserialization: pickle.load(..., encoding="bytes"). encoding="latin1"
     fails on these files with 'Bad pickle format: bad endian ID'; plain
     pickle.load fails likewise. This is the only working policy and is frozen.
  P3 Sanitization: molecules arrive pre-sanitized from the MONN pipeline. We
     attempt Chem.SanitizeMol on a copy; failures are recorded and the complex
     is quarantined, never silently dropped.
  P4 Exact ligand graph identity = SHA-256 of the canonical SMILES computed with
     isomericSmiles=False, canonical=True. Stereoisomers therefore collapse to
     one identity, which is the CONSERVATIVE direction for leakage closure.
  P5 Scaffold = Bemis-Murcko scaffold SMILES. An EMPTY scaffold (acyclic ligand)
     is a sentinel and generates NO closure edge; treating "" as a shared
     scaffold would merge every acyclic ligand into one artificial blob.
  P6 Homology: MMseqs2 generates CANDIDATE pairs only. The authoritative
     identity/coverage numbers come from parasail Smith-Waterman, BLOSUM62,
     gap_open=10, gap_extend=1, identity = matches / local alignment length,
     coverage = alignment length / min(len_a, len_b). Threshold: identity >= 0.40
     AND coverage >= 0.80.

No historical count is targeted. The census is reported as measured.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import pickle
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
CORPUS = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "r0r1_raw_corpus"
MONN = ROOT / "dataset" / "raw" / "monn" / "MONN" / "data"
OUT = ROOT / "report" / "s7_l2b_r0r"
CACHE = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "r0r2_closure"
MMSEQS = ROOT / "tools" / "mmseqs2" / "mmseqs" / "bin" / "mmseqs.exe"

MIN_IDENTITY = 0.40
MIN_COVERAGE = 0.80
CANDIDATE_PREFILTER_IDENTITY = 0.25   # permissive; parasail is the authority


def sha256_bytes(b: bytes) -> str:
    return hashlib.sha256(b).hexdigest()


def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def load_corpus(path: Path):
    rows = []
    with gzip.open(path, "rt", encoding="utf-8") as f:
        for line in f:
            rows.append(json.loads(line))
    return rows


def load_mol_dict(path: Path):
    with path.open("rb") as f:
        return pickle.load(f, encoding="bytes")


def ligand_identities(mol_dicts, needed_ccds):
    """CCD -> (graph_key, graph_smiles, scaffold_smiles) with a failure census."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")

    out, fail = {}, Counter()
    for ccd in sorted(needed_ccds):
        mol = None
        for d in mol_dicts:
            m = d.get(ccd.encode("ascii", "ignore"))
            if m is not None:
                mol = m
                break
        if mol is None:
            fail["ccd_absent_from_mol_dict"] += 1
            continue
        try:
            work = Chem.Mol(mol)
            Chem.SanitizeMol(work)
        except Exception:
            fail["sanitize_failed"] += 1
            continue
        try:
            smi = Chem.MolToSmiles(work, isomericSmiles=False, canonical=True)
        except Exception:
            fail["canonical_smiles_failed"] += 1
            continue
        if not smi:
            fail["empty_canonical_smiles"] += 1
            continue
        try:
            scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=work)
        except Exception:
            scaf = ""
        out[ccd] = {"graph_key": sha256_bytes(smi.encode()), "graph_smiles": smi,
                    "scaffold_smiles": scaf, "n_atoms": work.GetNumAtoms()}
    return out, fail


def mmseqs_candidate_pairs(seqs: dict[str, str], workdir: Path):
    """Candidate homologous sequence pairs. MMseqs2 is a CANDIDATE GENERATOR."""
    workdir.mkdir(parents=True, exist_ok=True)
    fasta = workdir / "seqs.fasta"
    with fasta.open("w") as f:
        for sid, s in sorted(seqs.items()):
            f.write(f">{sid}\n{s}\n")
    res = workdir / "hits.m8"
    tmp = workdir / "tmp"
    cmd = [str(MMSEQS), "easy-search", str(fasta), str(fasta), str(res), str(tmp),
           "--min-seq-id", str(CANDIDATE_PREFILTER_IDENTITY), "-c", "0.5",
           "--cov-mode", "0", "-s", "7.5", "--max-seqs", "4000", "-v", "1"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"mmseqs easy-search failed rc={r.returncode} "
                           f"stdout={r.stdout[-500:]} stderr={r.stderr[-500:]}")
    pairs = set()
    with res.open() as f:
        for line in f:
            p = line.split("\t")
            if len(p) < 2:
                continue
            a, b = p[0], p[1]
            if a != b:
                pairs.add((a, b) if a < b else (b, a))
    return sorted(pairs), sha256_file(res)


def parasail_verify(pairs, seqs):
    """Authoritative identity/coverage. Returns accepted edges + full stats."""
    import parasail
    matrix = parasail.blosum62
    accepted, stats = [], []
    for a, b in pairs:
        sa, sb = seqs[a], seqs[b]
        res = parasail.sw_trace_striped_16(sa, sb, 10, 1, matrix)
        tr = res.traceback
        aln_len = len(tr.comp)
        if aln_len == 0:
            continue
        matches = sum(1 for ch in tr.comp if ch == "|")
        identity = matches / aln_len
        coverage = aln_len / min(len(sa), len(sb))
        ok = identity >= MIN_IDENTITY and coverage >= MIN_COVERAGE
        stats.append((a, b, round(identity, 5), round(coverage, 5), bool(ok)))
        if ok:
            accepted.append((a, b))
    return accepted, stats


class UnionFind:
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-homology", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    import rdkit

    dev = load_corpus(CORPUS / "monn_development_edge_corpus.jsonl.gz")
    add = load_corpus(CORPUS / "monn_additional_pdb_edge_corpus.jsonl.gz")
    for r in dev:
        r["cohort"] = "development"
    for r in add:
        r["cohort"] = "additional_pdb"
    rows = dev + add
    print(f"loaded {len(dev)} development + {len(add)} additional_pdb records",
          flush=True)

    md_dev = load_mol_dict(MONN / "mol_dict")
    md_add = load_mol_dict(MONN / "independent_dataset_mol_dict")
    needed = {r["ligand_ccd"] for r in rows}
    lig, ligfail = ligand_identities([md_dev, md_add], needed)
    print(f"ligand identities: {len(lig)}/{len(needed)} resolved; failures={dict(ligfail)}",
          flush=True)

    kept, quarantined = [], []
    for r in rows:
        info = lig.get(r["ligand_ccd"])
        if info is None:
            quarantined.append({"source_key": r["source_key"],
                                "reason": "ligand_identity_unresolved"})
            continue
        r["graph_key"] = info["graph_key"]
        r["scaffold_smiles"] = info["scaffold_smiles"]
        r["seq_key"] = sha256_bytes(r["uniprot_sequence"].encode())
        kept.append(r)
    print(f"kept {len(kept)}; quarantined {len(quarantined)}", flush=True)

    seqs = {}
    for r in kept:
        seqs.setdefault(r["seq_key"], r["uniprot_sequence"])
    print(f"unique exact sequences: {len(seqs)}", flush=True)

    hom_edges, hom_stats, hits_hash = [], [], None
    if not args.skip_homology:
        print("running MMseqs2 candidate search ...", flush=True)
        cand, hits_hash = mmseqs_candidate_pairs(seqs, CACHE / "mmseqs")
        print(f"candidate sequence pairs: {len(cand)}; verifying with parasail ...",
              flush=True)
        hom_edges, hom_stats = parasail_verify(cand, seqs)
        print(f"alignment-verified homology edges: {len(hom_edges)}", flush=True)

    # ---------------- union closure ----------------
    ids = [r["source_key"] for r in kept]
    uf = UnionFind(ids)
    rel_counts = Counter()

    def link_by(keyfn, name):
        groups = defaultdict(list)
        for r in kept:
            k = keyfn(r)
            if k is None:
                continue
            groups[k].append(r["source_key"])
        for k, members in groups.items():
            for m in members[1:]:
                uf.union(members[0], m)
            if len(members) > 1:
                rel_counts[name] += len(members) - 1

    link_by(lambda r: r["pdb_id"], "exact_pdb")
    link_by(lambda r: r["seq_key"], "exact_sequence")
    link_by(lambda r: r["uniprot_id"], "uniprot")
    link_by(lambda r: r["graph_key"], "exact_ligand_graph")
    link_by(lambda r: r["scaffold_smiles"] or None, "murcko_scaffold")

    by_seq = defaultdict(list)
    for r in kept:
        by_seq[r["seq_key"]].append(r["source_key"])
    for a, b in hom_edges:
        if by_seq.get(a) and by_seq.get(b):
            uf.union(by_seq[a][0], by_seq[b][0])
            rel_counts["homology_40pct"] += 1

    comp_of = {i: uf.find(i) for i in ids}
    sizes = Counter(comp_of.values())
    n_comp = len(sizes)
    largest = max(sizes.values())

    def census(subset_rows):
        sk = [r["source_key"] for r in subset_rows]
        cs = Counter(comp_of[k] for k in sk)
        return {
            "complexes": len(sk),
            "exact_sequences": len({r["seq_key"] for r in subset_rows}),
            "exact_ligand_graphs": len({r["graph_key"] for r in subset_rows}),
            "scaffolds": len({r["scaffold_smiles"] for r in subset_rows if r["scaffold_smiles"]}),
            "acyclic_ligands_no_scaffold": sum(1 for r in subset_rows if not r["scaffold_smiles"]),
            "components_touched": len(cs),
            "largest_component_complexes": max(cs.values()) if cs else 0,
            "largest_component_fraction": (max(cs.values()) / len(sk)) if sk else 0.0,
            "positive_binary_edges": sum(len(r["positive_binary_edges"]) for r in subset_rows),
            "positive_typed_edges": sum(len(r["positive_typed_edges"]) for r in subset_rows),
        }

    dev_rows = [r for r in kept if r["cohort"] == "development"]
    add_rows = [r for r in kept if r["cohort"] == "additional_pdb"]
    dev_comps = {comp_of[r["source_key"]] for r in dev_rows}
    add_comps = {comp_of[r["source_key"]] for r in add_rows}

    report = {
        "schema": "MetaSieve.S7L2B.R0R2.ClosureAndSplitManifest.v1",
        "created_utc": "2026-08-09",
        "step": "R0R-2",
        "frozen_policies": {
            "rdkit_version": rdkit.__version__,
            "python_version": sys.version.split()[0],
            "pickle_encoding": "bytes",
            "pickle_encoding_alternatives_tested": {
                "latin1": "FAILED: RuntimeError Bad pickle format: bad endian ID",
                "default": "FAILED"},
            "sanitization": "Chem.SanitizeMol on a copy; failures quarantined, never dropped silently",
            "exact_graph_identity": "sha256(MolToSmiles(isomericSmiles=False, canonical=True))",
            "stereo_policy": "stereoisomers collapse to one identity (conservative for leakage)",
            "scaffold": "Bemis-Murcko scaffold SMILES; EMPTY scaffold generates NO closure edge",
            "homology_candidate_generator": "MMseqs2 easy-search, prefilter min-seq-id "
                                            f"{CANDIDATE_PREFILTER_IDENTITY}, -c 0.5",
            "homology_authority": "parasail sw_trace_striped_16, BLOSUM62, gap_open=10, gap_extend=1",
            "identity_definition": "matches / local alignment length",
            "coverage_definition": "alignment length / min(len_a, len_b)",
            "thresholds": {"identity": MIN_IDENTITY, "coverage": MIN_COVERAGE},
        },
        "inputs": {
            "development_corpus_sha256": sha256_file(
                CORPUS / "monn_development_edge_corpus.jsonl.gz"),
            "additional_pdb_corpus_sha256": sha256_file(
                CORPUS / "monn_additional_pdb_edge_corpus.jsonl.gz"),
            "mol_dict_sha256": sha256_file(MONN / "mol_dict"),
            "independent_mol_dict_sha256": sha256_file(MONN / "independent_dataset_mol_dict"),
            "mmseqs_hits_sha256": hits_hash,
            "affinity_tables_opened": 0,
        },
        "ligand_identity": {
            "ccds_needed": len(needed),
            "ccds_resolved": len(lig),
            "failures": dict(ligfail),
            "quarantined_complexes": len(quarantined),
        },
        "homology": {
            "unique_sequences": len(seqs),
            "mmseqs_candidate_pairs": len(hom_stats),
            "alignment_verified_edges": len(hom_edges),
            "rejected_by_parasail": len(hom_stats) - len(hom_edges),
        },
        "closure_relations_used": dict(rel_counts),
        "topology": {
            "total_complexes": len(kept),
            "union_components": n_comp,
            "largest_component_complexes": largest,
            "largest_component_fraction": largest / len(kept),
            "component_size_p50": float(np.median(list(sizes.values()))),
            "components_with_at_least_5_complexes": sum(1 for v in sizes.values() if v >= 5),
        },
        "cohort_census": {
            "development": census(dev_rows),
            "additional_pdb": census(add_rows),
        },
        "cohort_separation": {
            "components_shared_by_both_cohorts": len(dev_comps & add_comps),
            "development_only_components": len(dev_comps - add_comps),
            "additional_pdb_only_components": len(add_comps - dev_comps),
            "note": "a shared component means the additional-PDB cohort is NOT independent "
                    "of development under the frozen closure",
        },
    }

    (CACHE / "component_assignments.json").write_text(
        json.dumps({k: comp_of[k] for k in sorted(comp_of)}, sort_keys=True),
        encoding="utf-8")
    (CACHE / "homology_alignment_stats.json").write_text(
        json.dumps(hom_stats), encoding="utf-8")
    (CACHE / "quarantined.json").write_text(json.dumps(quarantined, indent=2),
                                            encoding="utf-8")
    report["outputs"] = {
        "component_assignments_sha256": sha256_file(CACHE / "component_assignments.json"),
        "homology_alignment_stats_sha256": sha256_file(CACHE / "homology_alignment_stats.json"),
    }
    (OUT / "NEW_CLOSURE_AND_SPLIT_MANIFEST.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({k: report[k] for k in
                      ("topology", "cohort_census", "cohort_separation",
                       "closure_relations_used", "homology")}, indent=2))
    print(f"\nwrote {OUT / 'NEW_CLOSURE_AND_SPLIT_MANIFEST.json'}")


if __name__ == "__main__":
    main()

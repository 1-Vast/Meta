"""S7_L2B — corpus assembly, atom/residue contract validation, features, split.

Registered by research/s7_l2b_r0r/PREREG_S7_L2B_UNIFIED.md
(SHA-256 2c333f223ae450c566cc62b1a3b276ff59c065c38348005ad9504ac1930b9a92,
committed as ce186f4 BEFORE any model existed).

ATOM CONTRACT (validated on every record, fail-closed):
  MONN `atom_name` lists the deposited ligand atoms INCLUDING hydrogens, while
  `mol_dict[ccd]` is a heavy-atom RDKit molecule. `atom_idx` is a plain identity
  range and is NOT a mapping into the molecule. The mapping is therefore:

      heavy_positions = [k for k, name in enumerate(atom_names)
                         if element_is_hydrogen(name)]        # complement
      atom_slot s  ->  molecule atom index = rank of s among heavy_positions

  Verified empirically: heavy-atom count equals mol.GetNumAtoms() and NO
  positive edge ever lands on a hydrogen. Both are asserted per record here.

No affinity table is opened. Label-blind with respect to affinity.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import pickle
import re
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
CORPUS = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "r0r1_raw_corpus"
MONN = ROOT / "dataset" / "raw" / "monn" / "MONN" / "data"
CACHE = ROOT / "dataset" / "processed" / "s7_l2b_r0r"
OUT = ROOT / "report" / "s7_l2b_r0r"

SEED_SPLIT = 20260810
TEST_TARGET_FRACTION = 0.20
AA = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {c: i for i, c in enumerate(AA)}
RADII = (3, 7, 15, 31)
RES_DIM = 20 + 20 * len(RADII) + 4          # 104
ELEMENTS = ["C", "N", "O", "S", "P", "F", "CL", "BR", "I", "B", "SE"]
ATOM_DIM = len(ELEMENTS) + 1 + 6 + 5 + 5 + 6 + 2 + 5   # 41


def sha(b) -> str:
    return hashlib.sha256(b).hexdigest()


def is_hydrogen_name(name: str) -> bool:
    """PDB atom-name convention: strip leading digits, then a leading H/D is
    hydrogen/deuterium. Validated by an exact heavy-count match against the
    molecule on every record."""
    s = re.sub(r"^[0-9]+", "", name.strip().upper())
    return bool(s) and s[0] in ("H", "D")


# ----------------------------------------------------------------- features
def residue_features(seq: str) -> np.ndarray:
    """Explicit non-PLM local sequence features. Deterministic."""
    L = len(seq)
    oh = np.zeros((L, 20), dtype=np.float32)
    for i, ch in enumerate(seq):
        j = AA_IDX.get(ch)
        if j is not None:
            oh[i, j] = 1.0
    cum = np.zeros((L + 1, 20), dtype=np.float32)
    np.cumsum(oh, axis=0, out=cum[1:])
    feats = [oh]
    for r in RADII:
        lo = np.clip(np.arange(L) - r, 0, L)
        hi = np.clip(np.arange(L) + r + 1, 0, L)
        win = cum[hi] - cum[lo]
        win /= np.maximum((hi - lo)[:, None], 1).astype(np.float32)
        feats.append(win)
    pos = np.arange(L, dtype=np.float32)
    extra = np.stack([
        pos / max(L - 1, 1),
        np.log1p(pos) / 10.0,
        np.log1p(L - 1 - pos) / 10.0,
        np.full(L, np.log1p(L) / 10.0, dtype=np.float32),
    ], axis=1)
    feats.append(extra)
    return np.concatenate(feats, axis=1).astype(np.float32)


def atom_features(mol) -> np.ndarray:
    from rdkit import Chem
    hyb = [Chem.rdchem.HybridizationType.SP, Chem.rdchem.HybridizationType.SP2,
           Chem.rdchem.HybridizationType.SP3, Chem.rdchem.HybridizationType.SP3D,
           Chem.rdchem.HybridizationType.SP3D2]
    n = mol.GetNumAtoms()
    X = np.zeros((n, ATOM_DIM), dtype=np.float32)
    ri = mol.GetRingInfo()
    for a in mol.GetAtoms():
        i, o = a.GetIdx(), 0
        sym = a.GetSymbol().upper()
        if sym in ELEMENTS:
            X[i, ELEMENTS.index(sym)] = 1.0
        else:
            X[i, len(ELEMENTS)] = 1.0
        o = len(ELEMENTS) + 1
        X[i, o + min(a.GetDegree(), 5)] = 1.0
        o += 6
        X[i, o + int(np.clip(a.GetFormalCharge(), -2, 2)) + 2] = 1.0
        o += 5
        X[i, o + min(a.GetTotalNumHs(), 4)] = 1.0
        o += 5
        h = a.GetHybridization()
        X[i, o + (hyb.index(h) if h in hyb else 5)] = 1.0
        o += 6
        X[i, o] = float(a.GetIsAromatic())
        X[i, o + 1] = float(a.IsInRing())
        o += 2
        for k, size in enumerate((3, 4, 5, 6, 7)):
            if ri.IsAtomInRingOfSize(i, size):
                X[i, o + k] = 1.0
    return X


# ----------------------------------------------------------------- assembly
def load_jsonl_gz(p):
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


def load_atom_quarantine():
    """Records failing the I-1 atom-correspondence contract are excluded here so
    the exclusion is part of the data contract, not a manual step."""
    p = ROOT / "report" / "s7_l2b_r0r" / "I1_ATOM_QUARANTINE.json"
    if not p.is_file():
        return set()
    return set(json.loads(p.read_text(encoding="utf-8"))["quarantine_keys"])


def build(limit=0):
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")
    atom_quarantine = load_atom_quarantine()

    dev = load_jsonl_gz(CORPUS / "monn_development_edge_corpus.jsonl.gz")
    add = load_jsonl_gz(CORPUS / "monn_additional_pdb_edge_corpus.jsonl.gz")
    for r in dev:
        r["cohort"] = "development"
    for r in add:
        r["cohort"] = "additional_pdb"
    rows = dev + add
    if limit:
        rows = rows[:limit]

    md = [pickle.load((MONN / "mol_dict").open("rb"), encoding="bytes"),
          pickle.load((MONN / "independent_dataset_mol_dict").open("rb"),
                      encoding="bytes")]

    mol_cache, feat_cache, ident = {}, {}, {}
    contract = Counter()
    kept, quarantine = [], []

    for r in rows:
        if r["source_key"] in atom_quarantine:
            contract["atom_correspondence_quarantine"] += 1
            quarantine.append({"key": r["source_key"],
                               "reason": "I1_atom_correspondence_quarantine"})
            continue
        ccd = r["ligand_ccd"]
        mol = mol_cache.get(ccd)
        if mol is None:
            for d in md:
                m = d.get(ccd.encode("ascii", "ignore"))
                if m is not None:
                    mol = m
                    break
            if mol is None:
                contract["mol_missing"] += 1
                quarantine.append({"key": r["source_key"], "reason": "mol_missing"})
                continue
            mol_cache[ccd] = mol

        names = r["atom_names"]
        heavy_pos = [k for k, nm in enumerate(names) if not is_hydrogen_name(nm)]
        if len(heavy_pos) != mol.GetNumAtoms():
            contract["heavy_count_mismatch"] += 1
            quarantine.append({"key": r["source_key"], "reason": "heavy_count_mismatch",
                               "heavy": len(heavy_pos), "mol": mol.GetNumAtoms()})
            continue
        slot_to_mol = {s: m for m, s in enumerate(heavy_pos)}

        L = len(r["uniprot_sequence"])
        bad = False
        edges = []
        for res, slot in r["positive_binary_edges"]:
            if res < 0 or res >= L:
                contract["residue_index_out_of_range"] += 1
                bad = True
                break
            if slot not in slot_to_mol:
                contract["positive_edge_on_hydrogen_or_unknown"] += 1
                bad = True
                break
            edges.append((res, slot_to_mol[slot]))
        if bad:
            quarantine.append({"key": r["source_key"], "reason": "edge_index_contract"})
            continue
        if not edges:
            contract["no_edges_after_map"] += 1
            quarantine.append({"key": r["source_key"], "reason": "no_edges"})
            continue
        contract["ok"] += 1

        if ccd not in ident:
            w = Chem.Mol(mol)
            try:
                Chem.SanitizeMol(w)
                smi = Chem.MolToSmiles(w, isomericSmiles=False, canonical=True)
                ident[ccd] = (sha(smi.encode()),
                              MurckoScaffold.MurckoScaffoldSmiles(mol=w))
                feat_cache[ccd] = atom_features(w)
            except Exception:
                contract["rdkit_failed"] += 1
                quarantine.append({"key": r["source_key"], "reason": "rdkit_failed"})
                continue
        r["graph_key"], r["scaffold"] = ident[ccd]
        r["seq_key"] = sha(r["uniprot_sequence"].encode())
        r["edges"] = edges
        r["n_atoms"] = mol.GetNumAtoms()
        r["n_res"] = L
        kept.append(r)

    return kept, quarantine, dict(contract), feat_cache


def protein_components(rows):
    hom = json.loads((CACHE / "r0r2_closure" / "homology_alignment_stats.json").read_text())
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


def make_split(rows, comp):
    dev = [r for r in rows if r["cohort"] == "development"]
    by_comp = defaultdict(list)
    for r in dev:
        by_comp[comp[r["source_key"]]].append(r)
    # Component iteration order must match the registered split census
    # (report/s7_l2b_r0r/SPLIT_CONSTRUCTION_CENSUS.json): components are first
    # ordered by descending size with the component id as tiebreak, then
    # shuffled with SEED_SPLIT. The pre-sort is part of the frozen construction
    # because rng.shuffle permutes the given order.
    rng = np.random.default_rng(SEED_SPLIT)
    order = sorted(by_comp, key=lambda c: (-len(by_comp[c]), c))
    rng.shuffle(order)
    target = int(TEST_TARGET_FRACTION * len(dev))
    test_comps, n = set(), 0
    for c in order:
        if n >= target:
            break
        test_comps.add(c)
        n += len(by_comp[c])
    train = [r for r in dev if comp[r["source_key"]] not in test_comps]
    held = [r for r in dev if comp[r["source_key"]] in test_comps]
    tg = {r["graph_key"] for r in train}
    ts = {r["scaffold"] for r in train if r["scaffold"]}
    held_A = [r for r in held if r["graph_key"] not in tg]
    held_B = [r for r in held_A if not r["scaffold"] or r["scaffold"] not in ts]
    return train, held, held_A, held_B


if __name__ == "__main__":
    kept, quarantine, contract, feats = build()
    print("atom/residue contract census:", json.dumps(contract, indent=2))
    print(f"kept={len(kept)} quarantined={len(quarantine)}")
    comp = protein_components(kept)
    train, held, hA, hB = make_split(kept, comp)
    print(f"train={len(train)} heldout_all={len(held)} heldoutA={len(hA)} heldoutB={len(hB)}")
    print(f"train positives={sum(len(r['edges']) for r in train)} "
          f"heldoutA positives={sum(len(r['edges']) for r in hA)}")

"""C0 closure/admissibility gate and C1 audit-only exact-coupling information test.

Registered by
  research/correspondence_router/PREREG_C0_C1_CORRESPONDENCE_INFORMATION_AUDIT.md
  (sha256 007f8439..., commit f844679) committed BEFORE any C0 or C1 code.

Audit only. Trains nothing, introduces zero parameters, opens no affinity field
and never references heldout-A. The fixed-degree rewire is an EVALUATION NULL,
not a biological non-binder.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
HERE = ROOT / "research" / "correspondence_router"
sys.path.insert(0, str(HERE))

from c0_corpus import (  # noqa: E402
    CONTACT_THRESHOLD, EXEC, OUT, PREREG, PREREG_SHA, RAW, SENSITIVITY_THRESHOLD,
    SLOTS, C0ContractError, git_head, read_jsonl, sha_file, sha_json, write_json,
)

# ---- frozen contract (prereg sections 6, 7, 8) -----------------------------
IDENTITY_EDGE = 0.40
KMER_PREFILTER = 0.03
KMER_K = 3
PREFILTER_VALIDATION_SAMPLE = 200
MIN_COMPONENTS = 60
MAX_COMPONENT_FRACTION = 0.25
MAX_MDE = 0.05
C1A_MARGIN = 0.05
MIN_POSITIVE_UNITS = 10_000
MIN_CHECKERBOARDS = 1_000
SEED_BOOT = 20260903
SEED_NULL = 20260904
N_BOOT = 10_000
# Monte-Carlo replicate count for the null arm. Not registered with a value;
# it controls only the precision of the null's expectation, never its location.
N_REWIRE = 10
SWAPS_PER_EDGE = 100

CCD_DIR = RAW / "ccd"


def bootstrap(left: dict, right: dict, seed: int = SEED_BOOT, n_boot: int = N_BOOT):
    keys = sorted(set(left) & set(right))
    if not keys:
        return {"delta": float("nan"), "lcb95_one_sided": float("nan"), "units": 0}
    a = np.array([left[k] for k in keys])
    b = np.array([right[k] for k in keys])
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(keys), (n_boot, len(keys)))
    draws = (a[idx] - b[idx]).mean(1)
    return {"delta": float((a - b).mean()),
            "lcb95_one_sided": float(np.percentile(draws, 5)),
            "ci95": [float(np.percentile(draws, 2.5)),
                     float(np.percentile(draws, 97.5))],
            "units": len(keys)}


# --------------------------------------------------------------- exact AP
_HARM = np.concatenate(([0.0], np.cumsum(1.0 / np.arange(1, 4097, dtype=np.float64))))


def ap_exact(scores: np.ndarray, labels: np.ndarray):
    """Exact E[AP] under uniformly random ordering inside each tied block."""
    labels = np.asarray(labels)
    positives = float(labels.sum())
    if positives == 0 or positives == labels.size:
        return None
    order = np.argsort(-np.asarray(scores, dtype=np.float64), kind="stable")
    s = np.asarray(scores, dtype=np.float64)[order]
    y = labels[order].astype(np.float64)
    bounds = np.flatnonzero(np.r_[True, s[1:] != s[:-1], True])
    csum = np.concatenate(([0.0], np.cumsum(y)))
    lo, hi = bounds[:-1].astype(np.int64), bounds[1:].astype(np.int64)
    k = csum[hi] - csum[lo]
    mask = k > 0
    if not mask.any():
        return 0.0
    n = (hi - lo).astype(np.float64)[mask]
    a = lo.astype(np.float64)[mask]
    b = csum[lo][mask]
    k = k[mask]
    s1 = _HARM[hi[mask]] - _HARM[lo[mask]]
    grow = np.where(n > 1.0, (k - 1.0) / np.maximum(n - 1.0, 1.0), 0.0)
    return float((((k / n) * ((b + 1.0) * s1 + grow * (n - (a + 1.0) * s1))).sum())
                 / positives)


# --------------------------------------------------------------- ligand chemistry
def murcko_scaffolds(comp_ids) -> dict:
    """Registered selection rule S8, applied at closure time."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem.Scaffolds import MurckoScaffold
    RDLogger.DisableLog("rdApp.*")
    import gemmi
    out = {}
    for comp in sorted(comp_ids):
        path = CCD_DIR / f"{comp}.cif"
        if not path.is_file():
            continue
        try:
            block = gemmi.cif.read(str(path)).sole_block()
            names = [gemmi.cif.as_string(v) for v in
                     block.find_values("_chem_comp_atom.atom_id")]
            elements = [gemmi.cif.as_string(v) for v in
                        block.find_values("_chem_comp_atom.type_symbol")]
            molecule = Chem.RWMol()
            index = {}
            for name, element in zip(names, elements):
                if element.upper() == "H":
                    continue
                index[name] = molecule.AddAtom(Chem.Atom(element.title()))
            order = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
                     "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
            for one, two, kind in zip(
                    block.find_values("_chem_comp_bond.atom_id_1"),
                    block.find_values("_chem_comp_bond.atom_id_2"),
                    block.find_values("_chem_comp_bond.value_order")):
                a = index.get(gemmi.cif.as_string(one))
                b = index.get(gemmi.cif.as_string(two))
                if a is not None and b is not None:
                    molecule.AddBond(a, b, order.get(
                        gemmi.cif.as_string(kind).upper(), Chem.BondType.SINGLE))
            mol = molecule.GetMol()
            Chem.SanitizeMol(mol)
            out[comp] = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        except Exception:
            continue
    return out


# --------------------------------------------------------------- closure
def kmer_set(sequence: str, k: int = KMER_K):
    return {sequence[i:i + k] for i in range(len(sequence) - k + 1)}


def identity(left: str, right: str) -> float:
    import parasail
    result = parasail.nw_trace_striped_16(left, right, 10, 1, parasail.blosum62)
    traceback = result.traceback
    matches = sum(1 for a, b in zip(traceback.query, traceback.ref) if a == b and a != "-")
    return matches / max(1, min(len(left), len(right)))


class UnionFind:
    def __init__(self, items):
        self.parent = {i: i for i in items}

    def find(self, a):
        while self.parent[a] != a:
            self.parent[a] = self.parent[self.parent[a]]
            a = self.parent[a]
        return a

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def protein_clusters(sequences: list) -> tuple[dict, dict]:
    """Exact brute-force all-vs-all identity clustering, plus the registered
    3-mer prefilter's recall measured against it."""
    kmers = [kmer_set(s) for s in sequences]
    union = UnionFind(range(len(sequences)))
    edges, prefilter_kept, missed = 0, 0, 0
    for i in range(len(sequences)):
        for j in range(i + 1, len(sequences)):
            inter = len(kmers[i] & kmers[j])
            containment = inter / max(1, min(len(kmers[i]), len(kmers[j])))
            passes = containment >= KMER_PREFILTER
            if passes:
                prefilter_kept += 1
            value = identity(sequences[i], sequences[j])
            if value >= IDENTITY_EDGE:
                edges += 1
                union.union(i, j)
                if not passes:
                    missed += 1
    labels = {i: union.find(i) for i in range(len(sequences))}
    diagnostics = {"true_identity_edges": edges,
                   "prefilter_candidate_pairs": prefilter_kept,
                   "prefilter_missed_true_edges": missed,
                   "prefilter_recall": 1.0 if edges == 0 else (edges - missed) / edges,
                   "method": "exact brute-force all-vs-all; the registered 3-mer "
                             "prefilter was measured against it rather than relied on"}
    return labels, diagnostics


# --------------------------------------------------------------- nulls
def curveball(contact: np.ndarray, rng, swaps_per_edge: int = SWAPS_PER_EDGE):
    """Degree-preserving rewire by Curveball trades. Rows and columns keep their
    exact degrees; this is an evaluation null, not a biological non-binder."""
    rows = [set(np.flatnonzero(contact[i]).tolist()) for i in range(contact.shape[0])]
    edges = int(contact.sum())
    if edges == 0 or contact.shape[0] < 2:
        return contact.copy()
    for _ in range(max(1, swaps_per_edge * edges // 2)):
        i, j = rng.integers(0, len(rows), 2)
        if i == j:
            continue
        left, right = rows[i], rows[j]
        shared = left & right
        free = list((left | right) - shared)
        if not free:
            continue
        rng.shuffle(free)
        take = len(left) - len(shared)
        rows[i] = shared | set(free[:take])
        rows[j] = shared | set(free[take:])
    out = np.zeros_like(contact)
    for i, columns in enumerate(rows):
        if columns:
            out[i, list(columns)] = 1
    return out


def within_slot_ap(contact: np.ndarray, slot_of: np.ndarray,
                   candidate_slots) -> tuple[list, list]:
    """AP of ranking a slot's candidate residues by their contact degree, for
    every positive unit. Returns (all positive units, informative subset)."""
    marginal = contact.sum(0).astype(np.float64)
    every, informative = [], []
    for slot, columns in candidate_slots:
        sub = contact[:, columns]
        counts = sub.sum(1)
        scores = marginal[columns]
        for row in np.flatnonzero(counts >= 1):
            value = ap_exact(scores, sub[row])
            if value is None:                       # all candidates positive
                every.append(1.0)
                continue
            every.append(value)
            informative.append(value)
    return every, informative


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--limit", type=int, default=0)
    args = parser.parse_args(argv)
    started = time.time()
    try:
        result = run(limit=args.limit)
        print(json.dumps({k: v for k, v in result.items()
                          if k not in {"per_component"}}, indent=2, default=str),
              flush=True)
        print(f"elapsed {time.time() - started:.1f}s", flush=True)
        return 0
    except Exception as exc:
        failure = {"schema": "MetaSieve.Correspondence.C1.FailClosed.v1",
                   "error_type": type(exc).__name__, "error": str(exc),
                   "TERMINAL_VERDICT": "SLOT_ROUTING_ESTIMAND_INVALID",
                   "affinity_value_reads": 0}
        write_json(OUT / "C1_FAIL_CLOSED.json", failure)
        print(json.dumps(failure, indent=2), file=sys.stderr, flush=True)
        return 1


def run(limit: int = 0) -> dict:
    if sha_file(PREREG) != PREREG_SHA:
        raise C0ContractError("C0/C1 preregistration hash mismatch")
    census = json.loads((OUT / "C0_CORPUS_AND_CENSUS.json").read_text())
    systems = list(read_jsonl(EXEC / "c0_systems.jsonl"))
    geometry = np.load(EXEC / "c0_geometry.npz")
    if limit:
        systems = systems[:limit]

    # ---- registered selection rule S8, applied at closure time
    scaffolds = murcko_scaffolds({s["ligand_comp_id"] for s in systems})
    dropped_s8 = [s["system_id"] for s in systems
                  if s["ligand_comp_id"] not in scaffolds]
    systems = [s for s in systems if s["ligand_comp_id"] in scaffolds]
    if not systems:
        raise C0ContractError("no systems survive selection rule S8")

    # ---- closure
    sequences = sorted({s["sequence"] for s in systems})
    seq_index = {s: i for i, s in enumerate(sequences)}
    cluster_of, prefilter = protein_clusters(sequences)
    union = UnionFind([s["system_id"] for s in systems])
    by_protein_cluster = defaultdict(list)
    by_scaffold = defaultdict(list)
    for s in systems:
        s["protein_cluster"] = int(cluster_of[seq_index[s["sequence"]]])
        s["scaffold"] = scaffolds[s["ligand_comp_id"]]
        by_protein_cluster[s["protein_cluster"]].append(s["system_id"])
        by_scaffold[(s["scaffold"], s["ligand_comp_id"])].append(s["system_id"])
    for group in by_protein_cluster.values():        # E1 and E3
        for other in group[1:]:
            union.union(group[0], other)
    for group in by_scaffold.values():               # E2
        for other in group[1:]:
            union.union(group[0], other)
    component_of = {s["system_id"]: union.find(s["system_id"]) for s in systems}
    sizes = Counter(component_of.values())
    largest_fraction = max(sizes.values()) / len(systems)

    # ---- registered DataSAIL-style fallback
    scaffold_cluster = {}
    for index, key in enumerate(sorted({(s["scaffold"], s["ligand_comp_id"])
                                        for s in systems})):
        scaffold_cluster[key] = index
    blocks = defaultdict(list)
    for s in systems:
        blocks[(s["protein_cluster"],
                scaffold_cluster[(s["scaffold"], s["ligand_comp_id"])])].append(
            s["system_id"])
    fallback_sizes = Counter({k: len(v) for k, v in blocks.items()})
    fallback = {
        "blocks": len(blocks),
        "largest_block_fraction": max(fallback_sizes.values()) / len(systems),
        "note": "two-dimensional (protein cluster x ligand scaffold) blocks; "
                "conflicting systems would be discarded at split assignment",
    }

    use_fallback = not (len(sizes) >= MIN_COMPONENTS and
                        largest_fraction <= MAX_COMPONENT_FRACTION)
    if use_fallback:
        component_of = {sid: f"block_{a}_{b}"
                        for (a, b), group in blocks.items() for sid in group}
        sizes = Counter(component_of.values())
        largest_fraction = max(sizes.values()) / len(systems)

    g0a = len(sizes) >= MIN_COMPONENTS
    g0b = largest_fraction <= MAX_COMPONENT_FRACTION
    print(f"systems={len(systems)} components={len(sizes)} "
          f"largest_fraction={largest_fraction:.4f} fallback={use_fallback}", flush=True)

    # ---- per-system geometry, slots and arms
    rng = np.random.default_rng(SEED_NULL)
    per_system = {}
    totals = Counter()
    for s in systems:
        distances = geometry[f"d_{s['row']}"]
        indices = geometry[f"i_{s['row']}"]
        slot_of = np.minimum(SLOTS - 1, indices * SLOTS // s["sequence_length"])
        contact = (distances <= CONTACT_THRESHOLD).astype(np.int8)
        per_slot = Counter(slot_of.tolist())
        candidate_slots = [(slot, np.flatnonzero(slot_of == slot))
                           for slot, n in sorted(per_slot.items()) if n >= 2]
        if not candidate_slots or contact.sum() == 0:
            continue
        empirical, empirical_informative = within_slot_ap(contact, slot_of, candidate_slots)
        if not empirical:
            continue

        rewired, rewired_informative, degree_ok, moved = [], [], True, 0
        for _ in range(N_REWIRE):
            null = curveball(contact, rng)
            if not (np.array_equal(null.sum(0), contact.sum(0)) and
                    np.array_equal(null.sum(1), contact.sum(1))):
                degree_ok = False
            moved += int(np.abs(null - contact).sum())
            every, informative = within_slot_ap(null, slot_of, candidate_slots)
            rewired.extend(every)
            rewired_informative.extend(informative)

        atom_permutation = rng.permutation(contact.shape[0])
        atom_shuffled = contact[atom_permutation]
        slot_permuted = slot_of[rng.permutation(len(slot_of))]
        per_slot_p = Counter(slot_permuted.tolist())
        candidate_p = [(slot, np.flatnonzero(slot_permuted == slot))
                       for slot, n in sorted(per_slot_p.items()) if n >= 2]

        marginal = contact.sum(0).astype(np.float64)
        atom_marginal = contact.sum(1).astype(np.float64)
        complete = ap_exact(np.outer(atom_marginal, marginal).ravel(), contact.ravel())

        per_system[s["system_id"]] = {
            "component": component_of[s["system_id"]],
            "ap_empirical": float(np.mean(empirical)),
            "ap_rewire": float(np.mean(rewired)) if rewired else float("nan"),
            "ap_empirical_informative": (float(np.mean(empirical_informative))
                                         if empirical_informative else float("nan")),
            "ap_rewire_informative": (float(np.mean(rewired_informative))
                                      if rewired_informative else float("nan")),
            "ap_atom_shuffle": float(np.mean(
                within_slot_ap(atom_shuffled, slot_of, candidate_slots)[0] or [np.nan])),
            "ap_geometry_shuffle": float(np.mean(
                within_slot_ap(contact, slot_permuted, candidate_p)[0] or [np.nan])),
            "ap_complete_edge_marginal": complete,
            "positive_units": len(empirical),
            "informative_units": len(empirical_informative),
            "degree_preserved": degree_ok,
            "rewire_moved_entries": moved,
            "sequence": s["sequence"], "ligand_comp_id": s["ligand_comp_id"],
            "pdb_id": s["pdb_id"],
            "contacted_residues": np.flatnonzero(contact.any(0)).tolist(),
            "slot_of": slot_of.tolist(),
            "candidate_slots": [int(slot) for slot, _ in candidate_slots],
        }
        totals["positive_units"] += len(empirical)
        totals["informative_units"] += len(empirical_informative)

    if not per_system:
        raise C0ContractError("no scorable systems after closure")
    if not all(v["degree_preserved"] for v in per_system.values()):
        raise C0ContractError("fixed-degree rewire failed its degree contract")

    def component_macro(field):
        grouped = defaultdict(list)
        for value in per_system.values():
            if np.isfinite(value[field]):
                grouped[value["component"]].append(value[field])
        return {k: float(np.mean(v)) for k, v in grouped.items()}

    comp_emp = component_macro("ap_empirical")
    comp_rew = component_macro("ap_rewire")
    comp_emp_i = component_macro("ap_empirical_informative")
    comp_rew_i = component_macro("ap_rewire_informative")
    comp_atom = component_macro("ap_atom_shuffle")
    comp_geom = component_macro("ap_geometry_shuffle")
    comp_complete = component_macro("ap_complete_edge_marginal")

    # ---- G0c power, from the NULL arm dispersion only
    null_values = np.array(sorted(comp_rew.values()))
    sigma = float(np.std(null_values, ddof=1)) if null_values.size > 1 else float("nan")
    mde = float((1.645 + 0.842) * sigma / np.sqrt(max(len(null_values), 1)))
    g0c = bool(np.isfinite(mde) and mde <= MAX_MDE)
    print(f"null sigma={sigma:.4f} MDE={mde:.4f} G0c={g0c}", flush=True)

    # ---- replicate reliability (C1c)
    replicate_groups = defaultdict(list)
    for sid, value in per_system.items():
        replicate_groups[(value["sequence"], value["ligand_comp_id"])].append(sid)
    replicate_pairs, jaccards, ceilings = 0, [], {}
    for group in replicate_groups.values():
        entries = sorted({per_system[s]["pdb_id"] for s in group})
        if len(group) < 2 or len(entries) < 2:
            continue
        for a_pos in range(len(group)):
            for b_pos in range(a_pos + 1, len(group)):
                a, b = per_system[group[a_pos]], per_system[group[b_pos]]
                if a["pdb_id"] == b["pdb_id"]:
                    continue
                left = set(a["contacted_residues"])
                right = set(b["contacted_residues"])
                if not (left | right):
                    continue
                replicate_pairs += 1
                jaccards.append(len(left & right) / len(left | right))
                ceilings.setdefault(a["component"], []).append(
                    len(left & right) / len(left | right))
    replicate = {
        "replicate_systems": sum(len(g) for g in replicate_groups.values() if len(g) > 1),
        "replicate_pairs": replicate_pairs,
        "mean_contacted_residue_jaccard": (float(np.mean(jaccards)) if jaccards
                                           else float("nan")),
        "components_with_replicates": len(ceilings),
    }

    c1a = bootstrap(comp_emp, comp_rew)
    c1a["margin"] = C1A_MARGIN
    c1a["pass"] = bool(c1a["delta"] >= C1A_MARGIN and c1a["lcb95_one_sided"] > 0)
    c1a_informative = bootstrap(comp_emp_i, comp_rew_i)
    checkerboards = census["deconvolution_census"]["valid_checkerboards"]
    c1b = bool(totals["positive_units"] >= MIN_POSITIVE_UNITS and
               checkerboards >= MIN_CHECKERBOARDS)
    c1c = bool(replicate_pairs > 0 and np.isfinite(replicate["mean_contacted_residue_jaccard"]))

    if not (g0a and g0b and g0c):
        verdict = "CORRESPONDENCE_DATA_OR_CLOSURE_NOT_IDENTIFIABLE"
    elif not (c1a["pass"] and c1b and c1c):
        verdict = "EXACT_EDGE_COUPLING_NOT_SUPPORTED_BY_TEACHER"
    else:
        verdict = "CORRESPONDENCE_INFORMATION_PRESENT_C2_AUTHORIZED"

    macro = lambda d: float(np.mean(list(d.values()))) if d else float("nan")
    result = {
        "schema": "MetaSieve.Correspondence.C1.Gate.v1",
        "created_utc": "2026-08-10", "execution_commit": git_head(),
        "preregistration_sha256": PREREG_SHA,
        "s8_dropped_systems": len(dropped_s8),
        "panel": {
            "systems": len(per_system),
            "components": len(sizes),
            "largest_component_fraction": largest_fraction,
            "closure": "datasail_2d_fallback" if use_fallback else "union_E1_E2_E3",
            "union_components": len(Counter(union.find(s["system_id"])
                                            for s in systems)),
            "protein_clusters": len(set(cluster_of.values())),
            "distinct_scaffolds": len({s["scaffold"] for s in systems}),
            "common_mask_sha256": sha_json(sorted(per_system)),
        },
        "prefilter_validation": prefilter,
        "datasail_fallback": fallback,
        "C0_gates": {
            "G0a_components": {"observed": len(sizes), "required_at_least": MIN_COMPONENTS,
                               "pass": g0a},
            "G0b_largest_fraction": {"observed": largest_fraction,
                                     "required_at_most": MAX_COMPONENT_FRACTION,
                                     "pass": g0b},
            "G0c_power": {"null_sigma": sigma, "mde_80pct": mde,
                          "required_at_most": MAX_MDE, "pass": g0c,
                          "note": "sigma estimated from the degree-preserving null "
                                  "arm only, never from the tested contrast"},
        },
        "macro_ap_within_slot": {
            "empirical": macro(comp_emp),
            "fixed_degree_rewire_null": macro(comp_rew),
            "atom_shuffle": macro(comp_atom),
            "geometry_shuffle": macro(comp_geom),
            "empirical_informative_units_only": macro(comp_emp_i),
            "rewire_informative_units_only": macro(comp_rew_i),
        },
        "macro_ap_complete_edge_additive_marginal": macro(comp_complete),
        "C1_gates": {
            "C1a_empirical_minus_fixed_degree_rewire": c1a,
            "C1b_unit_sufficiency": {
                "positive_units": totals["positive_units"],
                "required_at_least": MIN_POSITIVE_UNITS,
                "valid_checkerboards": checkerboards,
                "checkerboards_required_at_least": MIN_CHECKERBOARDS,
                "pass": c1b},
            "C1c_replicate_ceiling_defined": {"pass": c1c, **replicate},
        },
        "non_gating": {
            "C1a_informative_units_only": c1a_informative,
            "informative_units": totals["informative_units"],
            "rewire_replicates": N_REWIRE,
            "swaps_per_edge": SWAPS_PER_EDGE,
            "sensitivity_threshold_angstrom": SENSITIVITY_THRESHOLD,
        },
        "rewire_is_an_evaluation_null_not_a_biological_nonbinder": True,
        "affinity_value_reads": 0,
        "heldoutA_referenced": False,
        "trainable_parameters_introduced": 0,
        "TERMINAL_VERDICT": verdict,
        "authorized_next_action": (
            "write the C2 preregistration for the geometry-gated router"
            if verdict == "CORRESPONDENCE_INFORMATION_PRESENT_C2_AUTHORIZED"
            else "none; stop at the earliest failed boundary"),
    }
    write_json(OUT / "C1_INFORMATION_AUDIT.json", result)
    return result


if __name__ == "__main__":
    raise SystemExit(main())

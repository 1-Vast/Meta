"""S4R-A frozen ligand representation availability and information audit.

Registered by
  research/s7_l2b_r0r/PREREG_PHASE2B_S4R_LIGAND_REPRESENTATION_AUDIT.md
  (sha256 8c3be169..., commit 8a643e8)
  research/s7_l2b_r0r/PREREG_PHASE2B_S4R_AUDIT_AMENDMENT_01.md
  (sha256 5210197f..., commit 1d89971)
both committed BEFORE this file existed.

Label-blind. No residue mask, no heldout-B, no affinity value and no S3R metric
is opened. The only ligand encoder used is the frozen RDKit Morgan generator
already present in the repository; no graph network is trained.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(r"D:\MetaSieve")
HERE = ROOT / "research" / "s7_l2b_r0r"
sys.path.insert(0, str(HERE))

from p2b_residue_residual import D_ATOM, D_ESM, g_of, sha_file  # noqa: E402
from s7_dataset import build, make_split, protein_components  # noqa: E402
from s7_localizer import FeatureStore  # noqa: E402
from s7_run import load_mols  # noqa: E402

PREREG = HERE / "PREREG_PHASE2B_S4R_LIGAND_REPRESENTATION_AUDIT.md"
AMENDMENT = HERE / "PREREG_PHASE2B_S4R_AUDIT_AMENDMENT_01.md"
PREREG_SHA = "8c3be16973957c6d1e7e735a7c3214d1d7e4f3b5d59f791e3f72271894130138"
AMENDMENT_SHA = "5210197fee00d7b15288f89925d8f282791305bd35b082f17bb5bd38b96745f4"

OUT = ROOT / "report" / "s7_l2b_r0r"
EXEC = ROOT / "dataset" / "processed" / "s7_l2b_r0r" / "phase2b_s4r_audit"

# ---- frozen candidate family (prereg section 4, amendment 01 section 2) ----
RADII = (1, 2)
VOCAB_SIZES = (128, 256, 512)
COSINE_COLLAPSE_TOL = 0.999
ZERO_TOL = 1e-12

# ---- frozen A-gates (prereg section 6) -------------------------------------
A1_EFFRANK_MULTIPLE = 3.0
A2_MIN_INCREMENTAL = 0.25
A3_MAX_RETENTION_LOSS = 0.10
A4_MIN_COVERAGE = 0.99


class S4RAuditError(RuntimeError):
    pass


def write_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True, default=str),
                    encoding="utf-8")


def git_head() -> str:
    return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT,
                                   text=True).strip()


def require_absent(paths) -> None:
    existing = [str(p) for p in paths if p.exists()]
    if existing:
        raise S4RAuditError(f"no-clobber contract: outputs already exist: {existing}")


# --------------------------------------------------------------- statistics
def effective_rank(matrix: np.ndarray) -> float:
    """exp(Shannon entropy of the normalized squared singular spectrum) of the
    mean-centred matrix. 1.0 for a rank-one matrix, D for an isotropic one."""
    centred = matrix - matrix.mean(axis=0)
    singular = np.linalg.svd(centred, compute_uv=False)
    power = singular ** 2
    total = power.sum()
    if total <= 0:
        return 0.0
    p = power / total
    return float(np.exp(-(p * np.log(p + 1e-300)).sum()))


def conditioning(matrix: np.ndarray, k: int) -> float:
    centred = matrix - matrix.mean(axis=0)
    singular = np.linalg.svd(centred, compute_uv=False)
    k = min(k, singular.size)
    if k == 0 or singular[k - 1] <= 0:
        return float("inf")
    return float(singular[0] / singular[k - 1])


def residual_fraction(target: np.ndarray, regressors: np.ndarray) -> float:
    """1 - R^2 of the least-squares regression of `target` onto [1, regressors],
    in the Frobenius sense over all target columns jointly."""
    design = np.concatenate([np.ones((regressors.shape[0], 1)), regressors], axis=1)
    coef, *_ = np.linalg.lstsq(design, target, rcond=None)
    residual = target - design @ coef
    denominator = float((target ** 2).sum())
    if denominator <= 0:
        return 0.0
    return float((residual ** 2).sum() / denominator)


def cosine_similarity(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    ln = np.linalg.norm(left, axis=1)
    rn = np.linalg.norm(right, axis=1)
    denominator = ln * rn
    out = np.zeros(left.shape[0], dtype=np.float64)
    ok = denominator > 0
    out[ok] = np.sum(left[ok] * right[ok], axis=1) / denominator[ok]
    return out


def exact_collapse(matrix: np.ndarray, keys) -> dict:
    groups = defaultdict(list)
    for i, key in enumerate(keys):
        groups[matrix[i].tobytes()].append(key)
    duplicated = [v for v in groups.values() if len(v) > 1]
    return {"identical_groups": len(duplicated),
            "graphs_covered": int(sum(len(v) for v in duplicated))}


# --------------------------------------------------------------- pair superset
def label_blind_graph_pairs(records):
    """Unordered distinct ligand-graph pairs inside one exact protein construct
    with both Murcko scaffolds present and distinct. Uses no residue label."""
    by_construct = defaultdict(list)
    for record in records:
        by_construct[record["seq_key"]].append(record)
    pairs = []
    for construct in sorted(by_construct):
        recs = sorted(by_construct[construct], key=lambda r: r["source_key"])
        seen = set()
        for i in range(len(recs)):
            for j in range(i + 1, len(recs)):
                a, b = recs[i], recs[j]
                if a["graph_key"] == b["graph_key"]:
                    continue
                sa, sb = a["scaffold"], b["scaffold"]
                if not (sa and sb and sa != sb):
                    continue
                key = (a["graph_key"], b["graph_key"])
                if key in seen:
                    continue
                seen.add(key)
                pairs.append(key)
    return pairs


# --------------------------------------------------------------- encoders
def morgan_environments(radius, graph_to_ccd, molecules):
    """Unfolded 32-bit Weisfeiler-Lehman environment counts (amendment 01)."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius)
    folded_generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius, fpSize=1 << 20)
    counts, folded_ids, failures = {}, set(), 0
    for graph_key in sorted(graph_to_ccd):
        working = Chem.Mol(molecules[graph_to_ccd[graph_key]])
        try:
            Chem.SanitizeMol(working)
        except Exception:
            counts[graph_key] = {}
            failures += 1
            continue
        counts[graph_key] = dict(
            generator.GetSparseCountFingerprint(working).GetNonzeroElements())
        folded_ids |= set(
            folded_generator.GetCountFingerprint(working).GetNonzeroElements())
    unfolded_ids = set()
    for value in counts.values():
        unfolded_ids |= set(value)
    return counts, {
        "n_environments_unfolded": len(unfolded_ids),
        "n_environments_folded_2p20": len(folded_ids),
        "folding_collisions": len(unfolded_ids) - len(folded_ids),
        "sanitize_failures": failures,
    }


def build_vocabulary(counts, train_graphs, size):
    document_frequency = Counter()
    for graph_key in sorted(train_graphs):
        for environment in counts.get(graph_key, {}):
            document_frequency[environment] += 1
    ranked = sorted(document_frequency.items(), key=lambda kv: (-kv[1], kv[0]))
    return [int(environment) for environment, _n in ranked[:size]]


def embed(counts, vocabulary, graph_keys, heavy_atoms):
    index = {environment: i for i, environment in enumerate(vocabulary)}
    matrix = np.zeros((len(graph_keys), len(vocabulary)), dtype=np.float64)
    for i, graph_key in enumerate(graph_keys):
        n = max(int(heavy_atoms[graph_key]), 1)
        for environment, count in counts.get(graph_key, {}).items():
            j = index.get(environment)
            if j is not None:
                matrix[i, j] = count / n
    return matrix


def vocabulary_sha(vocabulary) -> str:
    payload = json.dumps(vocabulary, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


# --------------------------------------------------------------- audit
def describe(matrix, diff, pair_left, pair_right, size_delta, log_size_delta,
             baseline_diff=None):
    norms = np.linalg.norm(diff, axis=1)
    cosines = cosine_similarity(matrix[pair_left], matrix[pair_right])
    stats = {
        "dimension": int(matrix.shape[1]),
        "effective_rank_embedding": effective_rank(matrix),
        "numerical_rank_embedding": int(np.linalg.matrix_rank(
            matrix - matrix.mean(axis=0))),
        "conditioning_sigma1_over_sigma_k": conditioning(matrix, 20),
        "effective_rank_difference": effective_rank(diff),
        "numerical_rank_difference": int(np.linalg.matrix_rank(
            diff - diff.mean(axis=0))),
        "zero_difference_rate": float(np.mean(norms <= ZERO_TOL)),
        "coverage": float(np.mean(norms > ZERO_TOL)),
        "near_identical_scaffold_distinct_rate": float(
            np.mean(cosines > COSINE_COLLAPSE_TOL)),
        "median_pair_cosine": float(np.median(cosines)),
        "r2_norm_on_atom_count_delta": 1.0 - residual_fraction(
            norms[:, None], size_delta[:, None]),
    }
    if baseline_diff is not None:
        stats["incremental_beyond_baseline"] = residual_fraction(diff, baseline_diff)
        stats["baseline_retention_loss"] = residual_fraction(baseline_diff, diff)
        stats["incremental_beyond_baseline_and_size"] = residual_fraction(
            diff, np.concatenate(
                [baseline_diff, size_delta[:, None], log_size_delta[:, None]], axis=1))
    return stats


def run() -> dict:
    require_absent([OUT / "PHASE2B_S4R_REPRESENTATION_AUDIT.json",
                    OUT / "PHASE2B_S4R_REPRESENTATION_AUDIT.md"])
    if sha_file(PREREG) != PREREG_SHA or sha_file(AMENDMENT) != AMENDMENT_SHA:
        raise S4RAuditError("S4R audit preregistration hash mismatch")

    kept, quarantine, contract, _features = build()
    component_of = protein_components(kept)
    train, held_all, held_a, held_b = make_split(kept, component_of)

    molecules, store = load_mols(), FeatureStore()
    graph_to_ccd, heavy_atoms, baseline = {}, {}, {}
    for record in kept:
        graph_key = record["graph_key"]
        if graph_key in baseline:
            continue
        graph_to_ccd[graph_key] = record["ligand_ccd"]
        atom_matrix = store.atoms(record, molecules)
        baseline[graph_key] = g_of(atom_matrix)
        heavy_atoms[graph_key] = int(atom_matrix.shape[0])

    graph_keys = sorted(baseline)
    position = {graph_key: i for i, graph_key in enumerate(graph_keys)}
    B = np.stack([baseline[graph_key] for graph_key in graph_keys])
    if B.shape[1] != D_ATOM:
        raise S4RAuditError("baseline ligand dimension changed")

    train_graphs = {r["graph_key"] for r in train}
    held_graphs = {r["graph_key"] for r in held_a}
    if train_graphs & held_graphs:
        raise S4RAuditError("ligand graph firewall failed")

    held_pairs = label_blind_graph_pairs(held_a)
    train_pairs = label_blind_graph_pairs(train)
    if not held_pairs:
        raise S4RAuditError("empty heldout-A label-blind pair superset")
    left = np.array([position[a] for a, _b in held_pairs])
    right = np.array([position[b] for _a, b in held_pairs])

    atom_counts = np.array([heavy_atoms[k] for k in graph_keys], dtype=np.float64)
    size_delta = np.abs(atom_counts[left] - atom_counts[right])
    log_size_delta = np.abs(np.log(atom_counts[left]) - np.log(atom_counts[right]))

    dB = B[left] - B[right]
    baseline_stats = describe(B, dB, left, right, size_delta, log_size_delta)
    baseline_stats["exact_collapse"] = exact_collapse(B, graph_keys)
    baseline_effective_rank = baseline_stats["effective_rank_difference"]
    a1_threshold = A1_EFFRANK_MULTIPLE * baseline_effective_rank

    inventory = {
        "records_kept": len(kept),
        "records_quarantined": len(quarantine),
        "contract_census": contract,
        "records_train": len(train),
        "records_heldout_all": len(held_all),
        "records_heldoutA": len(held_a),
        "records_heldoutB": len(held_b),
        "distinct_ligand_graphs_all": len(graph_keys),
        "distinct_ligand_graphs_train": len(train_graphs),
        "distinct_ligand_graphs_heldoutA": len(held_graphs),
        "distinct_scaffolds_all": len({r["scaffold"] for r in kept if r["scaffold"]}),
        "distinct_constructs_all": len({r["seq_key"] for r in kept}),
        "closure_components_all": len(set(component_of.values())),
        "label_blind_pairs_train": len(train_pairs),
        "label_blind_pairs_heldoutA": len(held_pairs),
    }

    grid, admissible, environment_census = [], [], {}
    vocabularies = {}
    for radius in RADII:
        counts, census = morgan_environments(radius, graph_to_ccd, molecules)
        train_environments = set()
        for graph_key in train_graphs:
            train_environments |= set(counts.get(graph_key, {}))
        held_environments = set()
        for graph_key in held_graphs:
            held_environments |= set(counts.get(graph_key, {}))
        census["environments_train"] = len(train_environments)
        census["environments_heldoutA"] = len(held_environments)
        census["environments_heldoutA_absent_from_train"] = len(
            held_environments - train_environments)
        environment_census[f"radius_{radius}"] = census

        for size in VOCAB_SIZES:
            vocabulary = build_vocabulary(counts, train_graphs, size)
            if len(vocabulary) != size:
                raise S4RAuditError(
                    f"train vocabulary for radius={radius} has {len(vocabulary)} "
                    f"entries, expected {size}")
            M = embed(counts, vocabulary, graph_keys, heavy_atoms)
            dM = M[left] - M[right]
            stats = describe(M, dM, left, right, size_delta, log_size_delta,
                             baseline_diff=dB)
            stats["exact_collapse"] = exact_collapse(M, graph_keys)
            stats["radius"] = radius
            stats["vocabulary_size"] = size
            stats["vocabulary_sha256"] = vocabulary_sha(vocabulary)
            stats["implied_W_parameters"] = D_ESM * size
            gates = {
                "A1_effective_rank": {
                    "observed": stats["effective_rank_difference"],
                    "required_at_least": a1_threshold,
                    "pass": stats["effective_rank_difference"] >= a1_threshold},
                "A2_incremental_beyond_baseline": {
                    "observed": stats["incremental_beyond_baseline"],
                    "required_at_least": A2_MIN_INCREMENTAL,
                    "pass": stats["incremental_beyond_baseline"] >= A2_MIN_INCREMENTAL},
                "A3_baseline_retention_loss": {
                    "observed": stats["baseline_retention_loss"],
                    "required_at_most": A3_MAX_RETENTION_LOSS,
                    "pass": stats["baseline_retention_loss"] <= A3_MAX_RETENTION_LOSS},
                "A4_coverage": {
                    "observed": stats["coverage"],
                    "required_at_least": A4_MIN_COVERAGE,
                    "pass": stats["coverage"] >= A4_MIN_COVERAGE},
            }
            stats["A_gates"] = gates
            stats["admissible"] = all(g["pass"] for g in gates.values())
            grid.append(stats)
            vocabularies[(radius, size)] = vocabulary
            if stats["admissible"]:
                admissible.append((size, radius))
            print(f"radius={radius} d={size} effrank={stats['effective_rank_difference']:.2f} "
                  f"inc={stats['incremental_beyond_baseline']:.4f} "
                  f"ret={stats['baseline_retention_loss']:.4f} "
                  f"cov={stats['coverage']:.5f} admissible={stats['admissible']}",
                  flush=True)

    if not admissible:
        verdict = "GRAPH_LIGAND_REPRESENTATION_NOT_INFORMATIVE"
        selected = None
    else:
        verdict = "GRAPH_LIGAND_REPRESENTATION_AVAILABLE_AND_INFORMATIVE"
        size, radius = min(admissible)          # lexicographic (d, r), prereg s.7
        selected = {"radius": radius, "vocabulary_size": size}

    EXEC.mkdir(parents=True, exist_ok=True)
    selection_artifacts = {}
    if selected is not None:
        vocabulary = vocabularies[(selected["radius"], selected["vocabulary_size"])]
        path = EXEC / "selected_ligand_vocabulary.json"
        write_json(path, {
            "radius": selected["radius"],
            "vocabulary_size": selected["vocabulary_size"],
            "normalization": "count divided by heavy-atom count",
            "identifier_space": "unfolded 32-bit RDKit Morgan environment hash",
            "source": "top-d by distinct train-split ligand graphs, ties by ascending id",
            "vocabulary": vocabulary,
        })
        selected["vocabulary_sha256"] = vocabulary_sha(vocabulary)
        selected["implied_W_parameters"] = D_ESM * selected["vocabulary_size"]
        selected["baseline_W_parameters"] = D_ESM * D_ATOM
        selection_artifacts = {
            str(path.relative_to(ROOT)).replace("\\", "/"): sha_file(path)}

    result = {
        "schema": "MetaSieve.S7L2B.P2B.S4R.RepresentationAudit.v1",
        "created_utc": "2026-08-10",
        "execution_commit": git_head(),
        "preregistration_sha256": PREREG_SHA,
        "amendment_01_sha256": AMENDMENT_SHA,
        "inventory": inventory,
        "environment_census": environment_census,
        "baseline_41d_mean": baseline_stats,
        "A1_threshold_from_baseline": a1_threshold,
        "grid": grid,
        "selection_rule": "smallest (vocabulary_size, radius) among admissible",
        "selected": selected,
        "artifacts": selection_artifacts,
        "residue_label_reads": 0,
        "affinity_value_reads": 0,
        "heldoutB_reads": 0,
        "s3r_metric_reads": 0,
        "trainable_parameters_introduced": 0,
        "TERMINAL_VERDICT": verdict,
        "authorized_next_action": (
            "write the S4R training preregistration for the selected (radius, d) only"
            if selected is not None else
            "none; stop with GRAPH_LIGAND_REPRESENTATION_NOT_INFORMATIVE"),
    }
    write_json(OUT / "PHASE2B_S4R_REPRESENTATION_AUDIT.json", result)
    write_report(result)
    return result


def write_report(result) -> None:
    base = result["baseline_41d_mean"]
    lines = [
        "# Phase 2B S4R-A ligand representation audit", "",
        f"Terminal verdict: `{result['TERMINAL_VERDICT']}`", "",
        "Label-blind. Residue label reads: "
        f"{result['residue_label_reads']}. Affinity value reads: "
        f"{result['affinity_value_reads']}. Trainable parameters introduced: "
        f"{result['trainable_parameters_introduced']}.", "",
        "## Baseline mean-pooled 41-D reference", "",
        "| statistic | value |", "|---|---:|",
        f"| effective rank of the embedding | {base['effective_rank_embedding']:.3f} |",
        f"| numerical rank of the embedding | {base['numerical_rank_embedding']} |",
        f"| effective rank of pair differences | {base['effective_rank_difference']:.3f} |",
        f"| distinct graphs sharing an identical vector | "
        f"{base['exact_collapse']['graphs_covered']} in "
        f"{base['exact_collapse']['identical_groups']} groups |",
        f"| scaffold-distinct pairs with cosine > 0.999 | "
        f"{base['near_identical_scaffold_distinct_rate']:.5f} |",
        f"| R2 of difference norm on heavy-atom-count delta | "
        f"{base['r2_norm_on_atom_count_delta']:.4f} |",
        "",
        f"A1 therefore requires a candidate difference effective rank of at least "
        f"`{result['A1_threshold_from_baseline']:.3f}`.", "",
        "## Candidate grid", "",
        "| radius | d | eff.rank(dg) | INC | RET | coverage | INC_perp | W params | admissible |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|:---:|",
    ]
    for row in result["grid"]:
        lines.append(
            f"| {row['radius']} | {row['vocabulary_size']} | "
            f"{row['effective_rank_difference']:.2f} | "
            f"{row['incremental_beyond_baseline']:.4f} | "
            f"{row['baseline_retention_loss']:.4f} | "
            f"{row['coverage']:.5f} | "
            f"{row['incremental_beyond_baseline_and_size']:.4f} | "
            f"{row['implied_W_parameters']} | {row['admissible']} |")
    lines += ["", "`INC` is the fraction of candidate pair-difference energy that no",
              "linear function of the baseline pair difference can express. `RET` is",
              "the converse loss. `INC_perp` additionally removes heavy-atom-count and",
              "log-count differences. All are heldout-A, label-blind.", ""]
    if result["selected"] is not None:
        selected = result["selected"]
        lines += [
            "## Selection", "",
            f"Capacity-parsimony rule selects radius `{selected['radius']}`, "
            f"vocabulary `{selected['vocabulary_size']}`, "
            f"`W` parameters `{selected['implied_W_parameters']}` against the "
            f"baseline's `{selected['baseline_W_parameters']}`.", "",
            f"Frozen vocabulary SHA-256 `{selected['vocabulary_sha256']}`.", ""]
    else:
        lines += ["## Selection", "", "No candidate cleared every A-gate.", ""]
    lines += ["## Boundary", "",
              "This audit measures a ligand representation only. It admits no",
              "statistic to `z`, opens no residue label, and does not modify",
              "`A(F,z) = K(B(z)F(z))`.", ""]
    (OUT / "PHASE2B_S4R_REPRESENTATION_AUDIT.md").write_text(
        "\n".join(lines), encoding="utf-8")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", nargs="?", default="run", choices=("run",))
    parser.parse_args(argv)
    started = time.time()
    try:
        result = run()
        print(json.dumps({"TERMINAL_VERDICT": result["TERMINAL_VERDICT"],
                          "selected": result["selected"],
                          "elapsed_seconds": round(time.time() - started, 3)},
                         indent=2, default=str), flush=True)
        return 0
    except Exception as exc:
        failure = {
            "schema": "MetaSieve.S7L2B.P2B.S4R.AuditFailClosed.v1",
            "created_utc": "2026-08-10",
            "error_type": type(exc).__name__, "error": str(exc),
            "TERMINAL_VERDICT": "S4R_AUDIT_CONTRACT_FAIL_CLOSED",
            "residue_label_reads": 0, "affinity_value_reads": 0,
        }
        write_json(OUT / "PHASE2B_S4R_AUDIT_FAIL_CLOSED.json", failure)
        print(json.dumps(failure, indent=2), file=sys.stderr, flush=True)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

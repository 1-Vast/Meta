"""Stage W0 — soft controlled chemical-change surface census (Davis and KIBA).

Run:
    python -m tools.research.stageW_soft_mmp.w0_census
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from tools.research.stageU_mmp_interaction.mmp import fragment, transformation
from tools.research.stageW_soft_mmp.soft_mmp import soft_transformation

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DAVIS = ROOT / "dataset/raw/dta/davis.tab"
KIBA = ROOT / "dataset/raw/dta/kiba.tab"
CDHIT = ROOT / "tools/runtime/cdhit/4.8.1/cd-hit.exe"
PREREG_SHA = "ae96762e319521f30aa09eb1a79fb8bb0e3ea324b21d4b40868aa6826a45dc71"

GATES = {
    "same_target_observations": 1000,
    "targets": 30,
    "components": 10,
    "rich_families": 20,
    "cross_component_D_rows": 500,
    "top1_family_share": 0.10,
    "top10_family_share": 0.30,
    "top1_target_share": 0.20,
    "top5_target_share": 0.60,
    "same_core_residual_median": 1.00,
}


def file_sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_dataset(path: Path) -> tuple[dict, dict, dict]:
    """Return rows, targets {seq_hash: sequence}, ligands {smiles: id}."""
    targets: dict[str, str] = {}
    ligands: dict[str, str] = {}
    rows: dict[tuple[str, str], list[float]] = defaultdict(list)
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if len(header) < 5:
            raise ValueError(f"unexpected header in {path}: {header}")
        for line in reader:
            if len(line) < 5:
                continue
            smiles = line[1].strip()
            sequence = line[3].strip()
            try:
                affinity = float(line[4])
            except ValueError:
                continue
            if not smiles or not sequence:
                continue
            target_id = hashlib.sha256(sequence.encode("utf-8")).hexdigest()
            targets[target_id] = sequence
            ligand_id = hashlib.sha256(smiles.encode("utf-8")).hexdigest()
            ligands[ligand_id] = smiles
            rows[(target_id, ligand_id)].append(affinity)
    return dict(rows), targets, ligands


def cdhit_components(targets: dict[str, str], dataset: str) -> dict[str, str]:
    """CD-HIT 40 percent sequence components; component id = cluster id."""
    ordered = sorted(targets)
    short_ids = {target_id: f"s{index}" for index, target_id in enumerate(ordered)}
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        fasta = tmp_path / "targets.fasta"
        with fasta.open("w", encoding="utf-8") as handle:
            for target_id in ordered:
                handle.write(f">{short_ids[target_id]}\n{targets[target_id]}\n")
        out = tmp_path / "out"
        command = [str(CDHIT), "-i", str(fasta), "-o", str(out),
                   "-c", "0.4", "-n", "2", "-G", "1", "-g", "1",
                   "-M", "8000", "-T", "1"]
        subprocess.run(command, check=True, capture_output=True)
        cluster_file = Path(str(out) + ".clstr")
        short_assignment: dict[str, str] = {}
        current = None
        for line in cluster_file.read_text(encoding="utf-8").splitlines():
            if line.startswith(">"):
                current = line[1:].strip()
                continue
            if ">" in line:
                member = line.split(">", 1)[1].split()[0].replace("...", "").strip()
                short_assignment[member] = current
        missing = set(short_ids.values()) - set(short_assignment)
        if missing:
            raise RuntimeError(f"cd-hit missed {len(missing)} sequences")
        return {target_id: short_assignment[short_ids[target_id]]
                for target_id in ordered}


def build_observations(rows, targets, ligands, components):
    """Same-target MMP observations with soft family keys."""
    ligand_of: dict[str, str] = {v: k for k, v in ligands.items()}
    observations = []
    no_cut = 0
    for target in sorted(targets):
        by_target: dict[str, list[float]] = defaultdict(list)
        for (tid, lid), values in rows.items():
            if tid == target:
                by_target[lid].extend(values)
        fragments: dict[str, tuple] = {}
        for ligand_id, values in by_target.items():
            smiles = ligands[ligand_id]
            pieces = fragment(smiles)
            if not pieces:
                no_cut += 1
                continue
            fragments[ligand_id] = pieces
        by_core: dict[str, list[tuple[str, object]]] = defaultdict(list)
        for ligand_id, pieces in fragments.items():
            for piece in pieces:
                by_core[piece.core].append((ligand_id, piece))
        emitted: set[tuple[str, str, str]] = set()
        for core, entries in sorted(by_core.items()):
            for pos, (left_lid, left) in enumerate(entries):
                for right_lid, right in entries[pos + 1:]:
                    built = transformation(left, right)
                    if built is None:
                        continue
                    item, _exact_flip = built
                    signature = (core, item.r_a, item.r_b)
                    if signature in emitted:
                        continue
                    emitted.add(signature)
                    soft = soft_transformation(left, right)
                    if soft is None:
                        continue
                    # Exact canonical direction is from transformation(); use
                    # its flip for label direction and the soft family key.
                    left_y = float(np.median(by_target[left_lid]))
                    right_y = float(np.median(by_target[right_lid]))
                    delta = right_y - left_y
                    if soft.flipped:
                        delta = -delta
                    observations.append({
                        "target": target,
                        "component": components[target],
                        "exact_key": item.exact_key,
                        "family_key": soft.family_key,
                        "core": core,
                        "delta_y": delta,
                        "murcko_core": soft.murcko_core,
                        "category_a": list(soft.category_a),
                        "category_b": list(soft.category_b),
                    })
    return observations, no_cut


def family_census(observations) -> dict:
    key_targets: dict[str, set[str]] = defaultdict(set)
    key_components: dict[str, set[str]] = defaultdict(set)
    key_counts: Counter = Counter()
    target_counts: Counter = Counter()
    for item in observations:
        key = item["family_key"]
        key_targets[key].add(item["target"])
        key_components[key].add(item["component"])
        key_counts[key] += 1
        target_counts[item["target"]] += 1
    total = len(observations)
    rich = [key for key in key_counts
            if len(key_targets[key]) >= 3 and len(key_components[key]) >= 3]
    target_component = {item["target"]: item["component"]
                        for item in observations}
    cross_d = 0
    for key in key_targets:
        ordered = sorted(key_targets[key])
        for i, left in enumerate(ordered):
            for right in ordered[i + 1:]:
                if target_component[left] != target_component[right]:
                    cross_d += 1
    def share(counter, n):
        return sum(v for _, v in counter.most_common(n)) / total if total else 0.0
    return {
        "observations": total,
        "targets": len(target_counts),
        "components": len({o["component"] for o in observations}),
        "exact_keys": len({o["exact_key"] for o in observations}),
        "family_keys": len(key_counts),
        "rich_families": len(rich),
        "cross_component_D_rows": cross_d,
        "top1_family_share": share(key_counts, 1),
        "top10_family_share": share(key_counts, 10),
        "top1_target_share": share(target_counts, 1),
        "top5_target_share": share(target_counts, 5),
    }


def same_core_residual(observations) -> dict:
    groups: dict[tuple[str, str], list[float]] = defaultdict(list)
    for item in observations:
        groups[(item["target"], item["family_key"])].append(item["delta_y"])
    residuals = []
    multi_core_groups = 0
    for (target, key), values in groups.items():
        if len(values) >= 2:
            multi_core_groups += 1
            for i, left in enumerate(values):
                for right in values[i + 1:]:
                    residuals.append(abs(left - right))
    if not residuals:
        return {"groups": 0, "median": None, "p95": None}
    array = np.asarray(residuals)
    return {"groups": multi_core_groups, "pairs": len(residuals),
            "median": float(np.median(array)),
            "p95": float(np.quantile(array, 0.95)),
            "max": float(array.max())}


def random_protein_null(observations, draws=500, seed=20260820) -> dict:
    """Pooled between-component mean square, with a component-permuted null.

    The previous within-family total-MS null was permutation-invariant and
    therefore useless (observed == every draw). The corrected statistic is
    the pooled, within-family between-component variation of delta_y; the null
    permutes protein-component labels across observations, which destroys any
    protein x transformation association while preserving the family and
    marginal label distributions.
    """
    keys = sorted({item["family_key"] for item in observations})
    components = sorted({item["component"] for item in observations})
    key_index = {key: i for i, key in enumerate(keys)}
    component_index = {comp: i for i, comp in enumerate(components)}
    key_codes = np.asarray([key_index[item["family_key"]]
                            for item in observations], dtype=np.int64)
    comp_codes = np.asarray([component_index[item["component"]]
                             for item in observations], dtype=np.int64)
    delta = np.asarray([item["delta_y"] for item in observations],
                       dtype=np.float64)
    n_keys, n_components = len(keys), len(components)
    family_count = np.bincount(key_codes, minlength=n_keys)
    family_sum = np.bincount(key_codes, weights=delta, minlength=n_keys)
    with np.errstate(divide="ignore", invalid="ignore"):
        family_mean = family_sum / np.maximum(family_count, 1)

    def pooled_between_component_ms(comp_labels):
        flat = key_codes * n_components + comp_labels
        counts = np.bincount(flat, minlength=n_keys * n_components)
        sums = np.bincount(flat, weights=delta,
                           minlength=n_keys * n_components)
        with np.errstate(divide="ignore", invalid="ignore"):
            comp_mean = sums / np.maximum(counts, 1)
        cell = np.flatnonzero(counts > 0)
        diff = comp_mean[cell] - family_mean[cell // n_components]
        numerator = float((counts[cell] * diff * diff).sum())
        denominator = float(cell.size - n_keys)
        return numerator / max(denominator, 1.0)

    observed = pooled_between_component_ms(comp_codes)
    rng = np.random.default_rng(seed)
    nulls = np.empty(draws, dtype=np.float64)
    for draw in range(draws):
        nulls[draw] = pooled_between_component_ms(
            rng.permutation(comp_codes))
    return {
        "observed_between_component_MS": observed,
        "null_MS_median": float(np.median(nulls)),
        "null_MS_q05": float(np.quantile(nulls, 0.05)),
        "null_MS_q95": float(np.quantile(nulls, 0.95)),
        "p_permutation_upper": float((nulls >= observed).mean()),
        "draws": draws, "seed": seed,
        "definition": ("pooled within-family between-component variance of "
                       "delta_y; null permutes protein-component labels"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        default=HERE / "W0_SOFT_MMP_CENSUS.json")
    args = parser.parse_args()

    report = {
        "schema": "MetaSieve.StageW.W0SoftMMPCensus.v1",
        "stage": "stageW_soft_mmp",
        "preregistration_sha256": PREREG_SHA,
        "datasets": {},
        "gates": GATES,
    }
    for dataset, path in (("davis", DAVIS), ("kiba", KIBA)):
        rows, targets, ligands = read_dataset(path)
        components = cdhit_components(targets, dataset)
        observations, no_cut = build_observations(
            rows, targets, ligands, components)
        census = family_census(observations)
        residual = same_core_residual(observations)
        null = random_protein_null(observations)
        checks = {
            "same_target_observations": census["observations"],
            "targets": census["targets"],
            "components": census["components"],
            "rich_families": census["rich_families"],
            "cross_component_D_rows": census["cross_component_D_rows"],
            "top1_family_share": census["top1_family_share"],
            "top10_family_share": census["top10_family_share"],
            "top1_target_share": census["top1_target_share"],
            "top5_target_share": census["top5_target_share"],
            "same_core_residual_median": residual["median"],
        }
        gate = {}
        for name, value in GATES.items():
            measured = checks[name]
            share_gate = name in {
                "top1_family_share", "top10_family_share",
                "top1_target_share", "top5_target_share"} or                 name.startswith("same_core")
            if measured is None:
                gate[name] = {"measured": measured, "threshold": value,
                              "pass": False}
            else:
                gate[name] = {
                    "measured": measured, "threshold": value,
                    "pass": measured <= value if share_gate else
                            measured >= value}
        report["datasets"][dataset] = {
            "file": str(path),
            "file_sha256": file_sha(path),
            "unique_target_sequences": len(targets),
            "unique_ligand_smiles": len(ligands),
            "unique_rows": len(rows),
            "ligands_with_no_cut": no_cut,
            "census": census,
            "same_core_residual": residual,
            "random_protein_null": null,
            "gate": gate,
            "all_pass": all(item["pass"] for item in gate.values()),
        }

    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "output": str(args.output),
        "davis": {k: v["pass"] for k, v in
                  report["datasets"]["davis"]["gate"].items()},
        "davis_all_pass": report["datasets"]["davis"]["all_pass"],
        "kiba": {k: v["pass"] for k, v in
                 report["datasets"]["kiba"]["gate"].items()},
        "kiba_all_pass": report["datasets"]["kiba"]["all_pass"],
    }, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

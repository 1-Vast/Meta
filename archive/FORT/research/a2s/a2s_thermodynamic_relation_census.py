"""Source-only relation coverage gate for thermodynamic-cycle A2S adaptation.

This module does not train an adapter.  It asks whether a passive k-shot
episode contains enough fit-vocabulary matched molecular transformations to
support an explicit SAR grammar, and whether fit-estimated transformation
effects transfer to held-component probe targets.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator, rdMMPA

from research.a2s.a2s_assay_coherence_gate import K_VALUES, MIN_QUERY, build_episodes
from research.a2s.a2s_trace import DEFAULT_LOCK, DEFAULT_OOF, Substrate, load_substrate
from research.a2s.a2s_trace_stratum import paired_bootstrap


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "a2s_thermodynamic_relation_census_2026-08-02.json"
DEFAULT_FRAGMENTS = ROOT / "reports" / "active" / "a2s_thermodynamic_relation_fragments_2026-08-02.parquet"
DEFAULT_EDGES = ROOT / "reports" / "active" / "a2s_thermodynamic_relation_edges_2026-08-02.parquet"
DEFAULT_COVERAGE = ROOT / "reports" / "active" / "a2s_thermodynamic_relation_coverage_2026-08-02.parquet"

MMP_PATTERN = "[#6+0;!$(*=,#[!#6])]!@!=!#[*]"
MIN_CORE_HEAVY = 6
MIN_CORE_FRACTION = 0.67
MAX_SUBSTITUENT_HEAVY = 12
MIN_PRIOR_COMPONENTS = 3
MIN_DIRECTION_EFFECT = 0.25
LOW_SIMILARITY = 0.35
POWERED_COMPONENTS = 47
BOOTSTRAP_DRAWS = 2000


@dataclass(frozen=True, order=True)
class Fragment:
    core: str
    substituent: str
    core_heavy: int
    substituent_heavy: int


def canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def transform_key(substituent_a: str, substituent_b: str) -> str:
    ordered = sorted((substituent_a, substituent_b))
    return sha256(canonical(ordered).encode()).hexdigest()[:20]


def orient_replacement(
    substituent_a: str,
    value_a: float,
    substituent_b: str,
    value_b: float,
) -> tuple[str, float, str, float]:
    if substituent_a <= substituent_b:
        return substituent_a, value_a, substituent_b, value_b
    return substituent_b, value_b, substituent_a, value_a


def fragment_smiles(smiles: str) -> tuple[Fragment, ...]:
    mol = Chem.MolFromSmiles(str(smiles))
    if mol is None:
        return ()
    total_heavy = mol.GetNumHeavyAtoms()
    if total_heavy < MIN_CORE_HEAVY + 1:
        return ()
    records: set[Fragment] = set()
    cuts = rdMMPA.FragmentMol(mol, 1, 1, 100, MMP_PATTERN, False)
    for _, fragments in cuts:
        parts = fragments.split(".")
        if len(parts) != 2:
            continue
        part_mols = [Chem.MolFromSmiles(part) for part in parts]
        if any(part is None for part in part_mols):
            continue
        heavy = [part.GetNumHeavyAtoms() for part in part_mols]
        if heavy[0] == heavy[1]:
            continue
        core_index = int(np.argmax(heavy))
        substituent_index = 1 - core_index
        core_heavy = heavy[core_index]
        substituent_heavy = heavy[substituent_index]
        if core_heavy < MIN_CORE_HEAVY:
            continue
        if substituent_heavy < 1 or substituent_heavy > MAX_SUBSTITUENT_HEAVY:
            continue
        if core_heavy / total_heavy < MIN_CORE_FRACTION:
            continue
        records.add(
            Fragment(
                core=Chem.MolToSmiles(part_mols[core_index], canonical=True),
                substituent=Chem.MolToSmiles(part_mols[substituent_index], canonical=True),
                core_heavy=core_heavy,
                substituent_heavy=substituent_heavy,
            )
        )
    return tuple(sorted(records))


def build_fragment_frame(smiles_values: Iterable[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for smiles in sorted(set(str(value) for value in smiles_values)):
        for fragment in fragment_smiles(smiles):
            rows.append({"conn": smiles, **asdict(fragment)})
    return pd.DataFrame.from_records(
        rows,
        columns=("conn", "core", "substituent", "core_heavy", "substituent_heavy"),
    )


def aggregate_measurements(substrate: Substrate) -> pd.DataFrame:
    frame = substrate.labeled.copy()
    frame["base"] = substrate.base.detach().cpu().numpy().astype(np.float64)
    return (
        frame.groupby(["role", "component", "target", "assays", "conn"], as_index=False)
        .agg(affinity=("affinity", "median"), base=("base", "median"))
        .sort_values(["role", "component", "target", "assays", "conn"])
        .reset_index(drop=True)
    )


def fingerprint_cache(smiles_values: Iterable[str]) -> dict[str, object]:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    output: dict[str, object] = {}
    for smiles in sorted(set(str(value) for value in smiles_values)):
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            output[smiles] = generator.GetFingerprint(mol)
    return output


def build_mmp_edges(measurements: pd.DataFrame, fragments: pd.DataFrame) -> pd.DataFrame:
    expanded = measurements.merge(fragments, on="conn", how="inner", validate="many_to_many")
    fingerprints = fingerprint_cache(expanded.conn.unique())
    rows: list[dict[str, object]] = []
    group_columns = ["role", "component", "target", "assays", "core"]
    for keys, group in expanded.groupby(group_columns, sort=True):
        role, component, target, assay, core = keys
        group = group.sort_values(["substituent", "conn"]).drop_duplicates("conn")
        for (_, left), (_, right) in combinations(group.iterrows(), 2):
            if left.substituent == right.substituent:
                continue
            sub_a, affinity_a, sub_b, affinity_b = orient_replacement(
                str(left.substituent),
                float(left.affinity),
                str(right.substituent),
                float(right.affinity),
            )
            if str(left.substituent) == sub_a:
                row_a, row_b = left, right
            else:
                row_a, row_b = right, left
            fp_a = fingerprints.get(str(row_a.conn))
            fp_b = fingerprints.get(str(row_b.conn))
            if fp_a is None or fp_b is None:
                continue
            pair_id = sha256(canonical(sorted((str(row_a.conn), str(row_b.conn)))).encode()).hexdigest()[:20]
            rows.append(
                {
                    "role": str(role),
                    "component": str(component),
                    "target": str(target),
                    "assay": str(assay),
                    "core": str(core),
                    "core_heavy": int(min(row_a.core_heavy, row_b.core_heavy)),
                    "substituent_a": sub_a,
                    "substituent_b": sub_b,
                    "transformation": transform_key(sub_a, sub_b),
                    "pair_id": pair_id,
                    "conn_a": str(row_a.conn),
                    "conn_b": str(row_b.conn),
                    "delta_affinity": float(affinity_b - affinity_a),
                    "delta_base": float(row_b.base - row_a.base),
                    "delta_residual": float(
                        (affinity_b - row_b.base) - (affinity_a - row_a.base)
                    ),
                    "tanimoto": float(DataStructs.TanimotoSimilarity(fp_a, fp_b)),
                }
            )
    edges = pd.DataFrame.from_records(rows)
    if edges.empty:
        return edges
    return (
        edges.sort_values("core_heavy", ascending=False)
        .drop_duplicates(["role", "target", "assay", "pair_id"])
        .sort_values(["role", "component", "target", "assay", "transformation", "pair_id"])
        .reset_index(drop=True)
    )


def fit_transformation_priors(edges: pd.DataFrame) -> pd.DataFrame:
    fit = edges.loc[edges.role == "fit"]
    per_target = (
        fit.groupby(["transformation", "component", "target"], as_index=False)
        .agg(
            delta_residual=("delta_residual", "mean"),
            delta_affinity=("delta_affinity", "mean"),
            observations=("pair_id", "size"),
        )
    )
    return (
        per_target.groupby("transformation", as_index=False)
        .agg(
            fit_components=("component", "nunique"),
            fit_targets=("target", "nunique"),
            fit_target_units=("target", "size"),
            mean_residual_effect=("delta_residual", "mean"),
            mean_affinity_effect=("delta_affinity", "mean"),
            sd_affinity_effect=("delta_affinity", "std"),
        )
        .sort_values(["fit_components", "fit_targets"], ascending=False)
        .reset_index(drop=True)
    )


def softplus_loss(delta_affinity: np.ndarray, prediction: np.ndarray) -> np.ndarray:
    direction = np.sign(delta_affinity)
    return np.logaddexp(0.0, -direction * prediction)


def transfer_summary(edges: pd.DataFrame, priors: pd.DataFrame) -> dict[str, object]:
    robust = priors.loc[
        (priors.fit_components >= MIN_PRIOR_COMPONENTS)
        & (priors.fit_targets >= MIN_PRIOR_COMPONENTS)
    ]
    probe = edges.loc[edges.role == "probe"].merge(
        robust[["transformation", "mean_residual_effect"]],
        on="transformation",
        how="inner",
        validate="many_to_one",
    )
    output: dict[str, object] = {
        "fit_transformations": int(len(priors)),
        "robust_fit_transformations": int(len(robust)),
        "probe_edges_with_robust_prior": int(len(probe)),
        "probe_targets_with_robust_prior": int(probe.target.nunique()),
        "probe_components_with_robust_prior": int(probe.component.nunique()),
        "strata": {},
    }
    if probe.empty:
        return output
    probe = probe.copy()
    probe["prior_prediction"] = probe.delta_base + probe.mean_residual_effect
    probe["proper_gain"] = softplus_loss(
        probe.delta_affinity.to_numpy(), probe.delta_base.to_numpy()
    ) - softplus_loss(
        probe.delta_affinity.to_numpy(), probe.prior_prediction.to_numpy()
    )
    active = probe.delta_affinity.abs() >= MIN_DIRECTION_EFFECT
    probe["direction_gain"] = np.nan
    probe.loc[active, "direction_gain"] = (
        np.sign(probe.loc[active, "prior_prediction"])
        == np.sign(probe.loc[active, "delta_affinity"])
    ).astype(float) - (
        np.sign(probe.loc[active, "delta_base"])
        == np.sign(probe.loc[active, "delta_affinity"])
    ).astype(float)
    masks = {
        "all": np.ones(len(probe), dtype=bool),
        "low_similarity": probe.tanimoto.to_numpy() < LOW_SIMILARITY,
        "local_similarity": probe.tanimoto.to_numpy() >= LOW_SIMILARITY,
    }
    for name, mask in masks.items():
        frame = probe.loc[mask]
        output["strata"][name] = {
            "edges": int(len(frame)),
            "targets": int(frame.target.nunique()),
            "components": int(frame.component.nunique()),
            "proper_gain": paired_bootstrap(frame, "proper_gain", draws=BOOTSTRAP_DRAWS),
            "direction_gain": paired_bootstrap(
                frame, "direction_gain", draws=BOOTSTRAP_DRAWS
            ),
        }
    return output


def conn_fragment_map(fragments: pd.DataFrame) -> dict[str, dict[str, str]]:
    output: dict[str, dict[str, str]] = {}
    for row in fragments.itertuples(index=False):
        output.setdefault(str(row.conn), {})[str(row.core)] = str(row.substituent)
    return output


def nearest_tanimoto_bits(
    bits: np.ndarray,
    bit_counts: np.ndarray,
    query: np.ndarray,
    support: np.ndarray,
) -> np.ndarray:
    query_bits = bits[query].astype(np.int32, copy=False)
    support_bits = bits[support].astype(np.int32, copy=False)
    intersection = query_bits @ support_bits.T
    union = bit_counts[query, None] + bit_counts[support][None, :] - intersection
    return np.max(intersection / np.maximum(union, 1), axis=1)


def episode_coverage(
    substrate: Substrate,
    fragments: pd.DataFrame,
    known_transformations: set[str],
) -> pd.DataFrame:
    fragment_map = conn_fragment_map(fragments)
    conns = substrate.labeled.conn.astype(str).to_numpy()
    bits = substrate.bits.detach().cpu().numpy().astype(np.uint8)
    bit_counts = bits.sum(axis=1, dtype=np.int32)
    records: list[dict[str, object]] = []
    for episode in build_episodes(substrate, "probe"):
        query = np.asarray(episode.query, dtype=np.int64)
        for k in K_VALUES:
            support = np.asarray(episode.support5[:k], dtype=np.int64)
            nearest = nearest_tanimoto_bits(bits, bit_counts, query, support)
            support_cores: dict[str, set[str]] = {}
            for row_index in support:
                for core, substituent in fragment_map.get(conns[row_index], {}).items():
                    support_cores.setdefault(core, set()).add(substituent)
            for query_index, similarity in zip(query, nearest):
                connected = False
                known = False
                for core, query_substituent in fragment_map.get(conns[query_index], {}).items():
                    for support_substituent in support_cores.get(core, ()):
                        if query_substituent == support_substituent:
                            continue
                        connected = True
                        if transform_key(query_substituent, support_substituent) in known_transformations:
                            known = True
                records.append(
                    {
                        "episode_id": episode.episode_id,
                        "target": episode.target,
                        "component": episode.component,
                        "draw": episode.draw,
                        "k": k,
                        "query_row": int(query_index),
                        "nearest_tanimoto": float(similarity),
                        "mmp_connected": connected,
                        "fit_known_transformation": known,
                    }
                )
    return pd.DataFrame.from_records(records)


def coverage_summary(coverage: pd.DataFrame) -> dict[str, object]:
    output: dict[str, object] = {}
    for k in K_VALUES:
        k_frame = coverage.loc[coverage.k == k]
        output[f"k{k}"] = {}
        for name, frame in {
            "all": k_frame,
            "low_similarity": k_frame.loc[k_frame.nearest_tanimoto < LOW_SIMILARITY],
        }.items():
            episode_counts = (
                frame.groupby(["episode_id", "component", "target"], as_index=False)
                .agg(
                    queries=("query_row", "size"),
                    connected=("mmp_connected", "sum"),
                    known=("fit_known_transformation", "sum"),
                )
            )
            powered = episode_counts.loc[episode_counts.known >= MIN_QUERY]
            output[f"k{k}"][name] = {
                "queries": int(len(frame)),
                "targets": int(frame.target.nunique()),
                "components": int(frame.component.nunique()),
                "mmp_connected_fraction": float(frame.mmp_connected.mean()) if len(frame) else 0.0,
                "fit_known_fraction": float(frame.fit_known_transformation.mean()) if len(frame) else 0.0,
                "episodes_with_at_least_8_known_queries": int(len(powered)),
                "targets_with_at_least_8_known_queries": int(powered.target.nunique()),
                "components_with_at_least_8_known_queries": int(powered.component.nunique()),
            }
    return output


def decide(coverage: dict[str, object], transfer: dict[str, object]) -> dict[str, object]:
    k3_components = coverage.get("k3", {}).get("low_similarity", {}).get(
        "components_with_at_least_8_known_queries", 0
    )
    k5_components = coverage.get("k5", {}).get("low_similarity", {}).get(
        "components_with_at_least_8_known_queries", 0
    )
    low_transfer = transfer.get("strata", {}).get("low_similarity", {})
    proper_lower = low_transfer.get("proper_gain", {}).get("lower95", float("nan"))
    admitted = bool(
        k3_components >= POWERED_COMPONENTS
        and k5_components >= POWERED_COMPONENTS
        and np.isfinite(proper_lower)
        and proper_lower > 0.0
    )
    return {
        "k3_powered_components": int(k3_components),
        "k5_powered_components": int(k5_components),
        "required_components": POWERED_COMPONENTS,
        "low_similarity_prior_proper_gain_lower95": float(proper_lower),
        "verdict": (
            "EXPLICIT_TRANSFORMATION_GRAMMAR_ADMITTED"
            if admitted
            else "EXPLICIT_TRANSFORMATION_GRAMMAR_NOT_ADMITTED"
        ),
        "tcrs_relation_encoder_status": "REQUIRES_R1_LEARNED_RELATION_GATE",
    }


def run(
    lock_path: Path,
    oof_cache: Path,
    output_path: Path,
    fragments_path: Path,
    edges_path: Path,
    coverage_path: Path,
) -> dict[str, object]:
    substrate, context = load_substrate(lock_path, oof_cache)
    if set(substrate.labeled.role.unique()) - {"fit", "probe"}:
        raise AssertionError("an unauthorized role entered the relation census")
    measurements = aggregate_measurements(substrate)
    fragments = build_fragment_frame(measurements.conn)
    edges = build_mmp_edges(measurements, fragments)
    priors = fit_transformation_priors(edges)
    robust = priors.loc[
        (priors.fit_components >= MIN_PRIOR_COMPONENTS)
        & (priors.fit_targets >= MIN_PRIOR_COMPONENTS)
    ]
    coverage = episode_coverage(substrate, fragments, set(robust.transformation))
    transfer = transfer_summary(edges, priors)
    coverage_stats = coverage_summary(coverage)
    decision = decide(coverage_stats, transfer)

    for path in (fragments_path, edges_path, coverage_path, output_path):
        path.parent.mkdir(parents=True, exist_ok=True)
    fragments.to_parquet(fragments_path, index=False)
    edges.to_parquet(edges_path, index=False)
    coverage.to_parquet(coverage_path, index=False)
    result = {
        "schema": "a2s-thermodynamic-relation-census-v1",
        "status": "SOURCE_ONLY",
        "protocol": {
            "fragmentation": "single acyclic MMP cut; larger fragment is core",
            "min_core_heavy": MIN_CORE_HEAVY,
            "min_core_fraction": MIN_CORE_FRACTION,
            "max_substituent_heavy": MAX_SUBSTITUENT_HEAVY,
            "robust_prior_min_fit_components": MIN_PRIOR_COMPONENTS,
            "low_similarity_threshold": LOW_SIMILARITY,
            "support_policy": "nested random within exact assay",
            "k_values": list(K_VALUES),
        },
        "data": {
            "roles_opened": ["fit", "probe"],
            "locked_labels_requested": False,
            "recipient_labels_requested": False,
            "measurements": int(len(measurements)),
            "unique_molecules": int(measurements.conn.nunique()),
            "fragment_records": int(len(fragments)),
            "fragmented_molecules": int(fragments.conn.nunique()),
            "mmp_edges": int(len(edges)),
            "fit_edges": int((edges.role == "fit").sum()),
            "probe_edges": int((edges.role == "probe").sum()),
            "source_context": context,
        },
        "transfer": transfer,
        "episode_coverage": coverage_stats,
        "decision": decision,
        "artifacts": {
            "fragments": str(fragments_path.relative_to(ROOT)).replace("\\", "/"),
            "edges": str(edges_path.relative_to(ROOT)).replace("\\", "/"),
            "coverage": str(coverage_path.relative_to(ROOT)).replace("\\", "/"),
        },
    }
    result["content_sha256"] = sha256(canonical(result).encode()).hexdigest()
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--oof-cache", type=Path, default=DEFAULT_OOF)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--fragments", type=Path, default=DEFAULT_FRAGMENTS)
    parser.add_argument("--edges", type=Path, default=DEFAULT_EDGES)
    parser.add_argument("--coverage", type=Path, default=DEFAULT_COVERAGE)
    args = parser.parse_args()
    result = run(
        args.lock,
        args.oof_cache,
        args.output,
        args.fragments,
        args.edges,
        args.coverage,
    )
    print(json.dumps(result["decision"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

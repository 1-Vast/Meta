"""TRAIN-only protein-ligand interaction identifiability audit.

The audit is deliberately model-free.  It separates endpoint strata, builds
exact-ligand 2x2 rectangles within a source document, estimates replicate
noise from repeated raw measurements, and reports component-level summaries.
No development, confirmation, or sealed labels are loaded.
"""

from __future__ import annotations

import argparse
import gzip
import itertools
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "dataset/public/chembl_37/processed"
REGISTRY = SOURCE / "dualcold/registry.parquet"
RDLogger.DisableLog("rdApp.*")


def _json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def canonical_ligand(smiles: str) -> str | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def connectivity_ligand(smiles: str) -> str | None:
    """Return the registry's non-isomeric ligand-parent representation."""

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, isomericSmiles=False)


def _registry_train_keys() -> tuple[set[tuple[str, str, str]], pd.DataFrame]:
    """Return registry-closed target/ligand/endpoint keys and train metadata."""

    columns = [
        "target",
        "conn",
        "endpoint",
        "dual_cold_split",
        "hcluster",
        "replicate_sd",
    ]
    registry = pd.read_parquet(REGISTRY, columns=columns)
    train = registry.loc[registry.dual_cold_split == "train"].copy()
    keys = {
        (str(row.target), str(row.conn), str(row.endpoint))
        for row in train.itertuples(index=False)
    }
    return keys, train


def load_train_rows(
    *, registry_closed: bool = False, range_policy: str = "none"
) -> tuple[pd.DataFrame, dict[str, object]]:
    keys, registry_train = _registry_train_keys()
    split = pd.read_parquet(REGISTRY, columns=["target", "dual_cold_split", "hcluster"])
    train_split = split.loc[split.dual_cold_split == "train"]
    train_targets = set(train_split["target"])
    cluster = train_split.drop_duplicates("target").set_index("target")["hcluster"].to_dict()
    if range_policy not in {"none", "clip"}:
        raise ValueError("range_policy must be 'none' or 'clip'")
    rows: list[tuple] = []
    skipped = 0
    source_counts: Counter[str] = Counter()
    registry_counts: Counter[str] = Counter()
    excluded_registry_counts: Counter[str] = Counter()
    raw_range_counts: Counter[str] = Counter()
    for filename in ("chembl37_pKi.jsonl.gz", "chembl37_pKd.jsonl.gz"):
        with gzip.open(SOURCE / filename, "rt", encoding="utf-8") as handle:
            for line in handle:
                record = json.loads(line)
                target = record["target_chembl_id"]
                if target not in train_targets:
                    continue
                ligand = canonical_ligand(record["smiles"])
                if ligand is None:
                    skipped += 1
                    continue
                endpoint = str(record["endpoint"])
                connectivity = connectivity_ligand(record["smiles"])
                registry_member = (target, connectivity, endpoint) in keys
                if registry_closed and not registry_member:
                    excluded_registry_counts[endpoint] += 1
                    continue
                source_counts[endpoint] += 1
                registry_counts[endpoint] += int(registry_member)
                value = float(record["pK"])
                if not 0.0 <= value <= 14.0:
                    raw_range_counts[endpoint] += 1
                if range_policy == "clip":
                    value = float(np.clip(value, 0.0, 14.0))
                rows.append(
                    (
                        endpoint,
                        str(record["doc_id"]),
                        str(record["assay_id"]),
                        target,
                        ligand,
                        value,
                    )
                )
    frame = pd.DataFrame(rows, columns=["endpoint", "doc", "assay", "target", "ligand", "y"])
    audit = {
        "raw_train_rows": len(frame),
        "skipped_invalid_smiles": skipped,
        "train_targets": len(train_targets),
        "train_homology_components": len(set(cluster.values())),
        "endpoint_rows": dict(source_counts),
        "registry_train_rows": {
            str(endpoint): int(count)
            for endpoint, count in registry_train.groupby("endpoint").size().items()
        },
        "registry_member_rows": dict(registry_counts),
        "excluded_nonregistry_rows": dict(excluded_registry_counts),
        "registry_closed": registry_closed,
        "range_policy": range_policy,
        "raw_pK_out_of_range_rows": dict(raw_range_counts),
        "all_targets_train": bool(set(frame.target) <= train_targets),
        "development_labels_read": False,
        "confirmation_labels_read": False,
        "sealed_labels_read": False,
    }
    return frame, audit


def aggregate_cells(rows: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return cell medians and within-cell replicate differences."""
    keys = ["endpoint", "doc", "assay", "target", "ligand"]
    cells = rows.groupby(keys, as_index=False).agg(y=("y", "median"), n=("y", "size"))
    replicate_records: list[dict[str, object]] = []
    for key, group in rows.groupby(keys, sort=False):
        values = group.y.to_numpy(dtype=float)
        if len(values) < 2:
            continue
        # One variance estimate per source cell; it is not treated as an
        # independent rectangle and is only used to calibrate propagated noise.
        diffs = values[:, None] - values[None, :]
        upper = diffs[np.triu_indices(len(values), 1)]
        replicate_records.append({
            **dict(zip(keys, key)),
            "n": len(values),
            # For a two-replicate cell, |x1-x2|/sqrt(2) is the RMS error
            # estimate.  RMS pair differences remain defined for all n.
            "pair_sd": float(np.sqrt(np.mean(upper ** 2) / 2.0)),
        })
    return cells, pd.DataFrame(replicate_records)


def _component_bootstrap(values: pd.DataFrame, component_col: str, seed: int, draws: int = 2000) -> dict[str, float | None]:
    if values.empty:
        return {"mean": None, "q025": None, "q975": None, "components": 0}
    grouped = values.groupby(component_col).dd.agg(["mean", "count"]).reset_index()
    if grouped.empty:
        return {"mean": None, "q025": None, "q975": None, "components": 0}
    rng = np.random.default_rng(seed)
    weights = grouped["count"].to_numpy(dtype=float)
    means = grouped["mean"].to_numpy(dtype=float)
    sampled = rng.integers(0, len(means), size=(draws, len(means)))
    boot = (means[sampled] * weights[sampled]).sum(axis=1) / weights[sampled].sum(axis=1)
    return {
        "mean": float(np.average(means, weights=weights)),
        "q025": float(np.quantile(boot, 0.025)),
        "q975": float(np.quantile(boot, 0.975)),
        "components": int(len(means)),
    }


def _registry_replicates() -> pd.DataFrame:
    """Use registry-level replicate SD as a conservative noise sensitivity."""

    _, train = _registry_train_keys()
    values = train.loc[train.replicate_sd.notna(), ["endpoint", "replicate_sd"]]
    return values.rename(columns={"replicate_sd": "pair_sd"}).reset_index(drop=True)


def build_rectangles(
    cells: pd.DataFrame,
    target_component: dict[str, str],
    *,
    cross_component_only: bool = True,
) -> pd.DataFrame:
    """Build document-local target pairs and exact common-ligand pairs."""
    panel_rows: list[tuple] = []
    # A panel is one endpoint/document/target/assay.  Target pairs are only
    # compared when both panels are in the same document and share >=2 ligands.
    for (endpoint, doc), document in cells.groupby(["endpoint", "doc"], sort=True):
        panels = []
        for (panel_endpoint, panel_doc, target, assay), group in document.groupby(
            ["endpoint", "doc", "target", "assay"], sort=False
        ):
            values = group.set_index("ligand").y.to_dict()
            if len(values) >= 2:
                panels.append((target, assay, values))
        for (target_a, assay_a, values_a), (target_b, assay_b, values_b) in itertools.combinations(panels, 2):
            if target_a == target_b:
                continue
            if cross_component_only and target_component.get(target_a) == target_component.get(target_b):
                continue
            common = sorted(set(values_a) & set(values_b))
            if len(common) < 2:
                continue
            if target_b < target_a:
                target_a, target_b = target_b, target_a
                assay_a, assay_b = assay_b, assay_a
                values_a, values_b = values_b, values_a
            for ligand_i, ligand_j in itertools.combinations(common, 2):
                dd = (values_a[ligand_i] - values_a[ligand_j]) - (values_b[ligand_i] - values_b[ligand_j])
                sign_a = values_a[ligand_i] - values_a[ligand_j]
                sign_b = values_b[ligand_i] - values_b[ligand_j]
                target_pair = f"{target_a}|{target_b}"
                homology_pair = "|".join(sorted((target_component.get(target_a, ""), target_component.get(target_b, ""))))
                panel_rows.append(
                    (endpoint, doc, target_a, target_b, assay_a, assay_b, ligand_i, ligand_j,
                     target_pair, homology_pair, float(dd), int(sign_a * sign_b < 0),
                     int(assay_a == assay_b), abs(float(dd)))
                )
    return pd.DataFrame(
        panel_rows,
        columns=["endpoint", "doc", "target_a", "target_b", "assay_a", "assay_b", "ligand_i", "ligand_j",
                 "target_pair", "homology_pair", "dd", "order_reversal", "same_assay", "abs_dd"],
    )


def _summarize_flat(rectangles: pd.DataFrame, replicates: pd.DataFrame, seed: int) -> dict[str, object]:
    if rectangles.empty:
        return {"rectangles": 0, "target_pairs": 0, "homology_pairs": 0}
    # Independent unit is target-pair x document, not a four-row rectangle.
    grouped_columns = ["endpoint", "target_pair", "doc"]
    unit = rectangles.groupby(grouped_columns, as_index=False).agg(
        dd=("dd", "median"), order_reversal=("order_reversal", "mean"), n=("dd", "size"),
        homology_pair=("homology_pair", "first"), same_assay=("same_assay", "max"),
    )
    panel_counts = rectangles.assign(
        assay_pair=rectangles.assay_a.astype(str) + "|" + rectangles.assay_b.astype(str)
    ).groupby(["endpoint", "target_pair", "doc"], sort=False).assay_pair.nunique()
    # Approximate four-cell propagated noise from available replicate SDs.  The
    # denominator is reported explicitly; it is not used to inflate sample size.
    rep_sd = float(replicates.pair_sd.median()) if not replicates.empty else None
    noise_sd = math.sqrt(4.0) * rep_sd if rep_sd is not None else None
    noise_q90 = math.sqrt(4.0) * float(replicates.pair_sd.quantile(0.90)) if not replicates.empty else None
    noise_q95 = math.sqrt(4.0) * float(replicates.pair_sd.quantile(0.95)) if not replicates.empty else None
    summary: dict[str, object] = {
        "rectangles": int(len(rectangles)),
        "target_pairs": int(rectangles.target_pair.nunique()),
        "ligand_pairs": int(rectangles[["ligand_i", "ligand_j"]].drop_duplicates().shape[0]),
        "homology_pairs": int(rectangles.homology_pair.nunique()),
        "documents": int(rectangles.doc.nunique()),
        "same_assay_rectangles": int(rectangles.same_assay.sum()),
        "median_abs_dd": float(rectangles.abs_dd.median()),
        "mean_abs_dd": float(rectangles.abs_dd.mean()),
        "dd_sd": float(rectangles.dd.std(ddof=1)) if len(rectangles) > 1 else 0.0,
        "order_reversal_fraction": float(rectangles.order_reversal.mean()),
        "replicate_cells": int(len(replicates)),
        "replicate_pair_sd_median": rep_sd,
        "replicate_pair_sd_quantiles": (
            {str(q): float(replicates.pair_sd.quantile(q)) for q in (0.5, 0.9, 0.95, 0.99)}
            if not replicates.empty else None
        ),
        "replicate_nonzero_fraction": (
            float((replicates.pair_sd > 0).mean()) if not replicates.empty else None
        ),
        "propagated_four_cell_noise_sd": noise_sd,
        "propagated_four_cell_noise_q90": noise_q90,
        "propagated_four_cell_noise_q95": noise_q95,
        "unit_count": int(len(unit)),
        "panel_pair_count": int(panel_counts.sum()),
        "multi_panel_unit_count": int((panel_counts > 1).sum()),
        "multi_panel_unit_fraction": float((panel_counts > 1).mean()),
        "multi_panel_pair_fraction": float(
            panel_counts.loc[panel_counts > 1].sum() / panel_counts.sum()
        ),
        "unit_median_abs_dd": float(unit.dd.abs().median()),
        "unit_order_reversal_fraction": float(unit.order_reversal.mean()),
        "unit_component_bootstrap_dd": _component_bootstrap(unit, "homology_pair", seed),
        "unit_component_bootstrap_abs_dd": _component_bootstrap(unit.assign(dd=unit.dd.abs()), "homology_pair", seed + 1),
        "unit_component_bootstrap_reversal": _component_bootstrap(unit.assign(dd=unit.order_reversal), "homology_pair", seed + 2),
        "largest_target_pair_unit_fraction": float(unit.target_pair.value_counts(normalize=True).iloc[0]),
        "largest_homology_pair_unit_fraction": float(unit.homology_pair.value_counts(normalize=True).iloc[0]),
    }
    if noise_sd is not None:
        summary["dd_to_noise_ratio"] = float(summary["dd_sd"]) / noise_sd if noise_sd > 0 else None
        summary["dd_to_noise_ratio_q90"] = float(summary["dd_sd"]) / noise_q90 if noise_q90 and noise_q90 > 0 else None
        summary["dd_to_noise_ratio_q95"] = float(summary["dd_sd"]) / noise_q95 if noise_q95 and noise_q95 > 0 else None
    if "document_family" in rectangles.columns:
        family_counts = rectangles.groupby("document_family").doc.nunique()
        summary["document_families"] = int(len(family_counts))
        summary["largest_document_family_unit_fraction"] = float(
            unit.assign(
                document_family=unit.doc.map(
                    rectangles.drop_duplicates("doc").set_index("doc")["document_family"]
                ).fillna(unit.doc)
            ).document_family.value_counts(normalize=True).iloc[0]
        )
    return summary


def summarize(rectangles: pd.DataFrame, replicates: pd.DataFrame, seed: int) -> dict[str, object]:
    summary = _summarize_flat(rectangles, replicates, seed)
    summary["by_endpoint"] = {
        endpoint: _summarize_flat(
            group,
            replicates.loc[replicates.endpoint == endpoint],
            seed + (17 if endpoint == "pKd" else 31),
        )
        for endpoint, group in rectangles.groupby("endpoint")
    }
    return summary


def document_source_families() -> dict[str, str]:
    """Map documents to the coarsest source family available locally."""

    path = SOURCE / "dualcold/pcic_o0_document_metadata.json"
    if not path.is_file():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    mapping: dict[str, str] = {}
    for record in payload.get("documents", []):
        document = record.get("document_chembl_id")
        if document is None:
            continue
        # ChEMBL's src_id is a source/provider identifier, not an assay
        # comparability guarantee.  It is used only for dependence sensitivity.
        mapping[str(document)] = f"src:{record.get('src_id', 'unknown')}"
    return mapping


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "reports/active/interaction_identifiability_audit_2026-07-31.json")
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument(
        "--registry-closed",
        action="store_true",
        help="exclude raw target/ligand rows removed by the active registry collision policy",
    )
    parser.add_argument(
        "--range-policy",
        choices=("none", "clip"),
        default="none",
        help="optional TRAIN-only pK sensitivity; clip raw values to [0, 14]",
    )
    parser.add_argument(
        "--noise-source",
        choices=("cell", "registry"),
        default="cell",
        help="within-document replicate noise or conservative registry replicate_sd",
    )
    parser.add_argument(
        "--allow-same-component",
        action="store_true",
        help="retain same-homology target pairs as a non-strict sensitivity",
    )
    args = parser.parse_args()
    rows, source_audit = load_train_rows(
        registry_closed=args.registry_closed, range_policy=args.range_policy
    )
    cells, replicates = aggregate_cells(rows)
    if args.noise_source == "registry":
        replicates = _registry_replicates()
    split = pd.read_parquet(REGISTRY, columns=["target", "hcluster", "dual_cold_split"])
    components = split.loc[split.dual_cold_split == "train"].drop_duplicates("target").set_index("target")["hcluster"].to_dict()
    rectangles = build_rectangles(
        cells, components, cross_component_only=not args.allow_same_component
    )
    family_map = document_source_families()
    if not rectangles.empty:
        rectangles["document_family"] = rectangles.doc.map(family_map).fillna(
            rectangles.doc.map(lambda value: f"doc:{value}")
        )
    summary = summarize(rectangles, replicates, args.seed)
    result = {
        "schema": "protein-ligand-interaction-identifiability-audit-v1",
        "estimand": "endpoint-specific document-local exact-ligand 2x2 difference-in-differences",
        "source_audit": source_audit,
        "noise_source": args.noise_source,
        "cross_component_only": not args.allow_same_component,
        "range_policy": args.range_policy,
        "cell_audit": {"cells": int(len(cells)), "replicate_cells": int(len(replicates))},
        "summary": summary,
        "method_notes": [
            "pKi and pKd are never pooled.",
            "A rectangle is not treated as an independent observation; target-pair/document units and homology-pair bootstrap are reported.",
            (
                "Noise uses repeated raw within-cell measurements and is propagated over four cells."
                if args.noise_source == "cell"
                else "Noise uses registry target-ligand replicate_sd as a conservative sensitivity and is propagated over four cells."
            ),
            "Same-assay counts are reported separately; document-local rectangles may use distinct assays for the two targets.",
            "Registry-closed and [0,14]-clipped runs are sensitivity analyses, not additional independent evidence.",
            "Source-family labels are dependence blocks only; src_id does not establish assay comparability.",
            "Strict runs exclude target pairs sharing one registered homology component; --allow-same-component is sensitivity only.",
        ],
        "admission": {
            "exact_ligand_rectangles_exist": bool(len(rectangles)),
            "strict_same_assay_rectangles_exist": bool(len(rectangles) and rectangles.same_assay.sum() > 0),
            "cross_homology_units_exist": bool(len(rectangles) and rectangles.homology_pair.nunique() > 1),
            "component_concentration_below_25_percent": bool(
                len(rectangles) and summary["largest_homology_pair_unit_fraction"] < 0.25
            ),
        },
        "verdict": (
            "AUDIT_RECTANGLES_EXIST_STRICT_ASSAY_COMPARABILITY_LIMITED"
            if len(rectangles) and rectangles.same_assay.sum() == 0
            else ("AUDIT_DATA_SUPPORTS_FURTHER_PROBE" if len(rectangles) else "AUDIT_NO_RECTANGLES_STOP")
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, default=_json_default) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, default=_json_default))


if __name__ == "__main__":
    main()

"""Component-held-out probes for measured same-ligand cross-protein selectivity."""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import numpy as np
from sklearn.linear_model import Ridge

from research.e0_identifiability.run_tdir_pilot import _load_protein_rows
from research.meta_fewshot.a1_selectivity_dependency_audit import OUT as DEPENDENCY_RESULT
from research.meta_fewshot.seal_v1_development import read_gzip_jsonl, sha256

ROOT = Path(__file__).resolve().parents[2]
MAIN = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0"
DEV = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_v1_development"
PROTEIN_BANK = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_main_v0_protein_bank"
OUT = ROOT / "report/meta_fewshot/a1_selectivity_probe.json"
ALPHAS = (0.01, 0.1, 1.0, 10.0, 100.0)
AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"


def protein_nuisance(sequence: str) -> np.ndarray:
    counts = np.asarray([sequence.count(code) for code in AMINO_ACIDS], dtype=np.float64)
    return np.concatenate(([np.log1p(len(sequence))], counts / max(len(sequence), 1)))


def stable_rng(value: str) -> np.random.Generator:
    seed = int(hashlib.sha256(value.encode()).hexdigest()[:16], 16)
    return np.random.default_rng(seed)


def build_rows():
    dependency = json.loads(DEPENDENCY_RESULT.read_text())
    if not dependency["probe_identifiable"]:
        raise RuntimeError("selectivity dependency Gate is closed")
    groups = [row for row in read_gzip_jsonl(DEV / "source_contrast_groups.jsonl.gz")
              if row["kind"] == "measured_partner"]
    source_cells = list(read_gzip_jsonl(DEV / "source_cells.jsonl.gz"))
    source_position = {row["cell_id"]: index for index, row in enumerate(source_cells)}
    with np.load(DEV / "source_features.npz", allow_pickle=False) as stored:
        if stored["cell_id"].tolist() != [row["cell_id"] for row in source_cells]:
            raise ValueError("source feature order mismatch")
        tbasis = stored["correct"].astype(np.float64)
    proteins_json = [json.loads(line) for line in (MAIN / "proteins.jsonl").read_text().splitlines()]
    sequence = {row["sequence_sha256"]: row["sequence"] for row in proteins_json}
    targets = {member["target_id"] for group in groups for member in group["members"]}
    protein_bank = _load_protein_rows(PROTEIN_BANK, targets)
    pooled = {target: protein_bank[target]["pooled"].astype(np.float64) for target in targets}

    rows = []
    for group_index, group in enumerate(groups):
        group_id = f"{group_index:06d}"
        by_family = defaultdict(list)
        for member in group["members"]:
            by_family[member["protein_group_40"]].append(member)
        family_rows = []
        for family, members in sorted(by_family.items()):
            family_rows.append({
                "family": family,
                "y": float(np.median([member["pK"] for member in members])),
                "tbasis": np.median(np.stack([
                    tbasis[source_position[member["cell_id"]]] for member in members]), axis=0),
                "esm": np.median(np.stack([pooled[member["target_id"]] for member in members]), axis=0),
                "nuisance": np.median(np.stack([
                    protein_nuisance(sequence[member["target_id"]]) for member in members]), axis=0),
            })
        if len(family_rows) < 2:
            raise ValueError("measured partner group lost cross-family depth")
        y_mean = np.mean([row["y"] for row in family_rows])
        means = {name: np.mean(np.stack([row[name] for row in family_rows]), axis=0)
                 for name in ("tbasis", "esm", "nuisance")}
        permutation = stable_rng(f"coupling-null|{group_id}").permutation(len(family_rows))
        for index, row in enumerate(family_rows):
            rows.append({
                "group": group_id, "component": dependency["group_component_assignment"][group_id],
                "family": row["family"], "y": row["y"] - y_mean,
                "tbasis": row["tbasis"] - means["tbasis"],
                "tbasis_null": family_rows[int(permutation[index])]["tbasis"] - means["tbasis"],
                "esm": row["esm"] - means["esm"],
                "nuisance": row["nuisance"] - means["nuisance"],
            })
    return rows


def component_folds(rows, folds=5):
    counts = defaultdict(set)
    for row in rows:
        counts[row["component"]].add(row["group"])
    allocation = [[] for _ in range(folds)]
    loads = [0] * folds
    for component, groups in sorted(counts.items(), key=lambda item: (-len(item[1]), item[0])):
        fold = min(range(folds), key=lambda index: (loads[index], index))
        allocation[fold].append(component)
        loads[fold] += len(groups)
    return allocation


def standardize(train, test):
    mean, scale = train.mean(axis=0), train.std(axis=0)
    scale[scale < 1e-8] = 1.0
    return (train - mean) / scale, (test - mean) / scale


def component_losses(prediction, y, rows, indices):
    group_loss = defaultdict(list)
    group_component = {}
    for value, truth, index in zip(prediction, y, indices):
        row = rows[int(index)]
        group_loss[row["group"]].append(float((value - truth) ** 2))
        group_component[row["group"]] = row["component"]
    by_component = defaultdict(list)
    for group, values in group_loss.items():
        by_component[group_component[group]].append(float(np.mean(values)))
    return {component: float(np.mean(values)) for component, values in by_component.items()}


def contrast_sensitivity(left, right, component_group_counts):
    components = sorted(set(left) & set(right))
    delta = {key: left[key] - right[key] for key in components}
    giant = max(components, key=lambda key: component_group_counts[key])
    leave_one_out = [
        float(np.mean([delta[key] for key in components if key != omitted]))
        for omitted in components
    ]
    return {
        "components_favoring_correct": int(sum(value > 0 for value in delta.values())),
        "components_favoring_control": int(sum(value < 0 for value in delta.values())),
        "group_weighted_mean": float(sum(
            delta[key] * component_group_counts[key] for key in components
        ) / sum(component_group_counts[key] for key in components)),
        "giant_component": giant,
        "giant_component_delta": float(delta[giant]),
        "leave_giant_out_component_macro_mean": float(np.mean([
            delta[key] for key in components if key != giant
        ])),
        "leave_one_component_out_min": float(min(leave_one_out)),
        "leave_one_component_out_max": float(max(leave_one_out)),
    }


def fold_indices(rows, heldout_components):
    heldout = set(heldout_components)
    train = np.asarray([index for index, row in enumerate(rows) if row["component"] not in heldout])
    test = np.asarray([index for index, row in enumerate(rows) if row["component"] in heldout])
    return train, test


def select_alpha(X, y, rows, train_indices):
    train_components = sorted({rows[int(index)]["component"] for index in train_indices})
    inner = [train_components[index::3] for index in range(3)]
    values = []
    for alpha in ALPHAS:
        losses = []
        for heldout in inner:
            inner_train = np.asarray([index for index in train_indices
                                      if rows[int(index)]["component"] not in set(heldout)])
            inner_test = np.asarray([index for index in train_indices
                                     if rows[int(index)]["component"] in set(heldout)])
            train_x, test_x = standardize(X[inner_train], X[inner_test])
            model = Ridge(alpha=alpha, fit_intercept=False).fit(train_x, y[inner_train])
            losses.extend(component_losses(model.predict(test_x), y[inner_test], rows,
                                            inner_test).values())
        values.append((float(np.mean(losses)), alpha))
    return min(values)[1]


def cross_validate(X, y, rows, folds):
    prediction = np.empty(len(rows), dtype=np.float64)
    operators, selected = [], []
    for heldout in folds:
        train_index, test_index = fold_indices(rows, heldout)
        alpha = select_alpha(X, y, rows, train_index)
        train_x, test_x = standardize(X[train_index], X[test_index])
        covariance = train_x.T @ train_x + alpha * np.eye(train_x.shape[1])
        weights = np.linalg.solve(covariance, train_x.T)
        operator = test_x @ weights
        prediction[test_index] = operator @ y[train_index]
        operators.append((train_index, test_index, operator))
        selected.append(alpha)
    return prediction, operators, selected


def score(prediction, y, rows):
    indices = np.arange(len(rows))
    losses = component_losses(prediction, y, rows, indices)
    group_values = defaultdict(list)
    for index, row in enumerate(rows):
        group_values[row["group"]].append(index)
    correlations, signs = defaultdict(list), defaultdict(list)
    comparable_pairs = truth_tied_pairs = prediction_tied_pairs = 0
    for group, group_indices in group_values.items():
        truth, pred = y[group_indices], prediction[group_indices]
        component = rows[group_indices[0]]["component"]
        if np.std(truth) > 0 and np.std(pred) > 0:
            correlations[component].append(float(np.corrcoef(truth, pred)[0, 1]))
        pair_sign = []
        for left, right in combinations(range(len(group_indices)), 2):
            truth_sign = np.sign(truth[left] - truth[right])
            if truth_sign == 0:
                truth_tied_pairs += 1
                continue
            comparable_pairs += 1
            prediction_sign = np.sign(pred[left] - pred[right])
            prediction_tied_pairs += int(prediction_sign == 0)
            pair_sign.append(0.5 if prediction_sign == 0 else float(prediction_sign == truth_sign))
        if pair_sign:
            signs[component].append(float(np.mean(pair_sign)))
    zero = component_losses(np.zeros_like(y), y, rows, indices)
    return {
        "component_macro_mse": float(np.mean(list(losses.values()))),
        "component_macro_r2_vs_zero": float(np.mean([
            1.0 - losses[key] / zero[key] for key in losses if zero[key] > 0])),
        "component_macro_group_pearson": (
            float(np.mean([np.mean(values) for values in correlations.values() if values]))
            if any(correlations.values()) else None),
        "component_macro_pair_sign_accuracy": (
            float(np.mean([np.mean(values) for values in signs.values() if values]))
            if any(signs.values()) else None),
        "pair_sign_counts": {
            "comparable_truth_pairs": comparable_pairs,
            "truth_tied_pairs_exact": truth_tied_pairs,
            "prediction_tied_pairs_among_comparable": prediction_tied_pairs,
            "truth_tie_tolerance": 0.0,
        },
        "component_losses": losses,
    }


def bootstrap_contrast(left, right, draws=9999, seed=20260811):
    components = sorted(set(left) & set(right))
    delta = np.asarray([left[key] - right[key] for key in components])
    rng = np.random.default_rng(seed)
    sampled = delta[rng.integers(0, len(delta), size=(draws, len(delta)))].mean(axis=1)
    return {"components": len(components), "mean_loss_reduction": float(delta.mean()),
            "one_sided_95_lcb": float(np.quantile(sampled, 0.05)),
            "pass": bool(np.quantile(sampled, 0.05) > 0)}


def permuted_y(y, rows, iteration):
    result = y.copy()
    groups = defaultdict(list)
    for index, row in enumerate(rows):
        groups[row["group"]].append(index)
    rng = np.random.default_rng(20260811 + iteration)
    for indices in groups.values():
        result[indices] = result[np.asarray(indices)[rng.permutation(len(indices))]]
    return result


def operator_prediction(operators, y):
    prediction = np.empty(len(y), dtype=np.float64)
    for train_index, test_index, operator in operators:
        prediction[test_index] = operator @ y[train_index]
    return prediction


def run():
    rows = build_rows()
    folds = component_folds(rows)
    y = np.asarray([row["y"] for row in rows])
    matrices = {name: np.stack([row[name] for row in rows])
                for name in ("nuisance", "esm", "tbasis", "tbasis_null")}
    results, operators = {}, {}
    for name, matrix in matrices.items():
        prediction, value_operators, alphas = cross_validate(matrix, y, rows, folds)
        results[name] = {**score(prediction, y, rows), "selected_alpha_by_fold": alphas}
        operators[name] = value_operators
    results["ligand_only_zero"] = score(np.zeros_like(y), y, rows)
    coupling = bootstrap_contrast(results["tbasis_null"]["component_losses"],
                                  results["tbasis"]["component_losses"], seed=20260811)
    additive = bootstrap_contrast(results["esm"]["component_losses"],
                                  results["tbasis"]["component_losses"], seed=20260812)
    zero = bootstrap_contrast(results["ligand_only_zero"]["component_losses"],
                              results["tbasis"]["component_losses"], seed=20260813)
    component_group_counts = defaultdict(set)
    for row in rows:
        component_group_counts[row["component"]].add(row["group"])
    component_group_counts = {key: len(value) for key, value in component_group_counts.items()}
    observed = coupling["mean_loss_reduction"]
    permutation_values = []
    for iteration in range(999):
        shuffled = permuted_y(y, rows, iteration)
        correct_prediction = operator_prediction(operators["tbasis"], shuffled)
        null_prediction = operator_prediction(operators["tbasis_null"], shuffled)
        indices = np.arange(len(rows))
        correct_loss = component_losses(correct_prediction, shuffled, rows, indices)
        null_loss = component_losses(null_prediction, shuffled, rows, indices)
        permutation_values.append(np.mean([null_loss[key] - correct_loss[key] for key in correct_loss]))
    extreme_count = int(sum(value >= observed for value in permutation_values))
    permutation_p = (1 + extreme_count) / 1000
    permutation_array = np.asarray(permutation_values)

    null_equal = np.asarray([
        np.allclose(row["tbasis"], row["tbasis_null"], rtol=0.0, atol=0.0) for row in rows
    ])
    group_indices = defaultdict(list)
    for index, row in enumerate(rows):
        group_indices[row["group"]].append(index)
    identity_groups = sum(bool(null_equal[indices].all()) for indices in group_indices.values())

    # Same closure and probe path, with a known decodable signal.
    planted_y = matrices["tbasis"][:, int(np.argmax(matrices["tbasis"].var(axis=0)))]
    planted_prediction, _, _ = cross_validate(matrices["tbasis"], planted_y, rows, folds)
    planted = score(planted_prediction, planted_y, rows)
    passed = (coupling["pass"] and additive["pass"] and permutation_p <= 0.05
              and planted["component_macro_r2_vs_zero"] > 0.5)
    result = {
        "schema": "MetaSieve.A1SelectivityProbe.v1",
        "TERMINAL_VERDICT": ("TBASIS_SELECTIVITY_SIGNAL_IDENTIFIED" if passed
                             else "TBASIS_SELECTIVITY_SIGNAL_NOT_IDENTIFIED"),
        "source_only": True, "main_v0_test_values_used": 0,
        "groups": len({row["group"] for row in rows}),
        "family_centered_rows": len(rows),
        "components": len({row["component"] for row in rows}),
        "fold_components": folds,
        "representations": results,
        "registered_gate": {
            "status": "FAIL_CLOSED_STATISTICAL_PROTOCOL_NOT_CONFIRMATORY",
            "a2_authorized": False,
            "reasons": [
                "single rewiring realization leaves substantial coupling fixed",
                "label permutation does not preserve repeated-family incidence",
                "permutation refits ridge coefficients at hyperparameters selected on observed labels",
            ],
        },
        "diagnostics": {
            "coupling_null_minus_tbasis": {
                **coupling,
                **contrast_sensitivity(results["tbasis_null"]["component_losses"],
                                       results["tbasis"]["component_losses"],
                                       component_group_counts),
            },
            "esm_additive_minus_tbasis": {
                **additive,
                **contrast_sensitivity(results["esm"]["component_losses"],
                                       results["tbasis"]["component_losses"],
                                       component_group_counts),
            },
            "zero_minus_tbasis": {
                **zero,
                **contrast_sensitivity(results["ligand_only_zero"]["component_losses"],
                                       results["tbasis"]["component_losses"],
                                       component_group_counts),
            },
            "single_rewiring_realization": {
                "fixed_rows": int(null_equal.sum()),
                "total_rows": len(rows),
                "fixed_row_fraction": float(null_equal.mean()),
                "identity_groups": identity_groups,
                "total_groups": len(group_indices),
            },
            "fixed_hyperparameter_groupwise_label_permutation": {
                "refits": 999,
                "seed_base": 20260811,
                "observed_statistic": observed,
                "alternative": "coupling_null_loss_minus_tbasis_loss > 0",
                "extreme_count": extreme_count,
                "uncalibrated_p": permutation_p,
                "minimum": float(permutation_array.min()),
                "q05": float(np.quantile(permutation_array, 0.05)),
                "median": float(np.median(permutation_array)),
                "q95": float(np.quantile(permutation_array, 0.95)),
                "maximum": float(permutation_array.max()),
                "confirmatory": False,
            },
        },
        "planted_positive_control": {
            key: value for key, value in planted.items() if key != "component_losses"
        },
        "a2_authorized": passed,
        "interpretation_boundary": "source component-held-out decodability; not downstream use or confirmation",
        "inputs": {"dependency_result_sha256": sha256(DEPENDENCY_RESULT),
                   "source_features_sha256": sha256(DEV / "source_features.npz"),
                   "contrast_groups_sha256": sha256(DEV / "source_contrast_groups.jsonl.gz")},
    }
    result["representations"] = {
        name: {key: value for key, value in detail.items() if key != "component_losses"}
        for name, detail in result["representations"].items()
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, allow_nan=False, indent=2, sort_keys=True) + "\n",
                   encoding="utf-8")
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, sort_keys=True))

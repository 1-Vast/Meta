"""Stage B frozen evaluation, with every Stage A correction applied.

Differences from `stageA_innerloop/evaluate.py`, each traceable to
`CORRECTION_AUDIT.md`:

1. matched-wrong support is anchored to the **pre-adaptation** support
   prediction, so every arm faces an identically corrupted adversary;
2. **every** control is reported for **every** arm, at that arm's own operating
   condition, plus the incremental (difference-of-differences) label dependence;
3. conditioning uses the corrected `alpha = 2*lr*(||h||^2 + 1)` and reports the
   exact per-query step effect rather than inferring it from support contraction;
4. shape is the query correction with its within-episode mean removed, not
   "whatever the weight produced";
5. wrong-protein is labelled a full-system perturbation;
6. the meta_test block is the audited seal record, verbatim.

`meta_val` is read exactly once here, after the checkpoints are frozen.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import QPSMPData                              # noqa: E402
from scripts.stageR0_retrieval_falsification import (                 # noqa: E402
    component_bootstrap, component_target_mean,
)
from scripts.stageR6_compare_arms import load_arm                     # noqa: E402
from scripts.train_qpsmp import (                                     # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    compact_episode, normalized_episode, training_label_scale,
)
from tools.research.stageA_innerloop.evaluate import (                # noqa: E402
    concordance, metrics_for, pearson, r_squared, spearman,
)
from tools.research.stageA_innerloop.train_meta import (              # noqa: E402
    align_atoms, encode_parts, episode_tensors,
)
from tools.research.stageB_complementary.arms import (                # noqa: E402
    ADAPTED_BY_MODE, InnerStepSizes, StageBAdaptation, predict,
)
from tools.research.stageB_complementary.residual import (            # noqa: E402
    centered_shape, conditioning_alpha,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
SUPPORT_SIZES = (0, 1, 2, 3, 5)
EVALUATION_SEED = 73101
QUERY_SIZE = 16
DRAWS = 2
BOOTSTRAP_DRAWS = 9999
BOOTSTRAP_SEED = 20260816
CLIFF_SIMILARITY = 0.6
CLIFF_GAP = 1.0
NOVELTY_HIGH = 0.4          # max Tanimoto to meta_train above this = "recall risk"

CONTROLS = ("no_adaptation", "permuted_support", "matched_wrong_support",
            "wrong_protein", "bias_only", "weight_only")
DECLARED_DIFFERENCES = {"mode", "inner_steps", "inner_lr", "learned_step",
                        "max_step", "adapted", "first_order"}


def verify_arms_are_matched(stage: Path, arms: tuple[str, ...],
                            allow_mixed_selection: bool = False) -> dict:
    """Refuse to score arms that differ outside the declared change.

    The checkpoint-selection rule is checked separately and by default must be
    identical across arms: a `meta_val`-selected arm and an internally selected
    one are not comparable, and mixing them silently would reintroduce exactly
    the leak this stage exists to remove. `allow_mixed_selection` is for the
    leakage diagnostic, which compares them on purpose.
    """
    configs, selections = {}, {}
    for arm in arms:
        path = stage / arm / "RESULT.json"
        if not path.exists():
            raise FileNotFoundError(f"missing arm result: {path}")
        payload = json.loads(path.read_text(encoding="utf-8"))
        configs[arm] = payload["config"]
        selections[arm] = payload["report"].get("selection_rule", "internal")
    reference = configs[arms[0]]
    for arm, config in configs.items():
        differing = {key for key in set(reference) | set(config)
                     if reference.get(key) != config.get(key)}
        unexpected = differing - DECLARED_DIFFERENCES
        if unexpected:
            raise ValueError(f"arms are not matched: {arm} differs in "
                             f"{sorted(unexpected)}")
    distinct = set(selections.values())
    if len(distinct) > 1 and not allow_mixed_selection:
        raise ValueError(
            f"arms use different checkpoint-selection rules {selections}; "
            f"pass --allow-mixed-selection only for the leakage diagnostic")
    for arm, config in configs.items():
        config["_selection_rule"] = selections[arm]
    return configs


def cliff_sign_accuracy(prediction, truth, similarity) -> float:
    """Orientation-invariant sign accuracy on high-similarity, large-gap pairs."""
    rows, cols = np.triu_indices(len(truth), 1)
    keep = (similarity[rows, cols] >= CLIFF_SIMILARITY) & \
           (np.abs(truth[rows] - truth[cols]) >= CLIFF_GAP)
    if not keep.any():
        return float("nan")
    delta = truth[rows][keep] - truth[cols][keep]
    gap = prediction[rows][keep] - prediction[cols][keep]
    return float((np.sign(gap) == np.sign(delta)).mean())


def tanimoto_matrix(fingerprint: np.ndarray) -> np.ndarray:
    intersection = fingerprint @ fingerprint.T
    total = fingerprint.sum(1)
    union = total[:, None] + total[None, :] - intersection
    return intersection / np.maximum(union, 1e-9)


def extended_metrics(prediction, truth, similarity) -> dict:
    values = metrics_for(prediction, truth)
    centered_p = prediction - prediction.mean()
    centered_t = truth - truth.mean()
    values["centered_mse_pk"] = float(((centered_p - centered_t) ** 2).mean())
    values["cliff_sign"] = cliff_sign_accuracy(prediction, truth, similarity)
    return values


def evaluate_arm(model, data, banks, label_scale, adaptation, steps, device,
                 novelty) -> list[dict]:
    rows = []
    scale, mean = label_scale.scale, label_scale.mean
    adapts_bias = "interaction_head.2.bias" in ADAPTED_BY_MODE[adaptation.mode]
    for k, bank in banks.items():
        for episode in bank:
            spec = episode.spec
            parts = align_atoms(episode_tensors(model, episode, device,
                                                torch.float32))
            truth = parts["query_y"].squeeze(0).cpu().numpy() * scale + mean
            fingerprint = parts["query_fingerprint"].squeeze(0).cpu().numpy()
            similarity = tanimoto_matrix(fingerprint)
            ligands = [data.cells[i]["ligand_id"] for i in spec.query]
            novelties = np.asarray([novelty(l) for l in ligands])

            task = encode_parts(model, parts)

            def record(condition, prediction, **extra):
                values = (prediction.squeeze(0).detach().cpu().numpy() * scale
                          + mean)
                rows.append({
                    "k": k, "component": spec.component, "target": spec.target,
                    "condition": condition,
                    **extended_metrics(values, truth, similarity),
                    "max_train_tanimoto": float(novelties.max()),
                    "mean_train_tanimoto": float(novelties.mean()),
                    **extra})

            with torch.no_grad():
                output = predict(model, parts, task, adaptation, steps)
                record("correct", output["prediction"],
                       meta_abs=float(output["meta"].abs().mean()),
                       meta_shape_abs=float(centered_shape(
                           output["meta"]).abs().mean()),
                       transport_abs=float(output["transport"].abs().mean()),
                       level_abs=float(output["level"].abs().mean()),
                       complementary_abs=float(
                           output["complementary"].abs().mean())
                       if output["complementary"].numel() else 0.0,
                       inner_loss_first=(output["inner_trace"][0]
                                         if output["inner_trace"] else None))
                record("zero_shot", output["zero_shot"])

                if k == 0:
                    continue

                # No adaptation: the arm's transport path alone.
                record("no_adaptation",
                       predict(model, parts, task, adaptation, steps,
                               disable_meta=True)["prediction"])
                # Meta term alone, transport suppressed.
                record("no_transport",
                       predict(model, parts, task, adaptation, steps,
                               disable_transport=True)["prediction"])

                support_y = parts["support_y"]
                # CORRECTION 1: anchor to the PRE-adaptation support prediction,
                # which is shared-initialization and therefore arm-independent.
                pre_support = output["support_zero"]
                raw_residual = (support_y - pre_support).detach()
                record("matched_wrong_support",
                       predict(model, parts, task, adaptation, steps,
                               support_y_override=support_y - 2.0 * raw_residual
                               )["prediction"])
                if support_y.shape[-1] > 1:
                    record("permuted_support",
                           predict(model, parts, task, adaptation, steps,
                                   support_y_override=support_y.roll(1, dims=-1)
                                   )["prediction"])

                # CORRECTION 4: level and shape are read off the correction
                # itself, not off which parameter produced it.
                correction = output["prediction"] - output["zero_shot"]
                record("level_only",
                       output["zero_shot"]
                       + correction.mean(-1, keepdim=True).expand_as(correction))
                record("shape_only",
                       output["zero_shot"] + centered_shape(correction))

                # The parameter-level ablation, reported separately because it
                # answers a different question than level/shape: which adapted
                # tensor carries the effect. `bias_only` is undefined for arms
                # that do not adapt the bias.
                if adaptation.mode != "T":
                    record("weight_only",
                           predict(model, parts, task, adaptation, steps,
                                   keep="weight")["prediction"])
                    if adapts_bias:
                        record("bias_only",
                               predict(model, parts, task, adaptation, steps,
                                       keep="bias")["prediction"])

                if adaptation.mode != "T":
                    alpha = conditioning_alpha(
                        task.support_hidden,
                        float(steps.weight_step()) if steps is not None
                        else adaptation.inner_lr,
                        adapt_bias=adapts_bias)
                    rows[-1]["alpha_mean"] = float(alpha.mean())
    return rows


def wrong_protein_rows(model, data, banks, label_scale, adaptation, steps,
                       device) -> list[dict]:
    """CORRECTION 5: a full-system perturbation, not an inner-loop test."""
    rows = []
    scale, mean = label_scale.scale, label_scale.mean
    for k, bank in banks.items():
        if k == 0:
            continue
        for episode in bank:
            spec = episode.spec
            parts = align_atoms(episode_tensors(model, episode, device,
                                                torch.float32))
            pooled, tokens, mask = data.protein_for_target(spec.donor_target)
            chemistry = data.protein_chemistry_for_target(spec.donor_target)
            wrong = dict(parts)
            wrong["protein_pooled"] = pooled.unsqueeze(0).to(device, torch.float32)
            wrong["protein_tokens"] = tokens.unsqueeze(0).to(device, torch.float32)
            wrong["protein_mask"] = mask.unsqueeze(0).to(device, torch.float32)
            wrong["protein_chemistry"] = chemistry.unsqueeze(0).to(device, torch.float32)
            truth = parts["query_y"].squeeze(0).cpu().numpy() * scale + mean
            fingerprint = parts["query_fingerprint"].squeeze(0).cpu().numpy()
            with torch.no_grad():
                task = encode_parts(model, wrong)
                prediction = predict(model, wrong, task, adaptation,
                                     steps)["prediction"]
            values = prediction.squeeze(0).cpu().numpy() * scale + mean
            rows.append({"k": k, "component": spec.component,
                         "target": spec.target, "condition": "wrong_protein",
                         **extended_metrics(values, truth,
                                            tanimoto_matrix(fingerprint)),
                         "max_train_tanimoto": 0.0, "mean_train_tanimoto": 0.0})
    return rows


def summarize(rows, condition, k, field, predicate=None) -> float:
    selected = [(r["component"], r["target"], r[field]) for r in rows
                if r["k"] == k and r["condition"] == condition
                and (predicate is None or predicate(r))]
    return component_target_mean(selected)


def paired(left_rows, right_rows, left_condition, right_condition, k, field,
           predicate=None) -> dict:
    index: dict[str, list] = {}
    for row in right_rows:
        if row["k"] == k and row["condition"] == right_condition:
            index.setdefault(row["target"], []).append(row)
    counters: dict[str, int] = {}
    pairs = []
    for row in left_rows:
        if row["k"] != k or row["condition"] != left_condition:
            continue
        if predicate is not None and not predicate(row):
            continue
        seat = counters.get(row["target"], 0)
        matches = index.get(row["target"], [])
        if seat >= len(matches):
            continue
        counters[row["target"]] = seat + 1
        left, right = row[field], matches[seat][field]
        if left is None or right is None or not (np.isfinite(left) and np.isfinite(right)):
            continue
        pairs.append((row["component"], row["target"], left - right))
    return component_bootstrap(pairs, BOOTSTRAP_DRAWS, BOOTSTRAP_SEED)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--arms", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--allow-mixed-selection", action="store_true",
                        help="leakage diagnostic only")
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    arguments = parser.parse_args()
    arms = tuple(arguments.arms)

    configs = verify_arms_are_matched(arguments.stage, arms,
                                      arguments.allow_mixed_selection)
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    label_scale = training_label_scale(data)

    from tools.research.a2_readiness_v2._donors import novelty_and_scaffold_strata
    scaffolds: dict[str, str] = {}
    novelty, _ = novelty_and_scaffold_strata(data, scaffolds)

    specs = data.fixed_nested_episode_banks(
        "meta_val", SUPPORT_SIZES, QUERY_SIZE, DRAWS, EVALUATION_SEED, None)
    banks = {k: tuple(compact_episode(normalized_episode(
        data.materialize(spec), label_scale)) for spec in items)
        for k, items in specs.items()}
    print(f"meta_val bank: {len(banks[0])} episodes per k, "
          f"{len({e.spec.target for e in banks[0]})} targets, "
          f"{len({e.spec.component for e in banks[0]})} components")

    all_rows: dict[str, list[dict]] = {}
    for arm in arms:
        checkpoint = arguments.stage / arm / "checkpoint.pt"
        payload = torch.load(checkpoint, map_location="cpu", weights_only=False)
        adaptation = StageBAdaptation.from_dict(payload["adaptation"])
        model, _, _ = load_arm(checkpoint, data, arguments.device)
        model.eval()
        steps = None
        if payload.get("inner_step_state") is not None:
            steps = InnerStepSizes(adaptation.inner_lr, adaptation.max_step)
            steps.load_state_dict(payload["inner_step_state"])
            steps = steps.to(arguments.device)
        rows = evaluate_arm(model, data, banks, label_scale, adaptation, steps,
                            arguments.device, novelty)
        rows.extend(wrong_protein_rows(model, data, banks, label_scale,
                                       adaptation, steps, arguments.device))
        all_rows[arm] = rows
        print(f"{arm}: {len(rows)} rows, mode={adaptation.mode}")
        del model
        if arguments.device.startswith("cuda"):
            torch.cuda.empty_cache()

    fields = ("mse_pk", "rmse_pk", "pearson", "spearman", "ci", "r2",
              "centered_mse_pk", "cliff_sign")
    low_recall = lambda row: row["max_train_tanimoto"] < NOVELTY_HIGH   # noqa: E731

    payload = {
        "schema": "MetaSieve.StageB.Complementary.Evaluation.v1",
        "date": "2026-08-17",
        "evidence_grade": "development, single seed, directional screening",
        "corrections_applied": [
            "matched-wrong anchored to pre-adaptation support",
            "all controls reported for all arms",
            "alpha includes the adapted bias",
            "shape measured as the centered query correction",
            "wrong-protein labelled a full-system perturbation",
            "audited meta_test seal string"],
        "checkpoint_selection": "meta_train internal-validation components only",
        "population": {"split": "meta_val", "episodes_per_k": len(banks[0]),
                       "targets": len({e.spec.target for e in banks[0]}),
                       "components": len({e.spec.component for e in banks[0]}),
                       "query_size": QUERY_SIZE, "draws": DRAWS,
                       "seed": EVALUATION_SEED},
        "arm_configs": configs,
        "arm_reports": {arm: json.loads(
            (arguments.stage / arm / "RESULT.json").read_text(encoding="utf-8")
        )["report"] for arm in arms},
        "arm_metrics": {}, "contrasts": {}, "counterfactuals": {},
        "novelty_strata": {}, "meta_test": data.seal_record(),
    }

    for arm in arms:
        payload["arm_metrics"][arm] = {
            str(k): {condition: {field: summarize(all_rows[arm], condition, k, field)
                                 for field in fields}
                     for condition in sorted({r["condition"] for r in all_rows[arm]
                                              if r["k"] == k})}
            for k in SUPPORT_SIZES}

    for left in arms:
        for right in arms:
            if left == right:
                continue
            payload["contrasts"][f"{left}_vs_{right}"] = {
                str(k): {field: paired(all_rows[left], all_rows[right],
                                       "correct", "correct", k, field)
                         for field in fields}
                for k in SUPPORT_SIZES}

    # CORRECTION 2: every control, for every arm.
    for arm in arms:
        payload["counterfactuals"][arm] = {
            str(k): {control: paired(all_rows[arm], all_rows[arm], control,
                                     "correct", k, "mse_pk")
                     for control in CONTROLS + ("no_transport", "level_only",
                                                "shape_only")}
            for k in (1, 2, 3, 5)}

    # Novelty stratum: does any benefit survive outside high-similarity recall?
    for left in arms:
        for right in arms:
            if left == right:
                continue
            payload["novelty_strata"][f"{left}_vs_{right}"] = {
                str(k): {"low_recall_mse": paired(
                    all_rows[left], all_rows[right], "correct", "correct", k,
                    "mse_pk", predicate=low_recall)}
                for k in (1, 2, 3, 5)}

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    with arguments.output.with_suffix(".rows.jsonl").open(
            "w", encoding="utf-8") as handle:
        for arm, rows in all_rows.items():
            for row in rows:
                handle.write(json.dumps({"arm": arm, **row}) + "\n")
    print(f"wrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

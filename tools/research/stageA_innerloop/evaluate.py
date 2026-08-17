"""Stage A frozen evaluation: the three arms on `meta_val`, read once.

Everything the preregistration promised to report, and nothing selected after
the fact. The arms are verified matched before any metric is computed: if the
recorded configs differ outside the declared fields the evaluator raises rather
than producing an interpretable-looking contrast.

`meta_test` is never constructed.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import torch
import torch.nn.functional as F

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
from tools.research.stageA_innerloop.inner_loop import AdaptationConfig  # noqa: E402
from tools.research.stageA_innerloop.train_meta import (              # noqa: E402
    align_atoms, encode_parts, episode_tensors, predict,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
ARMS = ("A0", "A1", "A2")
SUPPORT_SIZES = (0, 1, 2, 3, 5)
EVALUATION_SEED = 73101
QUERY_SIZE = 16
DRAWS = 2
BOOTSTRAP_DRAWS = 9999
BOOTSTRAP_SEED = 20260816
INNER_STEP_SWEEP = (0, 1, 2, 3)

# Fields allowed to differ between arms. Anything else is an unmatched control.
DECLARED_DIFFERENCES = {"arm", "inner_steps", "inner_lr", "task_selection"}


def verify_arms_are_matched(stage: Path) -> dict:
    """Refuse to score arms that differ in anything but the declared change."""
    configs = {}
    for arm in ARMS:
        path = stage / arm / "RESULT.json"
        if not path.exists():
            raise FileNotFoundError(f"missing arm result: {path}")
        configs[arm] = json.loads(path.read_text(encoding="utf-8"))["config"]
    reference = configs["A0"]
    for arm, config in configs.items():
        differing = {key for key in set(reference) | set(config)
                     if reference.get(key) != config.get(key)}
        unexpected = differing - DECLARED_DIFFERENCES
        if unexpected:
            raise ValueError(
                f"arms are not matched: {arm} differs from A0 in "
                f"{sorted(unexpected)}")
    if configs["A1"]["inner_steps"] < 1 or configs["A2"]["inner_steps"] < 1:
        raise ValueError("A1/A2 must actually adapt; inner_steps < 1")
    if configs["A0"]["inner_steps"] != 0:
        raise ValueError("A0 must be the no-inner-loop baseline")
    if not configs["A2"]["task_selection"]:
        raise ValueError("A2 must have task selection enabled")
    if configs["A1"]["task_selection"]:
        raise ValueError("A1 must use uniform sampling")
    return configs


# --- metrics ---------------------------------------------------------------

def centered(values: np.ndarray) -> np.ndarray:
    return values - values.mean()


def pearson(prediction: np.ndarray, truth: np.ndarray) -> float:
    p, t = centered(prediction), centered(truth)
    denominator = float(np.sqrt((p ** 2).sum()) * np.sqrt((t ** 2).sum()))
    return float((p * t).sum() / denominator) if denominator > 1e-12 else 0.0


def spearman(prediction: np.ndarray, truth: np.ndarray) -> float:
    def rank(values):
        order = values.argsort()
        ranks = np.empty(len(values), dtype=np.float64)
        ranks[order] = np.arange(len(values), dtype=np.float64)
        return ranks
    return pearson(rank(prediction), rank(truth))


def concordance(prediction: np.ndarray, truth: np.ndarray) -> float:
    rows, cols = np.triu_indices(len(truth), 1)
    delta = truth[rows] - truth[cols]
    comparable = delta != 0
    if not comparable.any():
        return float("nan")
    gap = prediction[rows] - prediction[cols]
    agree = np.sign(gap[comparable]) == np.sign(delta[comparable])
    tied = gap[comparable] == 0
    return float((agree.sum() + 0.5 * tied.sum()) / comparable.sum())


def r_squared(prediction: np.ndarray, truth: np.ndarray) -> float:
    residual = float(((truth - prediction) ** 2).sum())
    total = float(((truth - truth.mean()) ** 2).sum())
    return 1.0 - residual / total if total > 1e-12 else float("nan")


def metrics_for(prediction: np.ndarray, truth: np.ndarray) -> dict:
    error = prediction - truth
    return {"mse_pk": float((error ** 2).mean()),
            "rmse_pk": float(np.sqrt((error ** 2).mean())),
            "pearson": pearson(prediction, truth),
            "spearman": spearman(prediction, truth),
            "ci": concordance(prediction, truth),
            "r2": r_squared(prediction, truth)}


# --- the sweep -------------------------------------------------------------

def evaluate_arm(model, data, banks, label_scale, adaptation, device) -> list[dict]:
    """One row per (k, episode, condition)."""
    rows = []
    scale = label_scale.scale
    for k, bank in banks.items():
        for episode in bank:
            parts = align_atoms(episode_tensors(
                model, episode, device, torch.float32))
            truth = (parts["query_y"].squeeze(0).cpu().numpy() * scale
                     + label_scale.mean)
            spec = episode.spec

            def record(condition: str, prediction: torch.Tensor, **extra):
                values = prediction.squeeze(0).detach().cpu().numpy() * scale + label_scale.mean
                rows.append({"k": k, "component": spec.component,
                             "target": spec.target, "condition": condition,
                             **metrics_for(values, truth), **extra})

            # One encoder pass serves every condition that keeps this episode's
            # protein and ligands. Only the wrong-protein control changes the
            # encoding, and it is handled separately. Re-encoding per condition
            # cost 10x for identical tensors.
            task = encode_parts(model, parts)

            with torch.no_grad():
                # Inner-step sweep, correct labels.
                for steps in INNER_STEP_SWEEP:
                    config = AdaptationConfig(
                        inner_steps=steps, inner_lr=adaptation.inner_lr)
                    output = predict(model, parts, config, task=task)
                    record(f"steps{steps}", output["prediction"],
                           inner_loss_first=(output["inner_trace"][0]
                                             if output["inner_trace"] else None),
                           inner_loss_last=(output["inner_trace"][-1]
                                            if output["inner_trace"] else None))
                    if steps == adaptation.inner_steps:
                        record("pre_adaptation", output["pre_adaptation_query"])
                        base_output = output

                if k == 0:
                    continue

                config = adaptation
                # Level-vs-shape decomposition of the adaptation.
                for keep in ("bias", "weight"):
                    record(f"keep_{keep}",
                           predict(model, parts, config, keep=keep,
                                   task=task)["prediction"])

                # Counterfactual supports, defined exactly as the accepted
                # recipe defines them. Each differs from the real arm in the
                # labels alone; everything else is the same episode.
                support_y = parts["support_y"]
                if support_y.shape[-1] > 1:
                    # A non-identity cyclic binding. A one-element permutation
                    # is the identity, so at k=1 the permuted control does not
                    # exist and only the magnitude-matched flip is meaningful.
                    record("permuted_support",
                           predict(model, parts, config, task=task,
                                   support_y_override=support_y.roll(1, dims=-1)
                                   )["prediction"])
                # Equal-magnitude residual flip: same label scale, wrong sign
                # of evidence. Defined at every k.
                raw_residual = (
                    support_y - base_output["post_adaptation_support"]).detach()
                record("matched_wrong_support",
                       predict(model, parts, config, task=task,
                               support_y_override=support_y - 2.0 * raw_residual
                               )["prediction"])
                record("no_adaptation",
                       predict(model, parts, AdaptationConfig(inner_steps=0),
                               task=task)["prediction"])
    return rows


def wrong_protein_rows(model, data, banks, label_scale, adaptation,
                       device) -> list[dict]:
    """Substitute the donor protein; keep the ligands and labels."""
    rows = []
    scale = label_scale.scale
    for k, bank in banks.items():
        for episode in bank:
            spec = episode.spec
            parts = align_atoms(episode_tensors(
                model, episode, device, torch.float32))
            pooled, tokens, mask = data.protein_for_target(spec.donor_target)
            chemistry = data.protein_chemistry_for_target(spec.donor_target)
            wrong = dict(parts)
            wrong["protein_pooled"] = pooled.unsqueeze(0).to(device, torch.float32)
            wrong["protein_tokens"] = tokens.unsqueeze(0).to(device, torch.float32)
            wrong["protein_mask"] = mask.unsqueeze(0).to(device, torch.float32)
            wrong["protein_chemistry"] = chemistry.unsqueeze(0).to(device, torch.float32)
            truth = (parts["query_y"].squeeze(0).cpu().numpy() * scale
                     + label_scale.mean)
            with torch.no_grad():
                prediction = predict(model, wrong, adaptation)["prediction"]
            values = prediction.squeeze(0).cpu().numpy() * scale + label_scale.mean
            rows.append({"k": k, "component": spec.component,
                         "target": spec.target, "condition": "wrong_protein",
                         **metrics_for(values, truth)})
    return rows


def summarize(rows: list[dict], condition: str, k: int, field: str) -> float:
    selected = [(r["component"], r["target"], r[field]) for r in rows
                if r["k"] == k and r["condition"] == condition]
    return component_target_mean(selected)


def paired(left_rows, right_rows, left_condition: str, right_condition: str,
           k: int, field: str) -> dict:
    """Component-paired `left - right` on matching (target, k) cells.

    The two conditions are named separately so the same helper serves both
    between-arm contrasts (same condition, different arm) and within-arm
    counterfactuals (same arm, different condition). Collapsing them into one
    argument silently made every counterfactual an exact zero.
    """
    index: dict[str, list[dict]] = {}
    for row in right_rows:
        if row["k"] == k and row["condition"] == right_condition:
            index.setdefault(row["target"], []).append(row)
    counters: dict[str, int] = {}
    pairs = []
    for row in left_rows:
        if row["k"] != k or row["condition"] != left_condition:
            continue
        seat = counters.get(row["target"], 0)
        matches = index.get(row["target"], [])
        if seat >= len(matches):
            continue
        counters[row["target"]] = seat + 1
        pairs.append((row["component"], row["target"],
                      row[field] - matches[seat][field]))
    return component_bootstrap(pairs, BOOTSTRAP_DRAWS, BOOTSTRAP_SEED)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    arguments = parser.parse_args()

    configs = verify_arms_are_matched(arguments.stage)
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    label_scale = training_label_scale(data)
    specs = data.fixed_nested_episode_banks(
        "meta_val", SUPPORT_SIZES, QUERY_SIZE, DRAWS, EVALUATION_SEED, None)
    banks = {k: tuple(compact_episode(normalized_episode(
        data.materialize(spec), label_scale)) for spec in items)
        for k, items in specs.items()}
    print(f"meta_val bank: {len(banks[0])} episodes per k, "
          f"{len({e.spec.target for e in banks[0]})} targets, "
          f"{len({e.spec.component for e in banks[0]})} components")

    all_rows: dict[str, list[dict]] = {}
    for arm in ARMS:
        checkpoint = arguments.stage / arm / "checkpoint.pt"
        payload = torch.load(checkpoint, map_location="cpu", weights_only=False)
        adaptation = AdaptationConfig.from_dict(payload["adaptation"])
        model, _, _ = load_arm(checkpoint, data, arguments.device)
        model.eval()
        rows = evaluate_arm(model, data, banks, label_scale, adaptation,
                            arguments.device)
        rows.extend(wrong_protein_rows(model, data, banks, label_scale,
                                       adaptation, arguments.device))
        all_rows[arm] = rows
        print(f"{arm}: {len(rows)} rows, adaptation={adaptation.to_dict()}")
        del model
        if arguments.device.startswith("cuda"):
            torch.cuda.empty_cache()

    fields = ("mse_pk", "rmse_pk", "pearson", "spearman", "ci", "r2")
    payload = {
        "schema": "MetaSieve.StageA.InnerLoop.Evaluation.v1",
        "date": "2026-08-17",
        "evidence_grade": "development, single seed, directional screening",
        "not_a_claim": ("not a performance claim, not a reproduction of "
                        "AdaMBind, not a confirmation; meta_test unread"),
        "population": {"split": "meta_val", "episodes_per_k": len(banks[0]),
                       "targets": len({e.spec.target for e in banks[0]}),
                       "components": len({e.spec.component for e in banks[0]}),
                       "query_size": QUERY_SIZE, "draws": DRAWS,
                       "seed": EVALUATION_SEED},
        "arm_configs": configs,
        "arm_reports": {arm: json.loads(
            (arguments.stage / arm / "RESULT.json").read_text(encoding="utf-8")
        )["report"] for arm in ARMS},
        "arm_metrics": {}, "contrasts": {}, "counterfactuals": {},
        "meta_test": data.seal_record(),
    }
    for arm in ARMS:
        payload["arm_metrics"][arm] = {
            str(k): {condition: {field: summarize(all_rows[arm], condition, k, field)
                                 for field in fields}
                     for condition in sorted({r["condition"] for r in all_rows[arm]
                                              if r["k"] == k})}
            for k in SUPPORT_SIZES}

    # Between-arm contrasts are scored on each arm's own operating condition:
    # A0 has no inner loop, so its k>0 row is `steps0`, while A1/A2 adapt.
    # Scoring A0 under `steps1` would credit it with a mechanism it was never
    # trained to use.
    def condition_for(arm: str, k: int) -> str:
        if k == 0:
            return "steps0"
        return "steps0" if configs[arm]["inner_steps"] == 0 else "steps1"

    for left, right in (("A1", "A0"), ("A2", "A1"), ("A2", "A0")):
        payload["contrasts"][f"{left}_vs_{right}"] = {
            str(k): {field: paired(all_rows[left], all_rows[right],
                                   condition_for(left, k),
                                   condition_for(right, k), k, field)
                     for field in fields}
            for k in SUPPORT_SIZES}

    # Counterfactuals: control minus correct, within one arm. Positive means
    # the control is worse, which is the direction a real mechanism predicts.
    for arm in ("A1", "A2"):
        payload["counterfactuals"][arm] = {
            str(k): {control: paired(all_rows[arm], all_rows[arm], control,
                                     "steps1", k, "mse_pk")
                     for control in ("permuted_support",
                                     "matched_wrong_support",
                                     "no_adaptation", "wrong_protein",
                                     "keep_bias", "keep_weight")}
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

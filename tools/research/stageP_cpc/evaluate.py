"""Stage P frozen evaluation: does centered supervision buy protein-conditioned
within-target ordering?

Reports the two protein conditions **separately** and decomposes the gap,
because a contrastive objective can widen correct-minus-wrong by damaging the
donor arm, which is an aversion rather than specificity. Gate P1 is about
`r_correct` alone; the gap is diagnostic only.

Everything here is frozen evaluation on `meta_val`, read once. Training donors
came from `meta_train` only (verified by construction in
`QPSMPData.draw_episode`); evaluation donors are the frozen `meta_val`
stratified rule whitened on `meta_train`.
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
from scripts.train_level_shape import normalized                      # noqa: E402
from scripts.train_qpsmp import (                                     # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    centered_protein_contrast, compact_episode, training_label_scale,
)
from tools.research.a2_readiness_v2 import _frozen                    # noqa: E402
from tools.research.a2_readiness_v2._donors import stratified_donors  # noqa: E402

SEEDS = (20260815, 20260816, 20260817)
ARMS = ("A0repro", "CPCoverdrive")
SUPPORT_SIZES = (0, 1, 2, 3, 5)
BRANCHES = ("protein_encoder", "grammar", "ligand_encoder", "embed",
            "interaction_head", "protein_head", "contact_weight", "transport")


def centered(values: np.ndarray) -> np.ndarray:
    return values - values.mean()


def correlation(prediction: np.ndarray, truth: np.ndarray) -> float:
    p, t = centered(prediction), centered(truth)
    denominator = float(np.sqrt((p ** 2).mean()) * np.sqrt((t ** 2).mean()))
    return float((p * t).mean() / denominator) if denominator > 1e-12 else 0.0


def concordance(prediction: np.ndarray, truth: np.ndarray) -> float:
    rows, cols = np.triu_indices(len(truth), 1)
    delta = truth[rows] - truth[cols]
    keep = delta != 0
    if not keep.any():
        return float("nan")
    signed = np.sign(delta[keep]) * (prediction[rows] - prediction[cols])[keep]
    return float((signed > 0).mean() + 0.5 * (signed == 0).mean())


def spearman(prediction: np.ndarray, truth: np.ndarray) -> float:
    def rank(values):
        order = values.argsort()
        ranks = np.empty(len(values), dtype=float)
        ranks[order] = np.arange(len(values), dtype=float)
        return ranks
    return correlation(rank(prediction), rank(truth))


def protein_inputs(data, target: str, device: str, dtype):
    pooled, tokens, mask = data.protein_for_target(target)
    chemistry = data.protein_chemistry_for_target(target)
    return [pooled.to(device, dtype).unsqueeze(0),
            tokens.to(device, dtype).unsqueeze(0),
            mask.to(device, dtype).unsqueeze(0),
            chemistry.to(device, dtype).unsqueeze(0)]


def predict(model, data, episode, target: str, support_size: int,
            device: str, dtype):
    """Prediction under an arbitrary protein, at the given support size."""
    pooled, tokens, mask, chemistry = protein_inputs(data, target, device, dtype)
    support_atoms = episode.support_atoms[:support_size].to(device, dtype).unsqueeze(0)
    support_bonds = episode.support_bonds[:support_size].to(device, dtype).unsqueeze(0)
    support_mask = episode.support_mask[:support_size].to(device, dtype).unsqueeze(0)
    support_y = episode.support_y[:support_size].to(device, dtype).unsqueeze(0)
    support_fp = episode.support_fingerprint[:support_size].to(device, dtype).unsqueeze(0)
    output = model(
        pooled, tokens, mask, support_atoms, support_bonds, support_mask,
        support_y,
        episode.query_atoms.to(device, dtype).unsqueeze(0),
        episode.query_bonds.to(device, dtype).unsqueeze(0),
        episode.query_mask.to(device, dtype).unsqueeze(0),
        adapt=support_size > 0, protein_chemistry=chemistry,
        support_fingerprint=support_fp,
        query_fingerprint=episode.query_fingerprint.to(device, dtype).unsqueeze(0))
    return (output.prediction.squeeze(0).detach().float().cpu().numpy(),
            output.zero_shot.squeeze(0).detach().float().cpu().numpy())


def branch_gradient_norms(model, data, episode, donor: str, scale,
                          device: str, dtype) -> dict:
    """Gradient of the centered protein contrast, by branch.

    `protein_head` must be identically zero: the whole point of centering is
    that the additive level branch cannot satisfy the term.
    """
    model.zero_grad(set_to_none=True)
    pooled, tokens, mask, chemistry = protein_inputs(
        data, episode.spec.target, device, dtype)
    empty = episode.query_atoms[:0].to(device, dtype).unsqueeze(0)

    def zero_shot(parts):
        p, t, m, c = parts
        return model(
            p, t, m, empty, episode.query_bonds[:0].to(device, dtype).unsqueeze(0),
            episode.query_mask[:0].to(device, dtype).unsqueeze(0),
            torch.zeros(1, 0, device=device, dtype=dtype),
            episode.query_atoms.to(device, dtype).unsqueeze(0),
            episode.query_bonds.to(device, dtype).unsqueeze(0),
            episode.query_mask.to(device, dtype).unsqueeze(0),
            adapt=False, protein_chemistry=c,
            support_fingerprint=episode.query_fingerprint[:0].to(device, dtype).unsqueeze(0),
            query_fingerprint=episode.query_fingerprint.to(device, dtype).unsqueeze(0)
        ).zero_shot

    truth = episode.query_y.to(device, dtype).unsqueeze(0)
    loss = centered_protein_contrast(
        zero_shot([pooled, tokens, mask, chemistry]),
        zero_shot(protein_inputs(data, donor, device, dtype)), truth, 0.1)
    loss.backward()
    out = {}
    for branch in BRANCHES:
        module = getattr(model, branch, None)
        if module is None:
            continue
        total = sum(float(p.grad.norm()) ** 2 for p in module.parameters()
                    if p.grad is not None)
        out[branch] = float(total ** 0.5)
    model.zero_grad(set_to_none=True)
    return out


EXPECTED_DIFFERENCES = {"protein_contrast_form", "protein_contrast_loss_weight"}


def verify_arms_are_matched(stage: Path) -> dict:
    """Refuse to score arms that differ in anything but the intended change.

    A "matched control" that quietly differs in seed, budget, architecture or
    split is not a control, and the resulting contrast would be uninterpretable.
    This runs before any metric is computed and raises rather than warning.
    """
    report = {}
    for seed in SEEDS:
        configs = {}
        for arm in ARMS:
            path = stage / f"{arm}_seed{seed}" / "RESULT.json"
            if not path.is_file():
                break
            configs[arm] = json.loads(path.read_text(encoding="utf-8"))["config"]
        if len(configs) != len(ARMS):
            continue
        left, right = configs["A0repro"], configs["CPCoverdrive"]
        differ = {key for key in left if left[key] != right.get(key)} - {"output"}
        if differ != EXPECTED_DIFFERENCES:
            raise ValueError(
                f"seed {seed}: arms differ in {sorted(differ)}, expected exactly "
                f"{sorted(EXPECTED_DIFFERENCES)}. The control is not matched.")
        if left["protein_contrast_form"] != "uncentered" or \
                right["protein_contrast_form"] != "centered":
            raise ValueError(f"seed {seed}: the contrast form flag did not take "
                             "effect in the recorded configs")
        report[str(seed)] = {
            "differing_fields": sorted(differ),
            "A0repro": {"form": left["protein_contrast_form"],
                        "weight": left["protein_contrast_loss_weight"]},
            "CPCoverdrive": {"form": right["protein_contrast_form"],
                             "weight": right["protein_contrast_loss_weight"]},
        }
    if not report:
        raise ValueError("no complete seed pair found under --stage")
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    matched = verify_arms_are_matched(arguments.stage)
    print(f"arm matching verified for {len(matched)} seed pair(s)")

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=_frozen.SPLIT_DIRECTORY)
    scale = training_label_scale(data)
    donors = stratified_donors(data, "meta_val", _frozen.DONOR_POOL,
                               _frozen.WHITENING_POOL)
    banks = data.fixed_nested_episode_banks(
        "meta_val", SUPPORT_SIZES, _frozen.QUERY_SIZE, 1,
        _frozen.EVALUATION_SEED, None)
    episodes = {spec.target: compact_episode(normalized(data.materialize(spec), scale))
                for spec in banks[max(SUPPORT_SIZES)]}

    rows: list[dict] = []
    shifts: dict[tuple[str, int], dict] = {}
    gradients: dict[tuple[str, int], dict] = {}

    for arm in ARMS:
        for seed in SEEDS:
            checkpoint = arguments.stage / f"{arm}_seed{seed}" / "checkpoint.pt"
            if not checkpoint.exists():
                print(f"  MISSING {checkpoint}")
                continue
            model, _, _ = load_arm(checkpoint, data, arguments.device)
            model.eval()
            dtype = next(model.parameters()).dtype
            per_target_shift = {}
            with torch.no_grad():
                for spec in banks[max(SUPPORT_SIZES)]:
                    episode = episodes[spec.target]
                    donor = donors[spec.target]["nearest"][0]
                    truth = (episode.query_y.numpy() * scale.scale + scale.mean)
                    for k in SUPPORT_SIZES:
                        correct, correct_zero = predict(
                            model, data, episode, spec.target, k,
                            arguments.device, dtype)
                        wrong, wrong_zero = predict(
                            model, data, episode, donor, k,
                            arguments.device, dtype)
                        correct = correct * scale.scale + scale.mean
                        wrong = wrong * scale.scale + scale.mean
                        row = {"arm": arm, "seed": seed, "k": k,
                               "target": spec.target, "component": spec.component}
                        for tag, values in (("correct", correct), ("wrong", wrong)):
                            row[f"mse_{tag}"] = float(((values - truth) ** 2).mean())
                            row[f"centered_mse_{tag}"] = float(
                                ((centered(values) - centered(truth)) ** 2).mean())
                            row[f"r_{tag}"] = correlation(values, truth)
                            row[f"ci_{tag}"] = concordance(values, truth)
                            row[f"spearman_{tag}"] = spearman(values, truth)
                        row["calibration"] = float(
                            (correct.mean() - truth.mean()) ** 2)
                        if k == 0:
                            shift = centered(correct_zero) - centered(wrong_zero)
                            row["shift_alignment"] = correlation(shift, centered(truth))
                            per_target_shift[spec.target] = (
                                shift * scale.scale).astype(float)
                        rows.append(row)
            shifts[(arm, seed)] = per_target_shift
            sample = banks[max(SUPPORT_SIZES)][:12]
            gradients[(arm, seed)] = {
                branch: float(np.mean([g[branch] for g in (
                    branch_gradient_norms(model, data, episodes[s.target],
                                          donors[s.target]["nearest"][0],
                                          scale, arguments.device, dtype)
                    for s in sample)]))
                for branch in BRANCHES if hasattr(model, branch)}
            del model
            if arguments.device.startswith("cuda"):
                torch.cuda.empty_cache()
            print(f"  scored {arm} seed {seed}")

    def weighted(subset, field):
        return component_target_mean(
            (r["component"], r["target"], r.get(field)) for r in subset)

    payload: dict = {
        "schema": "MetaSieve.StagePCPC.v1",
        "split": "meta_val", "seeds": list(SEEDS), "arms": list(ARMS),
        "frozen_design": _frozen.frozen_manifest(),
        "meta_test": data.seal_record(),
        "arm_matching": matched,
        "arm_metrics": {}, "gates": {}, "gradients": {},
    }
    for arm in ARMS:
        block = {}
        for k in SUPPORT_SIZES:
            subset = [r for r in rows if r["arm"] == arm and r["k"] == k]
            if not subset:
                continue
            block[str(k)] = {
                field: weighted(subset, field) for field in (
                    "mse_correct", "centered_mse_correct", "r_correct",
                    "ci_correct", "spearman_correct",
                    "mse_wrong", "centered_mse_wrong", "r_wrong",
                    "ci_wrong", "spearman_wrong", "calibration")}
            if k == 0:
                block[str(k)]["shift_alignment"] = weighted(
                    subset, "shift_alignment")
        payload["arm_metrics"][arm] = block
        payload["gradients"][arm] = {
            branch: float(np.mean([gradients[(arm, s)][branch] for s in SEEDS
                                   if (arm, s) in gradients]))
            for branch in BRANCHES
            if any((arm, s) in gradients and branch in gradients[(arm, s)]
                   for s in SEEDS)}

    # --- gates ------------------------------------------------------------
    def paired(field: str, k: int) -> dict:
        left = [r for r in rows if r["arm"] == "CPCoverdrive" and r["k"] == k]
        right = [r for r in rows if r["arm"] == "A0repro" and r["k"] == k]
        index = {(r["seed"], r["target"]): r for r in right}
        pairs = [(a["component"], a["target"],
                  a[field] - index[(a["seed"], a["target"])][field])
                 for a in left if (a["seed"], a["target"]) in index]
        return component_bootstrap(pairs, _frozen.BOOTSTRAP_DRAWS,
                                   _frozen.BOOTSTRAP_SEED)

    improvement = paired("r_correct", 0)
    degradation = paired("r_wrong", 0)
    payload["gates"] = {
        "P1_correct_ordering_improvement": dict(
            improvement,
            minimum_effect=_frozen.SMALLEST_EFFECT_OF_INTEREST_R,
            passes=bool(improvement["lo"] > 0
                        and improvement["mean"] >= _frozen.SMALLEST_EFFECT_OF_INTEREST_R)),
        "donor_change_r_wrong": degradation,
        "gap_decomposition": {
            "improvement_term": improvement["mean"],
            "degradation_term": -degradation["mean"],
            "degradation_share": (
                float(-degradation["mean"]
                      / (abs(improvement["mean"]) + abs(degradation["mean"])))
                if (abs(improvement["mean"]) + abs(degradation["mean"])) > 1e-12
                else float("nan")),
            "note": ("a wider correct-minus-wrong gap produced mainly by the "
                     "degradation term is donor damage, not specificity"),
        },
        "k0_mse_change": paired("mse_correct", 0),
        "k0_ci_change": paired("ci_correct", 0),
        "k0_spearman_change": paired("spearman_correct", 0),
        "k0_calibration_change": paired("calibration", 0),
    }

    # Seed-to-seed reproducibility of the protein-induced shift.
    for arm in ARMS:
        cosines = []
        available = [s for s in SEEDS if (arm, s) in shifts]
        for target in {t for s in available for t in shifts[(arm, s)]}:
            vectors = [shifts[(arm, s)][target] for s in available
                       if target in shifts[(arm, s)]]
            for i in range(len(vectors)):
                for j in range(i + 1, len(vectors)):
                    a, b = vectors[i], vectors[j]
                    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
                    if denominator > 1e-12:
                        cosines.append(float(a @ b / denominator))
        payload["gates"].setdefault("seed_cosine", {})[arm] = {
            "mean": float(np.mean(cosines)) if cosines else float("nan"),
            "sd": float(np.std(cosines)) if cosines else float("nan"),
            "pairs": len(cosines)}

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    with arguments.output.with_suffix(".rows.jsonl").open(
            "w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row) + "\n")
    report(payload)
    print(f"\nwrote {arguments.output}")
    return 0


def report(payload: dict) -> None:
    print(f"\n{'arm':<14}{'k':>3}{'MSE(c)':>9}{'r(c)':>8}{'CI(c)':>7}"
          f"{'r(w)':>8}{'CI(w)':>7}{'align':>8}")
    for arm, block in payload["arm_metrics"].items():
        for k, cell in block.items():
            # `shift_alignment` is a k=0 quantity; blank rather than NaN
            # elsewhere, so the table does not look broken.
            alignment = (f"{cell['shift_alignment']:>+8.3f}"
                         if "shift_alignment" in cell else f"{'':>8}")
            print(f"{arm if k == '0' else '':<14}{k:>3}"
                  f"{cell['mse_correct']:>9.4f}{cell['r_correct']:>+8.3f}"
                  f"{cell['ci_correct']:>7.3f}{cell['r_wrong']:>+8.3f}"
                  f"{cell['ci_wrong']:>7.3f}{alignment}")
    gates = payload["gates"]
    p1 = gates["P1_correct_ordering_improvement"]
    print(f"\nP1 correct-protein k=0 ordering: {p1['mean']:+.4f} "
          f"[{p1['lo']:+.4f},{p1['hi']:+.4f}]  "
          f"{'PASS' if p1['passes'] else 'FAIL'} "
          f"(needs lo>0 and mean>={p1['minimum_effect']})")
    decomposition = gates["gap_decomposition"]
    print(f"   improvement {decomposition['improvement_term']:+.4f} | "
          f"donor degradation {decomposition['degradation_term']:+.4f} | "
          f"degradation share {decomposition['degradation_share']:.2f}")
    for name in ("k0_mse_change", "k0_ci_change", "k0_spearman_change",
                 "k0_calibration_change"):
        interval = gates[name]
        print(f"   {name:<24}{interval['mean']:+.4f} "
              f"[{interval['lo']:+.4f},{interval['hi']:+.4f}]")
    for arm, cell in gates["seed_cosine"].items():
        print(f"   seed cosine {arm:<14}{cell['mean']:+.4f} "
              f"(sd {cell['sd']:.3f}, {cell['pairs']} pairs)")
    print("\ngradient norms of the centered contrast, by branch:")
    for arm, block in payload["gradients"].items():
        print(f"  {arm}")
        for branch, value in block.items():
            print(f"    {branch:<20}{value:.4e}")


if __name__ == "__main__":
    raise SystemExit(main())

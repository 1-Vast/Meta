"""Phase 1: is the protein-inert ordering real, or an artifact of the v1 design?

No training, no gradient, forward passes only. `meta_test` is unreachable (the
repaired fail-closed default).

The v1 probe (`tools/research/a2_readiness/branch_ordering_probe.py`) reported
that A0's interaction branch moves 0.2150 pK in level and 0.0007 pK in centered
ordering under a wrong protein, and that one randomly initialised model moved
110x more. Four things could have produced that pattern without the causal
story being true, and this probe measures each of them:

1. **donor choice.** v1 used only the nearest cross-component donor. Phase 1
   draws five donors per target at fixed quantiles of whitened protein
   similarity, so perturbation magnitude is a variable, not a confound.
2. **measurement floor.** v1 never measured how much the number moves when
   *nothing* changes. Phase 1 adds an identical-protein substitution (must be
   exactly zero) and a repeated forward pass (the numerical floor).
3. **metric construction.** Phase 1 adds shuffled-label, shuffled-protein and
   random-ligand-panel arms, so a collapse that the metric produces on its own
   is visible.
4. **random-init interpretation.** A nonzero shift at initialisation is not
   evidence of useful capacity. Phase 1 builds ten independent inits and
   separates four different things a shift could be:

   * *arbitrary* — magnitude of the centered shift;
   * *aligned* — correlation of the shift with the true within-target label
     differences;
   * *reproducible* — correlation between two inits' shift vectors on the same
     target;
   * *useful* — the branch's own within-target correlation with truth.

   Only the last two together would support "the architecture expresses
   protein-conditioned ordering and training removes it".

Seeds are aggregated *within target* before the component bootstrap, because
three checkpoints scored on the same 41 targets and the same query panel are
not three biological samples.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import QPSMPData                            # noqa: E402
from scripts.stageR0_retrieval_falsification import (               # noqa: E402
    component_bootstrap, component_target_mean, murcko_scaffolds,
)
from scripts.train_level_shape import normalized                    # noqa: E402
from scripts.train_qpsmp import (                                   # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    compact_episode, training_label_scale,
)
from tools.research.a2_readiness_v2 import _frozen                  # noqa: E402
from tools.research.a2_readiness_v2._arms import random_arm, trained_arm  # noqa: E402
from tools.research.a2_readiness_v2._donors import (                # noqa: E402
    novelty_and_scaffold_strata, stratified_donors,
)


# --------------------------------------------------------------------------
# metrics
# --------------------------------------------------------------------------

def centered(values: np.ndarray) -> np.ndarray:
    return values - values.mean()


def correlation(prediction: np.ndarray, truth: np.ndarray) -> float:
    p, y = centered(prediction), centered(truth)
    denominator = float(np.sqrt((p ** 2).mean()) * np.sqrt((y ** 2).mean()))
    return float((p * y).mean() / denominator) if denominator > 1e-12 else 0.0


def concordance(prediction: np.ndarray, truth: np.ndarray) -> float:
    rows, cols = np.triu_indices(len(truth), 1)
    delta = truth[rows] - truth[cols]
    comparable = delta != 0
    if not comparable.any():
        return float("nan")
    signed = np.sign(delta[comparable]) * (
        prediction[rows] - prediction[cols])[comparable]
    return float((signed > 0).mean() + 0.5 * (signed == 0).mean())


# --------------------------------------------------------------------------
# forward passes
# --------------------------------------------------------------------------

def protein_inputs(data, target: str, device: str, dtype):
    pooled, tokens, mask = data.protein_for_target(target)
    chemistry = data.protein_chemistry_for_target(target)
    return [pooled.to(device, dtype).unsqueeze(0),
            tokens.to(device, dtype).unsqueeze(0),
            mask.to(device, dtype).unsqueeze(0),
            chemistry.to(device, dtype).unsqueeze(0)]


def scramble_slots(parts: list, seed: int) -> list:
    """Permute the residue-slot axis: same composition, destroyed arrangement."""
    pooled, tokens, mask, chemistry = parts
    order = np.random.default_rng(seed).permutation(tokens.shape[1])
    index = torch.as_tensor(order, device=tokens.device, dtype=torch.long)
    return [pooled, tokens[:, index], mask[:, index], chemistry[:, index]]


def branches(model, parts: list, query_atoms, query_bonds, query_mask,
             query_fingerprint, device, dtype) -> dict:
    """Zero-shot branch decomposition. `adapt=False` reads no support label."""
    pooled, tokens, mask, chemistry = parts
    # An empty support keeps the batch axis: [1, 0, ...], not [0, Q, ...].
    output = model(
        pooled, tokens, mask,
        query_atoms[:, :0], query_bonds[:, :0], query_mask[:, :0],
        torch.zeros(1, 0, device=device, dtype=dtype),
        query_atoms, query_bonds, query_mask,
        adapt=False, protein_chemistry=chemistry,
        support_fingerprint=query_fingerprint[:, :0],
        query_fingerprint=query_fingerprint)
    # `.detach()` so the helper is usable outside `torch.no_grad()` too — the
    # structural tests call it directly.
    endpoint = output.zero_shot.squeeze(0).detach().float().cpu().numpy()
    ligand = output.ligand_only.squeeze(0).detach().float().cpu().numpy()
    additive = output.additive.squeeze(0).detach().float().cpu().numpy()
    return {"full": endpoint, "ligand_only": ligand,
            "protein_only": additive - ligand,
            "interaction": endpoint - additive}


def query_tensors(episode, device, dtype, indices=None):
    atoms = episode.query_atoms
    bonds = episode.query_bonds
    mask = episode.query_mask
    fingerprint = episode.query_fingerprint
    if indices is not None:
        atoms, bonds = atoms[indices], bonds[indices]
        mask, fingerprint = mask[indices], fingerprint[indices]
    return (atoms.to(device, dtype).unsqueeze(0),
            bonds.to(device, dtype).unsqueeze(0),
            mask.to(device, dtype).unsqueeze(0),
            fingerprint.to(device, dtype).unsqueeze(0))


# --------------------------------------------------------------------------
# per-episode measurement
# --------------------------------------------------------------------------

def measure(model, data, spec, episode, donors, scale, device, dtype,
            control_seed: int, foreign_episode) -> dict:
    """Every quantity Phase 1 needs from one (arm, seed, target)."""
    query = query_tensors(episode, device, dtype)
    truth = episode.query_y.numpy() * scale.scale + scale.mean
    correct_parts = protein_inputs(data, spec.target, device, dtype)
    correct = branches(model, correct_parts, *query, device, dtype)

    row: dict = {"target": spec.target, "component": spec.component,
                 "queries": int(len(truth)), "sd_truth": float(truth.std())}
    for branch in ("full", "ligand_only", "interaction"):
        row[f"r_{branch}"] = correlation(correct[branch], truth)
        row[f"ci_{branch}"] = concordance(correct[branch], truth)
        row[f"sd_{branch}"] = float((correct[branch] * scale.scale).std())

    # --- numerical floor: identical protein, and a repeated forward pass ---
    repeat = branches(model, protein_inputs(data, spec.target, device, dtype),
                      *query, device, dtype)
    row["floor_identical_protein_centered_pk"] = float(np.abs(
        centered(correct["interaction"] - repeat["interaction"])).max()
        * scale.scale)
    row["floor_identical_protein_level_pk"] = float(abs(
        (correct["protein_only"] - repeat["protein_only"]).mean()) * scale.scale)
    again = branches(model, correct_parts, *query, device, dtype)
    row["floor_repeated_forward_centered_pk"] = float(np.abs(
        centered(correct["interaction"] - again["interaction"])).max()
        * scale.scale)

    # --- wrong protein, one measurement per similarity stratum -------------
    truth_centered = centered(truth)
    for stratum, (donor_target, similarity) in donors[spec.target].items():
        wrong = branches(model, protein_inputs(data, donor_target, device, dtype),
                         *query, device, dtype)
        shift = correct["interaction"] - wrong["interaction"]
        row[f"donor_similarity__{stratum}"] = similarity
        row[f"r_interaction_wrong__{stratum}"] = correlation(
            wrong["interaction"], truth)
        row[f"r_full_wrong__{stratum}"] = correlation(wrong["full"], truth)
        row[f"level_shift_pk__{stratum}"] = float(abs(
            (correct["protein_only"] - wrong["protein_only"]).mean())
            * scale.scale)
        row[f"interaction_centered_shift_pk__{stratum}"] = float(
            centered(shift).std() * scale.scale)
        row[f"interaction_uncentered_shift_pk__{stratum}"] = float(
            abs(shift.mean()) * scale.scale)
        # Is the protein-driven change *aligned* with real affinity differences,
        # or is it arbitrary movement? A shift orthogonal to truth cannot help.
        row[f"shift_alignment__{stratum}"] = correlation(shift, truth_centered)
        row[f"shift_vector__{stratum}"] = (
            centered(shift) * scale.scale).astype(float).tolist()

    # --- shuffled protein: same composition, scrambled arrangement ---------
    scrambled = branches(model, scramble_slots(correct_parts, control_seed),
                         *query, device, dtype)
    row["level_shift_pk__scrambled"] = float(abs(
        (correct["protein_only"] - scrambled["protein_only"]).mean())
        * scale.scale)
    row["interaction_centered_shift_pk__scrambled"] = float(
        centered(correct["interaction"] - scrambled["interaction"]).std()
        * scale.scale)
    row["r_interaction_wrong__scrambled"] = correlation(
        scrambled["interaction"], truth)

    # --- shuffled labels: what the metric reports on pure noise ------------
    permutation = np.random.default_rng(
        control_seed + 1).permutation(len(truth))
    row["r_full__shuffled_label"] = correlation(correct["full"], truth[permutation])
    row["r_interaction__shuffled_label"] = correlation(
        correct["interaction"], truth[permutation])

    # --- foreign ligand panel: the recipient protein, someone else's ligands
    if foreign_episode is not None:
        foreign = query_tensors(foreign_episode, device, dtype)
        foreign_correct = branches(model, correct_parts, *foreign, device, dtype)
        foreign_wrong = branches(
            model,
            protein_inputs(data, donors[spec.target]["nearest"][0], device, dtype),
            *foreign, device, dtype)
        # The foreign panel may be a different size; compare on the overlap.
        # These ligands do not correspond to these labels, so the correlation
        # must collapse — that is the point of the control.
        width = min(len(truth), len(foreign_correct["full"]))
        row["r_full__foreign_ligands"] = correlation(
            foreign_correct["full"][:width], truth[:width])
        row["interaction_centered_shift_pk__foreign_ligands"] = float(
            centered(foreign_correct["interaction"]
                     - foreign_wrong["interaction"]).std() * scale.scale)
    return row


# --------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--random-inits", type=int,
                        default=len(_frozen.RANDOM_INIT_SEEDS))
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=_frozen.SPLIT_DIRECTORY)
    scale = training_label_scale(data)
    donors = stratified_donors(data, "meta_val", _frozen.DONOR_POOL,
                               _frozen.WHITENING_POOL)
    specs = data.fixed_nested_episode_banks(
        "meta_val", _frozen.SUPPORT_SIZES, _frozen.QUERY_SIZE, 1,
        _frozen.EVALUATION_SEED, None)[0]

    scaffolds = murcko_scaffolds(data._ligand_smiles)
    novelty, scaffold_seen = novelty_and_scaffold_strata(data, scaffolds)

    # Foreign panels are drawn once, before any model runs, so every arm sees
    # the identical control input.
    generator = np.random.default_rng(_frozen.CONTROL_SEED)
    foreign_of = {}
    for spec in specs:
        others = [s for s in specs if s.target != spec.target
                  and len(s.query) >= len(spec.query)]
        foreign_of[spec.target] = (
            others[int(generator.integers(len(others)))] if others else None)

    episodes = {}
    for spec in specs:
        episodes[spec.target] = compact_episode(
            normalized(data.materialize(spec), scale))
    foreign_episodes = {
        target: (episodes[other.target] if other is not None else None)
        for target, other in foreign_of.items()}

    arms: list[tuple[str, int, object]] = []
    for path in _frozen.A0_CHECKPOINTS:
        arms.append(("A0", int(str(path.parent.name).rsplit("seed", 1)[-1]), path))
    reference = _frozen.A0_CHECKPOINTS[0]
    for seed in _frozen.RANDOM_INIT_SEEDS[:arguments.random_inits]:
        arms.append(("randinit", int(seed), reference))

    rows: list[dict] = []
    for name, seed, path in arms:
        model, _, _ = (trained_arm(path, data, arguments.device) if name == "A0"
                       else random_arm(path, data, arguments.device, seed))
        model.eval()
        dtype = next(model.parameters()).dtype
        with torch.no_grad():
            for spec in specs:
                row = measure(model, data, spec, episodes[spec.target], donors,
                              scale, arguments.device, dtype,
                              _frozen.CONTROL_SEED,
                              foreign_episodes[spec.target])
                ligands = [data.cells[i]["ligand_id"] for i in spec.query]
                row["mean_novelty"] = float(np.mean([novelty(l) for l in ligands]))
                row["seen_scaffold_fraction"] = float(
                    np.mean([scaffold_seen(l) for l in ligands]))
                rows.append(dict(row, arm=name, seed=seed))
        del model
        if arguments.device.startswith("cuda"):
            torch.cuda.empty_cache()
        print(f"  measured {name} seed {seed}")

    payload = summarise(rows, data)
    payload["frozen_design"] = _frozen.frozen_manifest()
    payload["meta_test"] = data.seal_record()
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    with arguments.output.with_suffix(".rows.jsonl").open(
            "w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(
                {k: v for k, v in row.items()
                 if not k.startswith("shift_vector__")}) + "\n")
    report(payload)
    print(f"\nwrote {arguments.output}")
    return 0


def summarise(rows: list[dict], data) -> dict:
    """Aggregate seeds within target, then bootstrap over components."""
    def weighted(subset, field):
        return component_target_mean(
            (r["component"], r["target"], r.get(field)) for r in subset)

    def paired(subset, left, right):
        return component_bootstrap(
            [(r["component"], r["target"], r[left] - r[right]) for r in subset
             if r.get(left) is not None and r.get(right) is not None],
            _frozen.BOOTSTRAP_DRAWS, _frozen.BOOTSTRAP_SEED)

    payload: dict = {"schema": "MetaSieve.A2ReadinessV2.BranchOrdering.v1",
                     "split": "meta_val", "k": 0, "arms": {}}
    for arm in ("A0", "randinit"):
        subset = [r for r in rows if r["arm"] == arm]
        if not subset:
            continue
        cell: dict = {"seeds": sorted({r["seed"] for r in subset}),
                      "targets": len({r["target"] for r in subset}),
                      "components": len({r["component"] for r in subset})}
        for field in ("r_full", "r_ligand_only", "r_interaction",
                      "ci_full", "ci_interaction",
                      "sd_ligand_only", "sd_interaction", "sd_truth",
                      "floor_identical_protein_centered_pk",
                      "floor_identical_protein_level_pk",
                      "floor_repeated_forward_centered_pk",
                      "level_shift_pk__scrambled",
                      "interaction_centered_shift_pk__scrambled",
                      "r_full__shuffled_label", "r_interaction__shuffled_label",
                      "r_full__foreign_ligands",
                      "interaction_centered_shift_pk__foreign_ligands",
                      "mean_novelty", "seen_scaffold_fraction"):
            cell[field] = weighted(subset, field)

        cell["interaction_ordering_increment"] = paired(
            subset, "r_full", "r_ligand_only")
        cell["strata"] = {}
        for stratum in _frozen.DONOR_STRATA:
            contrast = paired(subset, "r_interaction",
                              f"r_interaction_wrong__{stratum}")
            cell["strata"][stratum] = {
                "donor_similarity": weighted(subset, f"donor_similarity__{stratum}"),
                "level_shift_pk": weighted(subset, f"level_shift_pk__{stratum}"),
                "interaction_centered_shift_pk": weighted(
                    subset, f"interaction_centered_shift_pk__{stratum}"),
                "interaction_uncentered_shift_pk": weighted(
                    subset, f"interaction_uncentered_shift_pk__{stratum}"),
                "shift_alignment_with_truth": weighted(
                    subset, f"shift_alignment__{stratum}"),
                "r_interaction_minus_wrong": contrast,
                "verdict": _frozen.verdict(contrast),
            }

        # Variance decomposition: how much of the spread in `r_interaction` is
        # between components, between targets within a component, and between
        # seeds on the same target.
        cell["variability"] = variance_components(subset, "r_interaction")
        cell["variability_shift"] = variance_components(
            subset, "interaction_centered_shift_pk__nearest")
        payload["arms"][arm] = cell

    payload["random_init_reproducibility"] = reproducibility(
        [r for r in rows if r["arm"] == "randinit"])
    payload["novelty_strata"] = novelty_strata(rows)
    return payload


def variance_components(subset: list[dict], field: str) -> dict:
    """Between-component / between-target / between-seed spread, separately."""
    by_target: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in subset:
        value = row.get(field)
        if value is not None and np.isfinite(value):
            by_target[(row["component"], row["target"])].append(float(value))
    if not by_target:
        return {}
    seed_sd = [float(np.std(v)) for v in by_target.values() if len(v) > 1]
    target_means = {key: float(np.mean(v)) for key, v in by_target.items()}
    by_component: dict[str, list[float]] = defaultdict(list)
    for (component, _), value in target_means.items():
        by_component[component].append(value)
    component_means = [float(np.mean(v)) for v in by_component.values()]
    within_component = [float(np.std(v)) for v in by_component.values() if len(v) > 1]
    return {
        "between_component_sd": float(np.std(component_means)),
        "between_target_within_component_sd": (
            float(np.mean(within_component)) if within_component else 0.0),
        "between_seed_within_target_sd": (
            float(np.mean(seed_sd)) if seed_sd else 0.0),
        "n_components": len(by_component), "n_targets": len(by_target),
    }


def reproducibility(rows: list[dict]) -> dict:
    """Do two independent initialisations move the *same way* on a target?

    A nonzero shift magnitude at random init says only that the untrained map is
    not constant in the protein. If the shift direction is uncorrelated across
    initialisations, it is a property of the draw, not of the architecture, and
    it cannot be evidence that training destroyed a structural capability.
    """
    if not rows:
        return {}
    out: dict = {}
    for stratum in ("nearest", "farthest"):
        key = f"shift_vector__{stratum}"
        by_target: dict[str, list[np.ndarray]] = defaultdict(list)
        for row in rows:
            vector = row.get(key)
            if vector:
                by_target[row["target"]].append(np.asarray(vector, dtype=float))
        pairwise, magnitudes = [], []
        for vectors in by_target.values():
            magnitudes.extend(float(np.std(v)) for v in vectors)
            for i in range(len(vectors)):
                for j in range(i + 1, len(vectors)):
                    a, b = vectors[i], vectors[j]
                    denominator = float(np.linalg.norm(a) * np.linalg.norm(b))
                    if denominator > 1e-12:
                        pairwise.append(float(a @ b / denominator))
        out[stratum] = {
            "mean_pairwise_shift_cosine": (
                float(np.mean(pairwise)) if pairwise else float("nan")),
            "sd_pairwise_shift_cosine": (
                float(np.std(pairwise)) if pairwise else float("nan")),
            "pairs": len(pairwise),
            "mean_shift_magnitude_pk": (
                float(np.mean(magnitudes)) if magnitudes else float("nan")),
            "note": ("cosine ~0 means each initialisation moves in its own "
                     "arbitrary direction; the sensitivity is a draw property, "
                     "not an architectural capability"),
        }
    return out


def novelty_strata(rows: list[dict]) -> dict:
    """Is any effect concentrated in the least novel ligands?"""
    out: dict = {}
    for arm in ("A0", "randinit"):
        subset = [r for r in rows if r["arm"] == arm]
        if not subset:
            continue
        values = sorted(r["mean_novelty"] for r in subset)
        if not values:
            continue
        low, high = np.quantile(values, [1 / 3, 2 / 3])
        buckets = {"low_novelty_most_similar_to_train": lambda v: v >= high,
                   "mid": lambda v: low <= v < high,
                   "high_novelty_least_similar": lambda v: v < low}
        out[arm] = {}
        for label, test in buckets.items():
            block = [r for r in subset if test(r["mean_novelty"])]
            if not block:
                continue
            out[arm][label] = {
                "targets": len({r["target"] for r in block}),
                "mean_max_tanimoto_to_meta_train": float(
                    np.mean([r["mean_novelty"] for r in block])),
                "seen_scaffold_fraction": float(
                    np.mean([r["seen_scaffold_fraction"] for r in block])),
                "r_interaction": component_target_mean(
                    (r["component"], r["target"], r["r_interaction"])
                    for r in block),
                "interaction_centered_shift_pk_nearest": component_target_mean(
                    (r["component"], r["target"],
                     r["interaction_centered_shift_pk__nearest"])
                    for r in block),
            }
    return out


def report(payload: dict) -> None:
    for arm, cell in payload["arms"].items():
        print(f"\n=== {arm}  ({len(cell['seeds'])} seeds, {cell['targets']} "
              f"targets, {cell['components']} components)")
        print(f"  r_full {cell['r_full']:+.4f} | r_ligand_only "
              f"{cell['r_ligand_only']:+.4f} | r_interaction "
              f"{cell['r_interaction']:+.4f}")
        increment = cell["interaction_ordering_increment"]
        print(f"  r(full)-r(ligand_only)  {increment['mean']:+.4f} "
              f"[{increment['lo']:+.4f},{increment['hi']:+.4f}] "
              f"{_frozen.verdict(increment)}")
        print(f"  floors: identical-protein centered "
              f"{cell['floor_identical_protein_centered_pk']:.2e} pK | "
              f"repeated forward "
              f"{cell['floor_repeated_forward_centered_pk']:.2e} pK")
        print(f"  shuffled-label r_full {cell['r_full__shuffled_label']:+.4f} "
              f"| scrambled-protein level shift "
              f"{cell['level_shift_pk__scrambled']:.4f} pK")
        print("  donor strata:")
        print(f"    {'stratum':<10}{'sim':>7}{'level':>9}{'centered':>10}"
              f"{'align':>8}{'r-diff':>9}   verdict")
        for name, block in cell["strata"].items():
            contrast = block["r_interaction_minus_wrong"]
            print(f"    {name:<10}{block['donor_similarity']:>7.3f}"
                  f"{block['level_shift_pk']:>9.4f}"
                  f"{block['interaction_centered_shift_pk']:>10.4f}"
                  f"{block['shift_alignment_with_truth']:>8.3f}"
                  f"{contrast['mean']:>+9.4f}   {block['verdict']}")
        variability = cell["variability"]
        print(f"  r_interaction spread: between-component "
              f"{variability['between_component_sd']:.4f} | between-target "
              f"{variability['between_target_within_component_sd']:.4f} | "
              f"between-seed {variability['between_seed_within_target_sd']:.4f}")

    print("\n=== random-init reproducibility (is the shift a capability or a draw?)")
    for stratum, block in payload.get("random_init_reproducibility", {}).items():
        print(f"  {stratum}: mean pairwise shift cosine "
              f"{block['mean_pairwise_shift_cosine']:+.4f} "
              f"(sd {block['sd_pairwise_shift_cosine']:.4f}, "
              f"{block['pairs']} pairs), mean |shift| "
              f"{block['mean_shift_magnitude_pk']:.4f} pK")


if __name__ == "__main__":
    raise SystemExit(main())

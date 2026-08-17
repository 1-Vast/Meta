"""Phase 2: does an *internal* A0 representation carry protein-conditioned SAR?

Two parts, and the second only matters if the first says there is something to
find.

**Part A, no training.** For every internal representation, how much of its
*ligand-differential* — the part left after removing the within-target mean —
changes when the recipient protein is replaced by a matched wrong one? The v1
probe asked this of the final scalar only. If an internal representation is
protein-conditioned while the scalar is not, A2 should be amended rather than
rejected, and the target is the readout.

**Part B, small trained probes.** Existence of a difference is not existence of
usable signal. A protein-conditioned representation earns the A2 premise only
if a *transferable* readout of it predicts within-target affinity differences
on held-out protein components. The probe is deliberately the smallest thing
that could work and is trained by ordinary gradient descent — no closed form,
no pseudoinverse, no inner loop:

    delta_hat(P, L_i, L_j) = s * < g(P), U (e_i - e_j) >

`U: R^D -> R^R` is the ligand-side projection A2-min calls `A_phi`; `g` maps the
pooled protein to the same R-dimensional space. The decisive contrast is
against an identical model whose `g(P)` is replaced by a learned **constant**:
if conditioning on the protein buys nothing, there is no protein-conditioned
SAR coordinate to identify, and A2's premise fails no matter how the operator
is parameterised.

Arms, all preregistered in `_frozen.py` before any result was seen:

| arm | what changes |
|---|---|
| `protein_conditioned` | the candidate |
| `shared_direction` | `g(P)` -> learned constant (no protein conditioning) |
| `wrong_protein` | evaluated with the donor's protein in `g` |
| `reference_protein` | features from one fixed protein for every target |
| `label_shuffled` | labels permuted within target, refit from scratch |
| `random_feature` | `e` replaced by frozen Gaussian noise of the same width |
| `random_projection` | `U` frozen at random init, only `g` and `s` trained |

Hyperparameters and early stopping use `meta_train` component folds only.
`meta_val` is scored once, with the probe frozen. `meta_test` is unreachable.
Query labels are targets of the loss and nothing else.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import torch
import torch.nn as nn

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.stageR0_retrieval_falsification import (                # noqa: E402
    component_bootstrap, component_target_mean,
)
from tools.research.a2_readiness_v2 import _frozen                   # noqa: E402
from tools.research.a2_readiness_v2._features import (               # noqa: E402
    PROTEIN_CONDITIONED, differential_geometry,
)

# Preregistered: the A2 plan's `e0`, plus the earliest pooled representation,
# so the result cannot be a consequence of the `embed` bottleneck alone. Both
# are named in advance; this is two arms, not a search over six.
#
# `ligand` was added **after** the first run, and is labelled as such. It is
# the protein-blind encoder output, so it cannot become a new candidate — it
# can only take credit away from `embed`. The first run found a resolved
# *protein-independent* SAR direction in `embed`; `embed` is built from
# `cat(ligand, mean_state, max_state, wide_summary, occupancy)`, so that
# direction might live entirely in its ligand component and add nothing to what
# the ligand encoder already provides. A post-hoc control that can only shrink
# one's own positive finding is the safe kind to add; it is reported as a
# control, never as a selected representation.
#
# `max_state` was added after Phase 3, which measured it as the representation
# that retains the most protein-differential of any stage (relative change
# 0.0503 under a protein swap, 16x `mean_state`'s 0.0031). It is therefore the
# strongest surviving candidate for A2, and testing it can only *help* A2 —
# adding it is the conservative direction for a post-hoc arm, exactly as
# `ligand` was the conservative direction against my own positive finding.
PROBE_REPRESENTATIONS = ("embed", "mean_state", "ligand", "max_state")

PROBE_RANKS = (4, 8, 16)          # the A2 plan's own search range
PROBE_STEPS = 400
PROBE_LR = 3e-3
PROBE_FOLDS = 5


# ---------------------------------------------------------------------------
# Part A — no training
# ---------------------------------------------------------------------------

def part_a(features: dict, targets: np.ndarray, components: np.ndarray) -> dict:
    """Per-representation ligand-differential geometry under a protein swap."""
    out: dict = {}
    for name in PROTEIN_CONDITIONED + ("ligand",):
        correct = features[f"correct__{name}"]
        rows = []
        for condition in ("wrong", "reference"):
            other = features[f"{condition}__{name}"]
            per_target = []
            for index in np.unique(targets):
                select = targets == index
                if select.sum() < 2:
                    continue
                geometry = differential_geometry(correct[select], other[select])
                geometry["component"] = int(components[index])
                per_target.append(geometry)
            block = {
                field: component_target_mean(
                    (row["component"], i, row[field])
                    for i, row in enumerate(per_target))
                for field in ("level_relative_shift", "differential_cosine",
                              "differential_relative_shift",
                              "differential_norm_ratio")}
            block["targets"] = len(per_target)
            block["differential_cosine_ci"] = component_bootstrap(
                [(row["component"], i, row["differential_cosine"])
                 for i, row in enumerate(per_target)],
                _frozen.BOOTSTRAP_DRAWS, _frozen.BOOTSTRAP_SEED)
            rows.append((condition, block))
        out[name] = dict(rows)
        out[name]["width"] = int(correct.shape[1])
    return out


# ---------------------------------------------------------------------------
# Part B — small trained probes
# ---------------------------------------------------------------------------

class DeltaProbe(nn.Module):
    """`delta_hat = s * <g(P), U (e_i - e_j)>`, trained end to end by SGD."""

    def __init__(self, width: int, protein_width: int, rank: int,
                 condition_on_protein: bool, freeze_projection: bool = False):
        super().__init__()
        self.projection = nn.Linear(width, rank, bias=False)
        self.condition_on_protein = bool(condition_on_protein)
        if condition_on_protein:
            self.protein = nn.Sequential(
                nn.Linear(protein_width, 64), nn.GELU(), nn.Linear(64, rank))
        else:
            # The control: the same capacity on the ligand side, but the
            # protein-side vector is a single learned constant. Anything this
            # arm achieves is transferable SAR that does not need the protein.
            self.constant = nn.Parameter(torch.zeros(rank))
            nn.init.normal_(self.constant, std=0.1)
        self.gain = nn.Parameter(torch.tensor(1.0))
        if freeze_projection:
            self.projection.weight.requires_grad_(False)

    def forward(self, delta_features: torch.Tensor,
                protein: torch.Tensor) -> torch.Tensor:
        ligand = self.projection(delta_features)
        side = (self.protein(protein) if self.condition_on_protein
                else self.constant.expand_as(ligand))
        return self.gain * (side * ligand).sum(-1)


def build_pairs(features: np.ndarray, labels: np.ndarray,
                targets: np.ndarray, rng: np.random.Generator,
                pairs_per_target: int = 48):
    """Within-target ligand pairs. Never crosses a target boundary."""
    left, right, delta, owner = [], [], [], []
    for index in np.unique(targets):
        select = np.flatnonzero(targets == index)
        if len(select) < 2:
            continue
        count = min(pairs_per_target, len(select) * (len(select) - 1) // 2)
        a = rng.choice(select, size=count * 2, replace=True)
        b = rng.choice(select, size=count * 2, replace=True)
        keep = a != b
        a, b = a[keep][:count], b[keep][:count]
        left.append(features[a])
        right.append(features[b])
        delta.append(labels[a] - labels[b])
        owner.append(np.full(len(a), index, dtype=np.int64))
    return (np.concatenate(left), np.concatenate(right),
            np.concatenate(delta), np.concatenate(owner))


def fit(probe: DeltaProbe, delta_features: torch.Tensor,
        protein: torch.Tensor, delta_y: torch.Tensor, steps: int,
        learning_rate: float, device: str) -> DeltaProbe:
    optimiser = torch.optim.Adam(
        [p for p in probe.parameters() if p.requires_grad], lr=learning_rate)
    probe.train()
    for _ in range(steps):
        optimiser.zero_grad(set_to_none=True)
        loss = (probe(delta_features, protein) - delta_y).square().mean()
        loss.backward()
        optimiser.step()
    probe.eval()
    return probe


def score(probe: DeltaProbe, delta_features: torch.Tensor,
          protein: torch.Tensor, delta_y: torch.Tensor,
          owner: np.ndarray, component_of_target: np.ndarray) -> dict:
    """Pair-level correlation, aggregated per target then per component."""
    with torch.no_grad():
        prediction = probe(delta_features, protein).cpu().numpy()
    truth = delta_y.cpu().numpy()
    rows = []
    for index in np.unique(owner):
        select = owner == index
        p, y = prediction[select], truth[select]
        if len(p) < 3 or p.std() < 1e-9 or y.std() < 1e-9:
            continue
        rows.append((int(component_of_target[index]), int(index),
                     float(np.corrcoef(p, y)[0, 1])))
    return {
        "delta_r": component_target_mean(rows),
        "delta_r_ci": component_bootstrap(
            rows, _frozen.BOOTSTRAP_DRAWS, _frozen.BOOTSTRAP_SEED),
        "targets": len(rows),
        "pair_level_r": float(np.corrcoef(prediction, truth)[0, 1])
        if prediction.std() > 1e-9 else 0.0,
        "sign_accuracy": float((np.sign(prediction) == np.sign(truth)).mean()),
    }


def tensors(left, right, protein_rows, delta, device):
    return (torch.as_tensor(left - right, dtype=torch.float32, device=device),
            torch.as_tensor(protein_rows, dtype=torch.float32, device=device),
            torch.as_tensor(delta, dtype=torch.float32, device=device))


def part_b(train: dict, val: dict, representation: str, device: str,
           rng_seed: int) -> dict:
    """Select the rank on meta_train folds, refit, score meta_val once."""
    train_features = train[f"correct__{representation}"]
    train_labels, train_targets = train["y"], train["target_index"]
    train_components = train["component_of_target"]
    train_protein = train["protein_pooled"]

    # Normalisation from meta_train only, and from the *features*, never labels.
    center = train_features.mean(0, keepdims=True)
    spread = train_features.std(0, keepdims=True) + 1e-6
    protein_center = train_protein.mean(0, keepdims=True)
    protein_spread = train_protein.std(0, keepdims=True) + 1e-6

    def prepare(block: dict, feature_key: str, protein_key: str,
                labels: np.ndarray | None = None, seed: int = 0):
        """Pairs, their label differences, and the protein row for each pair.

        `protein_key` selects *which* protein conditions the probe, so a
        wrong-protein arm can replace the protein on both sides — the features
        and the conditioning — instead of leaving one of them correct.
        """
        raw = (block[feature_key] - center) / spread
        y = block["y"] if labels is None else labels
        rng = np.random.default_rng(seed)
        left, right, delta, owner = build_pairs(
            raw, y, block["target_index"], rng)
        protein = (block[protein_key] - protein_center) / protein_spread
        return left, right, delta, owner, protein[owner]

    width = train_features.shape[1]
    protein_width = train_protein.shape[1]

    # --- rank selection on meta_train component folds --------------------
    # Selected *separately* for the protein-conditioned and the constant
    # configurations. Giving both the same rank would hand one of them a
    # hyperparameter chosen for the other, and the comparison between them is
    # the whole result.
    folds = PROBE_FOLDS
    unique_components = np.unique(train_components)
    assignment = {int(c): int(i % folds)
                  for i, c in enumerate(np.random.default_rng(
                      _frozen.MODEL_FOLD_SEED).permutation(unique_components))}
    fold_of_row = np.asarray(
        [assignment[int(train_components[t])] for t in train_targets])

    def fold_score(rank: int, condition: bool) -> float:
        held = []
        for fold in range(folds):
            fit_rows = fold_of_row != fold
            score_rows = ~fit_rows
            if score_rows.sum() < 40 or fit_rows.sum() < 40:
                continue
            block_fit = {k: (v[fit_rows] if k in ("y", "target_index")
                             or k.endswith(representation) else v)
                         for k, v in train.items()}
            block_score = {k: (v[score_rows] if k in ("y", "target_index")
                               or k.endswith(representation) else v)
                           for k, v in train.items()}
            left, right, delta, owner, protein = prepare(
                block_fit, f"correct__{representation}", "protein_pooled",
                seed=rng_seed + fold)
            features_t, protein_t, delta_t = tensors(
                left, right, protein, delta, device)
            probe = DeltaProbe(width, protein_width, rank, condition).to(device)
            fit(probe, features_t, protein_t, delta_t, PROBE_STEPS, PROBE_LR,
                device)
            left, right, delta, owner, protein = prepare(
                block_score, f"correct__{representation}", "protein_pooled",
                seed=rng_seed + 100 + fold)
            features_t, protein_t, delta_t = tensors(
                left, right, protein, delta, device)
            held.append(score(probe, features_t, protein_t, delta_t, owner,
                              train_components)["delta_r"])
        return float(np.mean(held)) if held else float("nan")

    selection = {
        "protein_conditioned": {rank: fold_score(rank, True)
                                for rank in PROBE_RANKS},
        "shared_direction": {rank: fold_score(rank, False)
                             for rank in PROBE_RANKS},
    }

    def best_of(curve: dict) -> int:
        return max(PROBE_RANKS, key=lambda r: (curve[r] if np.isfinite(curve[r])
                                               else -np.inf))

    best_rank = best_of(selection["protein_conditioned"])
    shared_rank = best_of(selection["shared_direction"])

    # --- refit on all of meta_train, score meta_val once ------------------
    results: dict = {"rank_selection_on_meta_train_folds": selection,
                     "selected_rank": best_rank,
                     "selected_rank_shared_direction": shared_rank,
                     "width": int(width), "arms": {}}

    shuffled_labels = train_labels.copy()
    shuffle_rng = np.random.default_rng(_frozen.CONTROL_SEED)
    for index in np.unique(train_targets):
        select = np.flatnonzero(train_targets == index)
        shuffled_labels[select] = shuffled_labels[
            shuffle_rng.permutation(select)]

    noise_rng = np.random.default_rng(_frozen.CONTROL_SEED + 7)
    random_train = noise_rng.normal(size=train_features.shape).astype(np.float32)
    random_val = noise_rng.normal(
        size=val[f"correct__{representation}"].shape).astype(np.float32)

    # Capacity-matched null. `g` has ~42k parameters against the constant's
    # `rank`, so a protein-conditioned arm that loses could be losing to
    # overfitting rather than to the absence of protein signal. This arm keeps
    # the identical architecture, parameter count and optimiser and only
    # destroys the *correspondence* between protein and target — during
    # training and evaluation alike. If the real protein does no better than a
    # permuted one, the protein input is carrying nothing.
    permute_rng = np.random.default_rng(_frozen.CONTROL_SEED + 11)
    train_protein_permuted = train["protein_pooled"][
        permute_rng.permutation(len(train["protein_pooled"]))]
    val_protein_permuted = val["protein_pooled"][
        permute_rng.permutation(len(val["protein_pooled"]))]

    arms = {
        "protein_conditioned": dict(condition=True, frozen=False,
                                    train_key=f"correct__{representation}",
                                    val_key=f"correct__{representation}",
                                    val_protein="protein_pooled", labels=None),
        "shared_direction": dict(condition=False, frozen=False,
                                 train_key=f"correct__{representation}",
                                 val_key=f"correct__{representation}",
                                 val_protein="protein_pooled", labels=None),
        "wrong_protein": dict(condition=True, frozen=False,
                              train_key=f"correct__{representation}",
                              val_key=f"wrong__{representation}",
                              val_protein="wrong_protein_pooled", labels=None),
        "reference_protein": dict(condition=True, frozen=False,
                                  train_key=f"correct__{representation}",
                                  val_key=f"reference__{representation}",
                                  val_protein="reference_protein_pooled",
                                  labels=None),
        "label_shuffled": dict(condition=True, frozen=False,
                               train_key=f"correct__{representation}",
                               val_key=f"correct__{representation}",
                               val_protein="protein_pooled",
                               labels=shuffled_labels),
        "random_feature": dict(condition=True, frozen=False,
                               train_key="__random__", val_key="__random__",
                               val_protein="protein_pooled", labels=None),
        "random_projection": dict(condition=True, frozen=True,
                                  train_key=f"correct__{representation}",
                                  val_key=f"correct__{representation}",
                                  val_protein="protein_pooled", labels=None),
        "protein_permuted": dict(condition=True, frozen=False,
                                 train_key=f"correct__{representation}",
                                 val_key=f"correct__{representation}",
                                 val_protein="protein_pooled", labels=None,
                                 permute_protein=True),
    }

    for name, spec in arms.items():
        train_block = dict(train)
        val_block = dict(val)
        if spec.get("permute_protein"):
            train_block["protein_pooled"] = train_protein_permuted
            val_block["protein_pooled"] = val_protein_permuted
        if spec["train_key"] == "__random__":
            train_block["__random__"] = random_train
            val_block["__random__"] = random_val
            feature_key_train = feature_key_val = "__random__"
        else:
            feature_key_train = spec["train_key"]
            feature_key_val = spec["val_key"]

        # Every arm is *trained* on the correct protein: the question is what a
        # readout learned under correct conditions does when the condition is
        # replaced at evaluation. Training each arm on its own corruption would
        # let it re-fit around the corruption and hide the dependence.
        left, right, delta, owner, protein = prepare(
            train_block, feature_key_train, "protein_pooled",
            labels=spec["labels"], seed=rng_seed)
        features_t, protein_t, delta_t = tensors(
            left, right, protein, delta, device)
        rank = best_rank if spec["condition"] else shared_rank
        probe = DeltaProbe(width, protein_width, rank,
                           spec["condition"], spec["frozen"]).to(device)
        fit(probe, features_t, protein_t, delta_t, PROBE_STEPS, PROBE_LR, device)

        left, right, delta, owner, protein = prepare(
            val_block, feature_key_val, spec["val_protein"], seed=rng_seed + 1)
        features_t, protein_t, delta_t = tensors(
            left, right, protein, delta, device)
        results["arms"][name] = dict(
            score(probe, features_t, protein_t, delta_t, owner,
                  val["component_of_target"]),
            rank=rank,
            trainable_parameters=int(sum(
                p.numel() for p in probe.parameters() if p.requires_grad)))

    # The contrast the A2 premise stands or falls on, as a paired interval
    # rather than two independent ones.
    results["protein_conditioning_gain"] = {
        "vs_shared_direction": (results["arms"]["protein_conditioned"]["delta_r"]
                                - results["arms"]["shared_direction"]["delta_r"]),
        "vs_permuted_protein": (results["arms"]["protein_conditioned"]["delta_r"]
                                - results["arms"]["protein_permuted"]["delta_r"]),
        "note": ("positive means knowing the recipient protein buys "
                 "transferable within-target SAR ordering that a "
                 "protein-independent direction does not already provide"),
    }
    return results


# ---------------------------------------------------------------------------

def load(path: Path) -> dict:
    with np.load(path) as stored:
        return {key: stored[key] for key in stored.files}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    val = load(arguments.features / "meta_val.npz")
    train = load(arguments.features / "meta_train.npz")

    payload: dict = {
        "schema": "MetaSieve.A2ReadinessV2.RepresentationProbe.v1",
        "frozen_design": _frozen.frozen_manifest(),
        "probe": {"representations": list(PROBE_REPRESENTATIONS),
                  "ranks": list(PROBE_RANKS), "steps": PROBE_STEPS,
                  "learning_rate": PROBE_LR, "folds": PROBE_FOLDS,
                  "optimiser": "Adam, ordinary gradient descent; no closed "
                               "form, no pseudoinverse, no inner loop"},
        "part_a_differential_geometry": part_a(
            val, val["target_index"], val["component_of_target"]),
        "part_b_trained_probes": {},
    }
    for representation in PROBE_REPRESENTATIONS:
        print(f"\n--- part B: {representation}")
        payload["part_b_trained_probes"][representation] = part_b(
            train, val, representation, arguments.device, _frozen.CONTROL_SEED)

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    report(payload)
    print(f"\nwrote {arguments.output}")
    return 0


def report(payload: dict) -> None:
    print("\n=== Part A: ligand-differential geometry under a protein swap")
    print(f"{'representation':<14}{'w':>5}{'lvl shift':>11}{'diff cos':>10}"
          f"{'diff shift':>12}{'[95% CI on cos]':>24}")
    for name, block in payload["part_a_differential_geometry"].items():
        wrong = block["wrong"]
        interval = wrong["differential_cosine_ci"]
        print(f"{name:<14}{block['width']:>5}"
              f"{wrong['level_relative_shift']:>11.4f}"
              f"{wrong['differential_cosine']:>10.4f}"
              f"{wrong['differential_relative_shift']:>12.4f}"
              f"   [{interval['lo']:+.4f},{interval['hi']:+.4f}]")

    print("\n=== Part B: transferable delta-affinity readout (meta_val, once)")
    for representation, block in payload["part_b_trained_probes"].items():
        print(f"\n  {representation} (width {block['width']}); rank selected on "
              f"meta_train folds: protein-conditioned {block['selected_rank']}, "
              f"shared-direction {block['selected_rank_shared_direction']}")
        for configuration, curve in block[
                "rank_selection_on_meta_train_folds"].items():
            print("      fold curve "
                  f"{configuration:<20}{ {k: round(v, 4) for k, v in curve.items()} }")
        print(f"    {'arm':<22}{'rank':>5}{'params':>9}{'delta_r':>9}"
              f"{'[95% CI]':>22}{'sign':>8}")
        for arm, cell in block["arms"].items():
            interval = cell["delta_r_ci"]
            print(f"    {arm:<22}{cell['rank']:>5}{cell['trainable_parameters']:>9}"
                  f"{cell['delta_r']:>+9.4f}"
                  f"   [{interval['lo']:+.4f},{interval['hi']:+.4f}]"
                  f"{cell['sign_accuracy']:>8.3f}")
        gain = block["protein_conditioning_gain"]
        print(f"    protein conditioning gain: vs shared direction "
              f"{gain['vs_shared_direction']:+.4f} | vs permuted protein "
              f"{gain['vs_permuted_protein']:+.4f}")


if __name__ == "__main__":
    raise SystemExit(main())

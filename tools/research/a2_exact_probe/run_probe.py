"""Stage R: the exact episodic A2 operator on real episodes, k=0/1/2/3/5.

Trains `A_φ` and the two shrinkage scalars by ordinary gradient descent on
`meta_train` components only, selects the rank on `meta_train` component folds,
and scores `meta_val` **once** after the choice is frozen. The frozen A0 trunk
supplies `e0` and `f0`; nothing in the trunk trains.

The objective is the episodic one A2 actually proposes:

    minimise  mean_q ( f0(P,L_q) + δ_q  −  y_q )²
    where     δ_q = η(k) ⟨ (1/k) Σ_i r_i z_i , z_q ⟩,  r_i = stopgrad(y_i − f0_i)

Because `f0` is frozen this is exactly a regression of `δ_q` onto the query
residual `y_q − f0(P,L_q)`, which is the quantity A2 exists to predict.

Metrics are reported in pK on the *full* prediction `f0 + δ`, so they are
comparable across arms and across k:

* `mse_pk` — the primary quantity;
* `r` — within-target Pearson correlation, the R14 ordering lever;
* `ci` — concordance;
* `spearman`;
* `query_spread_pk` — the standard deviation of `δ` across the queries of one
  episode. **At k=1 this is the gate that separates A2 from A0**, whose k=1
  transport is provably a pure level shift.

Query labels appear only as the loss target and in the metrics. They enter no
encoder, no selector, no donor rule and no normalisation statistic.
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

from scripts.qpsmp_data import stable_seed                          # noqa: E402
from scripts.stageR0_retrieval_falsification import (               # noqa: E402
    component_bootstrap, component_target_mean,
)
from tools.research.a2_readiness_v2 import _frozen                  # noqa: E402
from tools.research.a2_exact_probe.operator import (                # noqa: E402
    A2MomentOperator, ScalarLevelOperator, SharedMomentOperator,
    tanimoto_transport,
)

SUPPORT_SIZES = (0, 1, 2, 3, 5)
RANKS = (4, 8, 16)
STEPS = 600
LEARNING_RATE = 3e-3
FOLDS = 5
TRAIN_EPISODES_PER_TARGET = 8
EVAL_DRAWS_PER_TARGET = 4
QUERY_SIZE = 16


# ---------------------------------------------------------------------------
# episodes
# ---------------------------------------------------------------------------

class Corpus:
    """Per-target ligand blocks, sliced into episodes without touching labels."""

    def __init__(self, path: Path, device: str):
        with np.load(path) as stored:
            self.data = {key: stored[key] for key in stored.files}
        self.device = device
        owner = self.data["target_index"]
        self.blocks = [np.flatnonzero(owner == index)
                       for index in range(int(owner.max()) + 1)]
        self.components = self.data["component_of_target"]
        self.target_names = [str(name) for name in self.data["target_names"]]
        self.ligand_ids = [str(name) for name in self.data["ligand_ids"]]

    def episode_identity(self, episodes) -> list:
        """Episodes named by biology, for cross-process reproducibility checks."""
        return [(self.target_names[target],
                 tuple(self.ligand_ids[int(i)] for i in support),
                 tuple(self.ligand_ids[int(i)] for i in query))
                for target, support, query in episodes]

    def tensor(self, key: str) -> torch.Tensor:
        return torch.as_tensor(self.data[key], dtype=torch.float32,
                               device=self.device)

    def nested_banks(self, draws: int, seed: int,
                     support_sizes=SUPPORT_SIZES,
                     query_size: int = QUERY_SIZE) -> dict:
        """Nested banks: one query panel per (target, draw), shared across k.

        Reproduces `QPSMPData.fixed_nested_episode_banks` exactly. Two defects
        in the superseded version are repaired here:

        1. **Seeding.** It used Python's built-in `hash()` over a tuple
           containing a string, which is salted per process, so the "fixed"
           banks differed between runs. `stable_seed` is sha256-based, and
           keying on the target *name* rather than its positional index also
           makes the bank immune to renumbering by a changed ligand cap.
        2. **Nesting.** It sliced the query panel as
           `order[k : k + query_size]`, so the panel *moved* with the support
           size and the k-curve compared different queries. The nested contract
           takes the panel once, after the largest support, and gives each k a
           prefix of the same support ordering. Without this, a change across k
           mixes the support effect with a change of population.

        Within a single k the superseded run was still internally paired — every
        arm saw the same panel — so its per-k arm rankings stand. Its k-curve
        does not.
        """
        sizes = tuple(sorted(set(int(k) for k in support_sizes)))
        max_support = sizes[-1]
        banks = {size: [] for size in sizes}
        for target, block in enumerate(self.blocks):
            if len(block) < max_support + 2:      # largest support + 2 queries
                continue
            for draw in range(draws):
                rng = np.random.default_rng(stable_seed(
                    "a2-exact", seed, self.target_names[target], draw)
                    % (2 ** 32))
                order = rng.permutation(block)
                query = order[max_support:max_support + query_size]
                if len(query) < 2:
                    continue
                for size in sizes:
                    banks[size].append((target, order[:size], query))
        return banks


def episode_tensors(corpus: Corpus, feature_key: str, target: int,
                    support: np.ndarray, query: np.ndarray, f0_key: str,
                    labels: np.ndarray | None = None):
    features = corpus.data[feature_key]
    f0 = corpus.data[f0_key]
    y = corpus.data["y"] if labels is None else labels
    device = corpus.device
    to = lambda a: torch.as_tensor(a, dtype=torch.float32, device=device)
    support_residual = to(y[support] - f0[support]) if len(support) else \
        torch.zeros(0, device=device)
    return (to(features[support]) if len(support)
            else torch.zeros(0, features.shape[1], device=device),
            support_residual,
            to(features[query]),
            to(y[query] - f0[query]),          # the target of δ
            to(f0[query]), to(y[query]))


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------

def correlation(prediction: np.ndarray, truth: np.ndarray) -> float:
    p, y = prediction - prediction.mean(), truth - truth.mean()
    denominator = float(np.sqrt((p ** 2).mean()) * np.sqrt((y ** 2).mean()))
    return float((p * y).mean() / denominator) if denominator > 1e-12 else 0.0


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


# ---------------------------------------------------------------------------
# training and scoring
# ---------------------------------------------------------------------------

def train_operator(operator, corpus: Corpus, feature_key: str, f0_key: str,
                   episodes_by_k: dict, steps: int, learning_rate: float,
                   seed: int, labels: np.ndarray | None = None):
    trainable = [p for p in operator.parameters() if p.requires_grad]
    if not trainable:
        return operator
    optimiser = torch.optim.Adam(trainable, lr=learning_rate)
    rng = np.random.default_rng(seed)
    sizes = [k for k in SUPPORT_SIZES if k > 0 and episodes_by_k.get(k)]
    operator.train()
    for step in range(steps):
        k = sizes[step % len(sizes)]
        batch = episodes_by_k[k]
        picks = rng.choice(len(batch), size=min(16, len(batch)), replace=False)
        optimiser.zero_grad(set_to_none=True)
        loss = 0.0
        for pick in picks:
            target, support, query = batch[int(pick)]
            sf, sr, qf, target_delta, _, _ = episode_tensors(
                corpus, feature_key, target, support, query, f0_key, labels)
            delta = operator(sf, sr, qf)
            loss = loss + (delta - target_delta).square().mean()
        (loss / len(picks)).backward()
        optimiser.step()
    operator.eval()
    return operator


def score(operator, corpus: Corpus, feature_key: str, f0_key: str,
          episodes, fingerprint_arm: bool = False) -> dict:
    rows, spreads = [], []
    with torch.no_grad():
        for target, support, query in episodes:
            sf, sr, qf, _, f0_query, y_query = episode_tensors(
                corpus, feature_key, target, support, query, f0_key)
            if fingerprint_arm:
                fp = corpus.tensor("fingerprint")
                delta = tanimoto_transport(
                    fp[support] if len(support) else fp[:0], fp[query], sr)
            else:
                delta = operator(sf, sr, qf)
            prediction = (f0_query + delta).cpu().numpy()
            truth = y_query.cpu().numpy()
            component = int(corpus.components[target])
            rows.append({
                "component": component, "target": target,
                "mse_pk": float(((prediction - truth) ** 2).mean()),
                "r": correlation(prediction, truth),
                "ci": concordance(prediction, truth),
                "spearman": spearman(prediction, truth),
                "query_spread_pk": float(delta.cpu().numpy().std()),
            })
            spreads.append(float(delta.cpu().numpy().std()))
    if not rows:
        return {}
    out = {field: component_target_mean(
        (r["component"], r["target"], r[field]) for r in rows)
        for field in ("mse_pk", "r", "ci", "spearman", "query_spread_pk")}
    out["episodes"] = len(rows)
    out["mse_ci"] = component_bootstrap(
        [(r["component"], r["target"], r["mse_pk"]) for r in rows],
        _frozen.BOOTSTRAP_DRAWS, _frozen.BOOTSTRAP_SEED)
    out["r_ci"] = component_bootstrap(
        [(r["component"], r["target"], r["r"]) for r in rows],
        _frozen.BOOTSTRAP_DRAWS, _frozen.BOOTSTRAP_SEED)
    out["_rows"] = rows
    return out


def build_operator(kind: str, width: int, rank: int):
    if kind == "scalar_level":
        return ScalarLevelOperator()
    if kind == "shared_moment":
        return SharedMomentOperator(width, rank)
    if kind == "random_projection":
        return A2MomentOperator(width, rank, learn_projection=False)
    return A2MomentOperator(width, rank)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    train = Corpus(arguments.features / "meta_train.npz", arguments.device)
    val = Corpus(arguments.features / "meta_val.npz", arguments.device)

    train_episodes = train.nested_banks(
        TRAIN_EPISODES_PER_TARGET, _frozen.MODEL_FOLD_SEED)
    val_episodes = val.nested_banks(
        EVAL_DRAWS_PER_TARGET, _frozen.EVALUATION_SEED)

    # --- rank selection on meta_train component folds --------------------
    components = np.unique(train.components)
    assignment = {int(c): int(i % FOLDS) for i, c in enumerate(
        np.random.default_rng(_frozen.MODEL_FOLD_SEED).permutation(components))}

    def fold_score(rank: int, representation: str) -> float:
        width = train.data[f"correct__{representation}"].shape[1]
        held = []
        for fold in range(FOLDS):
            fit = {k: [e for e in v
                       if assignment[int(train.components[e[0]])] != fold]
                   for k, v in train_episodes.items()}
            check = {k: [e for e in v
                         if assignment[int(train.components[e[0]])] == fold]
                     for k, v in train_episodes.items()}
            if not check.get(1) or not fit.get(1):
                continue
            operator = A2MomentOperator(width, rank).to(arguments.device)
            train_operator(operator, train, f"correct__{representation}",
                           "correct__f0", fit, STEPS // 2, LEARNING_RATE,
                           _frozen.CONTROL_SEED + fold)
            block = [score(operator, train, f"correct__{representation}",
                           "correct__f0", check[k])
                     for k in (1, 2, 3, 5) if check.get(k)]
            held.append(float(np.mean([b["mse_pk"] for b in block if b])))
        return float(np.mean(held)) if held else float("nan")

    selection = {}
    for representation in ("embed", "max_state"):
        curve = {rank: fold_score(rank, representation) for rank in RANKS}
        best = min(RANKS, key=lambda r: (curve[r] if np.isfinite(curve[r])
                                         else np.inf))
        selection[representation] = {"curve": curve, "selected_rank": best}
        print(f"  rank curve {representation}: "
              f"{ {k: round(v, 4) for k, v in curve.items()} } -> {best}")

    # --- arms, all frozen before meta_val is read ------------------------
    noise = np.random.default_rng(_frozen.CONTROL_SEED + 7)
    for corpus in (train, val):
        corpus.data["random__feature"] = noise.normal(
            size=corpus.data["correct__embed"].shape).astype(np.float32)

    # Label-shuffled: permute pK *within a target*, so the support residuals
    # no longer correspond to their ligands while the target's mean and spread
    # are preserved exactly.
    shuffle = np.random.default_rng(_frozen.CONTROL_SEED + 11)
    shuffled = train.data["y"].copy()
    for block in train.blocks:
        shuffled[block] = shuffled[shuffle.permutation(block)]

    arms = {
        "scalar_level":       ("scalar_level", "embed", "correct__embed", "correct__f0", None),
        "tanimoto":           ("tanimoto", "embed", "correct__embed", "correct__f0", None),
        "a2_embed":           ("a2", "embed", "correct__embed", "correct__f0", None),
        "a2_max_state":       ("a2", "max_state", "correct__max_state", "correct__f0", None),
        "a2_ligand_only":     ("a2", "embed", "correct__ligand", "correct__f0", None),
        "a2_wrong_protein":   ("a2", "embed", "wrong__embed", "wrong__f0", None),
        "a2_reference_protein": ("a2", "embed", "reference__embed", "reference__f0", None),
        "a2_label_shuffled":  ("a2", "embed", "correct__embed", "correct__f0", shuffled),
        "a2_random_feature":  ("a2", "embed", "random__feature", "correct__f0", None),
        "a2_random_projection": ("random_projection", "embed", "correct__embed", "correct__f0", None),
        "shared_moment":      ("shared_moment", "embed", "correct__embed", "correct__f0", None),
    }

    results: dict = {}
    per_episode: dict = {}
    for name, (kind, representation, feature_key, f0_key, labels) in arms.items():
        rank = selection[representation]["selected_rank"]
        width = train.data[feature_key if feature_key in train.data
                           else "correct__embed"].shape[1]
        operator = build_operator(
            "a2" if kind in {"a2", "tanimoto"} else kind, width, rank
        ).to(arguments.device)
        # Every arm trains on the CORRECT protein. Corruptions are applied at
        # evaluation, so an arm cannot re-fit around its own corruption and
        # hide the dependence. The exceptions are the two whose corruption is
        # of the training signal itself (label shuffle) or of the feature space
        # (random features), which must be trained on what they are given.
        train_key = feature_key
        train_f0 = f0_key
        if name in {"a2_wrong_protein", "a2_reference_protein"}:
            train_key, train_f0 = "correct__embed", "correct__f0"
        if kind != "tanimoto":
            train_operator(operator, train, train_key, train_f0,
                           train_episodes, STEPS, LEARNING_RATE,
                           _frozen.CONTROL_SEED, labels)
        cell = {}
        for k in SUPPORT_SIZES:
            block = score(operator, val, feature_key, f0_key,
                          val_episodes[k], fingerprint_arm=(kind == "tanimoto"))
            per_episode[(name, k)] = block.pop("_rows", [])
            cell[str(k)] = block
        cell["rank"] = rank
        cell["trainable_parameters"] = int(sum(
            p.numel() for p in operator.parameters() if p.requires_grad))
        results[name] = cell
        zero, one, five = cell["0"], cell["1"], cell["5"]
        print(f"  {name:<22} k0 MSE {zero['mse_pk']:.4f} | k1 MSE "
              f"{one['mse_pk']:.4f} spread {one['query_spread_pk']:.4f} | "
              f"k5 MSE {five['mse_pk']:.4f} r {five['r']:+.3f}")

    # --- the preregistered gates, as paired component bootstraps ----------
    # Episodes are identical across arms (same bank, same draws), so the
    # contrasts are genuinely paired per (target, draw).
    def paired(left: str, right: str, k: int, field: str = "mse_pk") -> dict:
        a, b = per_episode.get((left, k), []), per_episode.get((right, k), [])
        if not a or len(a) != len(b):
            return {}
        # Positive = `left` is better. MSE is a loss, so the sign flips.
        sign = -1.0 if field == "mse_pk" else 1.0
        return component_bootstrap(
            [(x["component"], f"{x['target']}:{i}",
              sign * (x[field] - y[field]))
             for i, (x, y) in enumerate(zip(a, b))],
            _frozen.BOOTSTRAP_DRAWS, _frozen.BOOTSTRAP_SEED)

    gates = {
        "beats_scalar_level": {
            str(k): paired("a2_embed", "scalar_level", k)
            for k in SUPPORT_SIZES if k > 0},
        "beats_tanimoto": {
            str(k): paired("a2_embed", "tanimoto", k)
            for k in SUPPORT_SIZES if k > 0},
        "degrades_under_wrong_protein": {
            str(k): paired("a2_embed", "a2_wrong_protein", k)
            for k in SUPPORT_SIZES if k > 0},
        "degrades_under_label_shuffle": {
            str(k): paired("a2_embed", "a2_label_shuffled", k)
            for k in SUPPORT_SIZES if k > 0},
        "beats_shared_moment": {
            str(k): paired("a2_embed", "shared_moment", k)
            for k in SUPPORT_SIZES if k > 0},
        "beats_random_projection": {
            str(k): paired("a2_embed", "a2_random_projection", k)
            for k in SUPPORT_SIZES if k > 0},
    }

    payload = {
        "schema": "MetaSieve.A2ExactProbe.Episodic.v1",
        "gates": gates,
        "split": "meta_val",
        "frozen_design": _frozen.frozen_manifest(),
        "probe": {"ranks": list(RANKS), "steps": STEPS,
                  "learning_rate": LEARNING_RATE, "folds": FOLDS,
                  "train_episodes_per_target": TRAIN_EPISODES_PER_TARGET,
                  "eval_draws_per_target": EVAL_DRAWS_PER_TARGET,
                  "optimiser": "Adam; no closed form, no pseudoinverse, "
                               "no inner loop, no test-time gradient"},
        "rank_selection_on_meta_train_folds": selection,
        "arms": results,
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    report(payload)
    print(f"\nwrote {arguments.output}")
    return 0


def report(payload: dict) -> None:
    print(f"\n{'arm':<24}{'k':>3}{'MSE pK':>10}{'r':>8}{'CI':>7}"
          f"{'spearman':>10}{'spread pK':>11}")
    for name, cell in payload["arms"].items():
        for k in ("0", "1", "2", "3", "5"):
            block = cell.get(k) or {}
            if not block:
                continue
            print(f"{name if k == '0' else '':<24}{k:>3}"
                  f"{block['mse_pk']:>10.4f}{block['r']:>+8.3f}"
                  f"{block['ci']:>7.3f}{block['spearman']:>+10.3f}"
                  f"{block['query_spread_pk']:>11.4f}")

    print("\npreregistered gates — paired component bootstrap on k=0 MSE pK, "
          "positive = a2_embed better:")
    for gate, block in payload["gates"].items():
        cells = []
        for k in ("1", "2", "3", "5"):
            interval = block.get(k) or {}
            if not interval:
                continue
            mark = ("PASS" if interval["lo"] > 0
                    else "fail" if interval["hi"] < 0 else "unres")
            cells.append(f"k{k} {interval['mean']:+.4f} {mark}")
        print(f"  {gate:<30}{'  '.join(cells)}")


if __name__ == "__main__":
    raise SystemExit(main())

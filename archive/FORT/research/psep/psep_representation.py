"""R2: is the surviving gain the objective or the representation?

R1 produced one result that overturns how the previous gate was read.  A plain
multitask MLP on 1024 Morgan bits + 10 descriptors, trained with **MSE**, scores
**-0.0025** within-document concordance against the frozen linear base.  Yet the
`cnp` arm of the operator gate -- which was trained with a **bounded
smoothed-concordance surrogate** and which ignored its support entirely -- scored
**+0.0274** on the same substrate.

So the "+0.0274 target-agnostic representation" was probably not a representation
effect at all.  It was an *objective* effect.  This gate separates them, and fixes
R1's second defect at the same time: R1 tested on 20 % of components (76), which
left its largest arm (`centred`, +0.0195 vs ERM) unresolved.

**Factorial.**  Two objectives crossed with three context treatments, everything
else identical -- same encoder, same capacity, same optimiser, same budget:

    objective   mse      standard multitask DTA regression
                rank     bounded smoothed-concordance surrogate on
                         within-document training pairs
    context     raw      no context handling
                centred  target is the residual from its document mean
                fixed    y ~ f(x) + b_c, per-document intercept dropped at test

**Power.**  5-fold cross-fitting over homology components: every component serves
as test exactly once, so the estimate uses all ~379 components instead of 76, and
no component is ever trained on when it is scored.

**Speed.**  R1 took 75 minutes because it densified sparse batches on the CPU.
The full design is 263 318 x 1034 float32 = 1.09 GB, which fits on the GPU, so it
is materialised once and indexed there.

Primary metric is unchanged and therefore comparable to every earlier gate:
within-document pair concordance on the scaffold/document/assay-separated rows,
gain over the frozen cross-fitted linear base, bootstrapped over components.

Reads the `discover` role only.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F

from research.psep.psep_d0 import DEFAULT_SUBSTRATE, SEED, build_splits, paired_bootstrap
from research.psep.psep_operator import concordance
from research.psep.psep_transfer import MIN_HALF_PAIRS, within_document_pairs

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "psep_representation_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "psep_representation_records_2026-08-02.parquet"

FOLDS = 5
HIDDEN = (512, 256)
EMBED = 256
DROPOUT = 0.1
EPOCHS = 30
PATIENCE = 5
LEARNING_RATE = 1e-3
WEIGHT_DECAY = 1e-5
MSE_BATCH = 4096
RANK_DOCUMENTS = 96          # documents per ranking batch
RANK_STEPS = 150             # ranking steps per epoch
TAU = 0.5
MIN_DOC_ROWS = 2
ARMS = (("mse", "raw"), ("mse", "centred"), ("mse", "fixed"),
        ("rank", "raw"), ("rank", "centred"))
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class Encoder(nn.Module):
    def __init__(self, dimension: int):
        super().__init__()
        layers: list[nn.Module] = []
        previous = dimension
        for width in HIDDEN:
            layers += [nn.Linear(previous, width), nn.GELU(), nn.Dropout(DROPOUT)]
            previous = width
        layers += [nn.Linear(previous, EMBED), nn.GELU()]
        self.body = nn.Sequential(*layers)
        self.head = nn.Linear(EMBED, 1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.head(self.body(x)).squeeze(-1)


def component_fold(component) -> int:
    return int(sha256(f"{SEED}:r2fold:{component}".encode()).hexdigest()[:8], 16) % FOLDS


def smoothed_concordance(prediction, label, left, right):
    gap = prediction[left] - prediction[right]
    sign = torch.sign(label[left] - label[right])
    return torch.sigmoid(-gap * sign / TAU).mean()


def train(
    objective: str, context_mode: str, design: torch.Tensor, target: torch.Tensor,
    affinity: torch.Tensor, context: torch.Tensor, train_rows: np.ndarray,
    val_rows: np.ndarray, document_batches: list[np.ndarray], seed: int,
) -> tuple[Encoder, dict[str, object]]:
    torch.manual_seed(seed)
    model = Encoder(design.shape[1]).to(DEVICE)
    parameters = list(model.parameters())
    effects = None
    if context_mode == "fixed":
        effects = nn.Embedding(int(context.max().item()) + 1, 1).to(DEVICE)
        nn.init.zeros_(effects.weight)
        parameters += list(effects.parameters())
    optimiser = torch.optim.Adam(parameters, lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY)

    rng = np.random.default_rng(seed)
    train_tensor = torch.as_tensor(train_rows, device=DEVICE)
    history, best, best_state, stale = [], -np.inf, None, 0
    for epoch in range(EPOCHS):
        model.train()
        total, steps = 0.0, 0
        if objective == "mse":
            order = train_tensor[torch.randperm(len(train_tensor), device=DEVICE)]
            for start in range(0, len(order), MSE_BATCH):
                batch = order[start : start + MSE_BATCH]
                optimiser.zero_grad()
                prediction = model(design[batch])
                if context_mode == "fixed":
                    prediction = prediction + effects(context[batch]).squeeze(-1)
                loss = F.mse_loss(prediction, target[batch])
                loss.backward()
                torch.nn.utils.clip_grad_norm_(parameters, 5.0)
                optimiser.step()
                total += float(loss.item())
                steps += 1
        else:
            for _ in range(RANK_STEPS):
                chosen = rng.choice(len(document_batches),
                                    size=min(RANK_DOCUMENTS, len(document_batches)), replace=False)
                rows, left, right, offset = [], [], [], 0
                for index in chosen:
                    group = document_batches[index]
                    count = len(group)
                    a, b = np.triu_indices(count, k=1)
                    rows.append(group)
                    left.append(a + offset)
                    right.append(b + offset)
                    offset += count
                if not rows:
                    continue
                batch = torch.as_tensor(np.concatenate(rows), device=DEVICE)
                left_t = torch.as_tensor(np.concatenate(left), device=DEVICE)
                right_t = torch.as_tensor(np.concatenate(right), device=DEVICE)
                optimiser.zero_grad()
                prediction = model(design[batch])
                loss = smoothed_concordance(prediction, target[batch], left_t, right_t)
                loss.backward()
                torch.nn.utils.clip_grad_norm_(parameters, 5.0)
                optimiser.step()
                total += float(loss.item())
                steps += 1

        score = validation_concordance(model, design, affinity, val_rows)
        history.append({"epoch": epoch, "loss": total / max(steps, 1), "val_ci": score})
        if score > best + 1e-5:
            best, stale = score, 0
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            stale += 1
            if stale >= PATIENCE:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    return model, {"best_val_ci": best, "epochs_run": len(history), "history": history}


@torch.no_grad()
def validation_concordance(model, design, affinity, val_rows) -> float:
    """Model selection on within-document concordance -- the deployed metric --
    so no arm is penalised for being tuned against a different objective."""

    model.eval()
    scores = []
    for group in val_rows:
        if len(group) < 2:
            continue
        batch = torch.as_tensor(group, device=DEVICE)
        prediction = model(design[batch]).cpu().numpy()
        label = affinity[batch].cpu().numpy()
        left, right = np.triu_indices(len(group), k=1)
        keep = label[left] != label[right]
        if keep.sum() < 4:
            continue
        left, right = left[keep], right[keep]
        agree = np.sign(prediction[left] - prediction[right]) == np.sign(label[left] - label[right])
        scores.append(float(agree.mean()))
    return float(np.mean(scores)) if scores else float("nan")


def run(substrate_dir: Path, role: str, output: Path, records_path: Path, seed: int) -> dict[str, object]:
    started = time.time()
    from scipy.sparse import load_npz

    rows = pd.read_parquet(substrate_dir / "rows.parquet")
    rows = rows.loc[rows.role == role].reset_index(drop=True)
    bits = load_npz(substrate_dir / "morgan.npz").tocsr()
    descriptors = np.load(substrate_dir / "descriptors.npy")
    index = rows.structure_row.to_numpy()
    scale = np.maximum(descriptors.std(axis=0), 1e-6)
    standardised = ((descriptors - descriptors.mean(axis=0)) / scale)[index]
    dense = np.hstack([np.asarray(bits[index].todense(), dtype=np.float32),
                       standardised.astype(np.float32)])
    print(f"design {dense.shape} = {dense.nbytes / 1e9:.2f} GB -> {DEVICE}", flush=True)
    design = torch.as_tensor(dense, device=DEVICE)
    del dense

    affinity_np = rows.affinity.to_numpy(dtype=np.float64)
    base_np = rows.base.to_numpy(dtype=np.float64)
    documents = rows.docs.astype(str).to_numpy()
    context_np = pd.factorize(documents)[0]
    document_mean = pd.DataFrame({"c": context_np, "y": affinity_np}).groupby("c")["y"].transform("mean").to_numpy()

    affinity = torch.as_tensor(affinity_np, dtype=torch.float32, device=DEVICE)
    centred = torch.as_tensor(affinity_np - document_mean, dtype=torch.float32, device=DEVICE)
    context = torch.as_tensor(context_np, dtype=torch.long, device=DEVICE)

    folds = rows.component.map(component_fold).to_numpy()
    splits = [s for s in build_splits(rows) if s.regime == "separated"]
    units = []
    for split in splits:
        left, right = within_document_pairs(documents[split.evaluation], affinity_np[split.evaluation])
        if len(left) < 2 * MIN_HALF_PAIRS:
            continue
        split.pair_left, split.pair_right = left, right
        units.append(split)
    print(f"{len(rows)} rows | {len(units)} units | "
          f"{rows.component.nunique()} components | folds {FOLDS}", flush=True)

    document_groups: dict[int, list[np.ndarray]] = {}
    frame = pd.DataFrame({"doc": context_np, "fold": folds, "row": np.arange(len(rows))})
    for (fold, _), group in frame.groupby(["fold", "doc"], sort=False):
        if len(group) >= MIN_DOC_ROWS:
            document_groups.setdefault(int(fold), []).append(group.row.to_numpy())

    records: list[dict[str, object]] = []
    training: dict[str, object] = {}
    for objective, context_mode in ARMS:
        name = f"{objective}_{context_mode}"
        print(f"=== {name} ===", flush=True)
        target = centred if context_mode == "centred" else affinity
        per_fold = []
        for fold in range(FOLDS):
            test_mask = folds == fold
            validation_fold = (fold + 1) % FOLDS
            val_mask = folds == validation_fold
            train_mask = ~test_mask & ~val_mask
            train_rows = np.flatnonzero(train_mask)
            train_documents = [g for f, groups in document_groups.items() if f not in (fold, validation_fold)
                               for g in groups]
            val_groups = [g for g in document_groups.get(validation_fold, [])][:2000]
            model, info = train(objective, context_mode, design, target, affinity, context,
                                train_rows, val_groups, train_documents, seed + fold)
            per_fold.append({"fold": fold, **{k: v for k, v in info.items() if k != "history"}})
            print(f"  fold {fold}: val_ci={info['best_val_ci']:.4f} epochs={info['epochs_run']}", flush=True)

            model.eval()
            with torch.no_grad():
                for unit in units:
                    if not test_mask[unit.evaluation[0]]:
                        continue
                    batch = torch.as_tensor(unit.evaluation, device=DEVICE)
                    prediction = model(design[batch]).cpu().numpy().astype(np.float64)
                    label = affinity_np[unit.evaluation]
                    reference = base_np[unit.evaluation]
                    records.append({
                        "arm": name, "objective": objective, "context": context_mode,
                        "fold": fold, "unit": unit.unit, "component": unit.component,
                        "endpoint": unit.endpoint,
                        "ci_within": concordance(prediction, label, unit.pair_left, unit.pair_right),
                        "base_ci_within": concordance(reference, label, unit.pair_left, unit.pair_right),
                    })
        training[name] = per_fold

    table = pd.DataFrame.from_records(records)
    table["gain"] = table.ci_within - table.base_ci_within
    summary: dict[str, object] = {}
    for arm, part in table.groupby("arm"):
        summary[str(arm)] = {
            "gain_vs_base": paired_bootstrap(part, "gain"),
            "by_endpoint": {e: paired_bootstrap(p, "gain") for e, p in part.groupby("endpoint")},
        }
    reference_arm = "mse_raw"
    for arm in summary:
        left = table.loc[table.arm == arm, ["component", "unit", "gain"]]
        right = table.loc[table.arm == reference_arm, ["component", "unit", "gain"]]
        merged = left.merge(right, on=["component", "unit"], suffixes=("", "_ref"))
        merged["delta"] = merged.gain - merged.gain_ref
        summary[arm]["versus_mse_raw"] = paired_bootstrap(merged, "delta")

    def mean_of(arm: str) -> float:
        return summary[arm]["gain_vs_base"]["mean"]

    objective_effect = mean_of("rank_raw") - mean_of("mse_raw")
    context_effect = mean_of("mse_centred") - mean_of("mse_raw")
    combined = mean_of("rank_centred") - mean_of("mse_raw")
    winners = [a for a in summary if a != reference_arm
               and summary[a]["versus_mse_raw"]["lower95"] > 0.005]
    verdict = ("REPRESENTATION_IMPROVEMENT_CONFIRMED" if winners
               else "NO_ARM_BEATS_STANDARD_MULTITASK_DTA")

    payload = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "device": str(DEVICE),
        "firewall": {"evaluated_role": role, "validate_read": False, "confirm_read": False,
                     "scheme": "5-fold component cross-fitting; every component tested once"},
        "protocol": {
            "seed": seed, "folds": FOLDS, "arms": [f"{o}_{c}" for o, c in ARMS],
            "features": "1024 Morgan bits + 10 standardised descriptors (dense, on device)",
            "model_selection": "within-document concordance on a held-out fold",
            "metric": "within-document pair concordance vs the frozen cross-fitted linear base",
            "capacity_matched": True,
            "reference_arm": reference_arm,
        },
        "counts": {"rows": int(len(rows)), "units": len(units),
                   "components": int(rows.component.nunique())},
        "training": training,
        "summary": summary,
        "decomposition": {
            "objective_effect_rank_minus_mse": objective_effect,
            "context_effect_centred_minus_raw": context_effect,
            "combined_rank_centred_minus_mse_raw": combined,
            "interaction": combined - objective_effect - context_effect,
        },
        "winners": winners,
        "verdict": verdict,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_parquet(records_path, index=False)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="R2 objective vs invariance, cross-fitted")
    parser.add_argument("--substrate", type=Path, default=DEFAULT_SUBSTRATE)
    parser.add_argument("--role", default="discover")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    parser.add_argument("--seed", type=int, default=SEED)
    arguments = parser.parse_args()
    payload = run(arguments.substrate, arguments.role, arguments.output, arguments.records, arguments.seed)
    print(json.dumps({"verdict": payload["verdict"], "winners": payload["winners"],
                      "decomposition": payload["decomposition"]}, indent=2))
    print(f"\n{'arm':<16s}{'gain vs base':>26s}{'vs mse_raw':>26s}")
    for arm, cell in payload["summary"].items():
        gain, versus = cell["gain_vs_base"], cell["versus_mse_raw"]
        print(f"{arm:<16s}{gain['mean']:+.4f} [{gain['lower95']:+.4f},{gain['upper95']:+.4f}]"
              f"{versus['mean']:+.4f} [{versus['lower95']:+.4f},{versus['upper95']:+.4f}]"
              f"  n={gain['components']}")


if __name__ == "__main__":
    main()

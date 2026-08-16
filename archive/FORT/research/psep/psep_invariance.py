"""R1: can context-invariant training improve the transferable representation?

Three adaptation hypotheses are closed (latent task state, support-conditioned
operators, source routing): every target-conditional signal collapsed under
provenance separation.  What survived is *unconditional* -- a target-agnostic
nonlinear chemistry head gained +0.0274 within-document concordance on held-out
components.  So the question becomes representation learning under the
decomposition the programme measured:

    y = f(x, p) + g(c) + eps,     Var(g(c)) = 68 % of residual variance,
                                  SD(g(c)) ~ 1.7 pKi, group is affine.

**The design constraint that makes this falsifiable.**  Within-document
concordance is *already exactly invariant* to an additive `g(c)`: any per-document
constant cancels in every within-document pair.  An invariance mechanism therefore
cannot help by cleaning up the metric.  Its only route to a gain is indirect and
substantive -- by not spending model capacity on `c`, it should learn a better
`f`.  That is the hypothesis, and it is refutable at matched capacity.

Arms.  Identical encoder, identical capacity, identical optimiser and budget;
only the treatment of context changes.

  erm            plain MSE on y                              (standard multitask DTA)
  fixed_effects  y ~ f(x) + b_c with a free per-document intercept, discarded at
                 test time.  This is the decomposition written down literally --
                 a mixed-effects/within estimator, not an adversary.
  centred        regress on y minus its document mean         (G0-invariant target)
  irm            ERM + IRMv1 gradient penalty, documents as environments
  dro            GroupDRO over document groups                (worst-case context)
  adversarial    gradient-reversal head predicting the document's mean residual
                 (removes context-predictive directions from the representation)

Evaluation, on held-out `meta_test` components, on exactly the rows every earlier
gate used (scaffold + document + assay separated):

  ci_within        primary; gain over the frozen cross-fitted linear base
  rmse, spearman   absolute and monotone accuracy
  context_r2       R^2 of a linear probe from the learned representation to the
                   row's document-mean residual.  This is requirement "remove
                   document dependence", measured rather than asserted.

`meta_test` components are never trained on.  `validate` and `confirm` roles are
never opened.
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
from research.psep.psep_operator import assign_meta_role, concordance, unit_metrics
from research.psep.psep_transfer import MIN_HALF_PAIRS, within_document_pairs

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "reports" / "active" / "psep_invariance_2026-08-02.json"
DEFAULT_RECORDS = ROOT / "reports" / "active" / "psep_invariance_records_2026-08-02.parquet"

HIDDEN = (512, 256)
EMBED = 256
DROPOUT = 0.1
BATCH = 1024
EPOCHS = 25
PATIENCE = 5
LEARNING_RATE = 1e-3
WEIGHT_DECAY = 1e-5
IRM_WEIGHT = 1.0
ADVERSARY_WEIGHT = 1.0
DRO_STEP = 0.01
MIN_CONTEXT_ROWS = 5
ARMS = ("erm", "fixed_effects", "centred", "irm", "dro", "adversarial")
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


class Encoder(nn.Module):
    """Shared trunk.  Every arm gets exactly this, so capacity is matched."""

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

    def represent(self, x: torch.Tensor) -> torch.Tensor:
        return self.body(x)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.head(self.represent(x)).squeeze(-1)


class GradientReversal(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, weight):
        ctx.weight = weight
        return x.view_as(x)

    @staticmethod
    def backward(ctx, grad):
        return -ctx.weight * grad, None


def irm_penalty(prediction: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    """IRMv1: squared gradient of the risk w.r.t. a dummy unit scale."""

    scale = torch.tensor(1.0, device=prediction.device, requires_grad=True)
    loss = F.mse_loss(prediction * scale, target)
    grad = torch.autograd.grad(loss, [scale], create_graph=True)[0]
    return (grad ** 2).sum()


def build_features(substrate_dir: Path, rows: pd.DataFrame):
    from scipy.sparse import load_npz, hstack, csr_matrix

    bits = load_npz(substrate_dir / "morgan.npz").tocsr()
    descriptors = np.load(substrate_dir / "descriptors.npy")
    index = rows.structure_row.to_numpy()
    scale = np.maximum(descriptors.std(axis=0), 1e-6)
    standardised = (descriptors - descriptors.mean(axis=0)) / scale
    design = hstack([bits, csr_matrix(standardised.astype(np.float32))]).tocsr()
    return design[index]


def train_arm(
    arm: str, design, target: np.ndarray, context: np.ndarray, context_offset: np.ndarray,
    train_mask: np.ndarray, val_mask: np.ndarray, seed: int,
) -> tuple[Encoder, dict[str, object]]:
    torch.manual_seed(seed)
    model = Encoder(design.shape[1]).to(DEVICE)
    parameters = list(model.parameters())

    n_contexts = int(context.max()) + 1
    effects = None
    adversary = None
    if arm == "fixed_effects":
        effects = nn.Embedding(n_contexts, 1).to(DEVICE)
        nn.init.zeros_(effects.weight)
        parameters += list(effects.parameters())
    if arm == "adversarial":
        adversary = nn.Sequential(nn.Linear(EMBED, 128), nn.GELU(), nn.Linear(128, 1)).to(DEVICE)
        parameters += list(adversary.parameters())

    optimiser = torch.optim.Adam(parameters, lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY)
    train_rows = np.flatnonzero(train_mask)
    group_weights = None
    if arm == "dro":
        group_weights = torch.ones(n_contexts, device=DEVICE) / n_contexts

    rng = np.random.default_rng(seed)
    history, best, best_state, stale = [], np.inf, None, 0
    for epoch in range(EPOCHS):
        model.train()
        order = rng.permutation(train_rows)
        total, steps = 0.0, 0
        for start in range(0, len(order), BATCH):
            batch = order[start : start + BATCH]
            x = torch.as_tensor(np.asarray(design[batch].todense(), dtype=np.float32), device=DEVICE)
            y = torch.as_tensor(target[batch], dtype=torch.float32, device=DEVICE)
            c = torch.as_tensor(context[batch], dtype=torch.long, device=DEVICE)
            optimiser.zero_grad()

            representation = model.represent(x)
            prediction = model.head(representation).squeeze(-1)

            if arm == "fixed_effects":
                loss = F.mse_loss(prediction + effects(c).squeeze(-1), y)
            elif arm == "irm":
                loss = F.mse_loss(prediction, y)
                penalty = torch.zeros((), device=DEVICE)
                unique = torch.unique(c)
                counted = 0
                for value in unique[:32]:
                    mask = c == value
                    if int(mask.sum()) >= 4:
                        penalty = penalty + irm_penalty(prediction[mask], y[mask])
                        counted += 1
                if counted:
                    loss = loss + IRM_WEIGHT * penalty / counted
            elif arm == "dro":
                errors = (prediction - y) ** 2
                loss = torch.zeros((), device=DEVICE)
                unique = torch.unique(c)
                with torch.no_grad():
                    for value in unique:
                        mask = c == value
                        group_weights[value] = group_weights[value] * torch.exp(
                            DRO_STEP * errors[mask].mean().detach()
                        )
                    group_weights /= group_weights.sum()
                for value in unique:
                    mask = c == value
                    loss = loss + group_weights[value] * errors[mask].mean()
                loss = loss * len(unique)
            elif arm == "adversarial":
                loss = F.mse_loss(prediction, y)
                reversed_representation = GradientReversal.apply(representation, ADVERSARY_WEIGHT)
                offset = torch.as_tensor(context_offset[batch], dtype=torch.float32, device=DEVICE)
                loss = loss + F.mse_loss(adversary(reversed_representation).squeeze(-1), offset)
            else:                                   # erm, centred
                loss = F.mse_loss(prediction, y)

            loss.backward()
            torch.nn.utils.clip_grad_norm_(parameters, 5.0)
            optimiser.step()
            total += float(loss.item())
            steps += 1

        score = evaluate_mse(model, design, target, val_mask)
        history.append({"epoch": epoch, "train_loss": total / max(steps, 1), "val_mse": score})
        if score < best - 1e-5:
            best, stale = score, 0
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            stale += 1
            if stale >= PATIENCE:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    return model, {"history": history, "best_val_mse": best, "epochs_run": len(history)}


@torch.no_grad()
def evaluate_mse(model: Encoder, design, target: np.ndarray, mask: np.ndarray) -> float:
    model.eval()
    rows = np.flatnonzero(mask)
    total, count = 0.0, 0
    for start in range(0, len(rows), 4096):
        batch = rows[start : start + 4096]
        x = torch.as_tensor(np.asarray(design[batch].todense(), dtype=np.float32), device=DEVICE)
        prediction = model(x).cpu().numpy()
        total += float(((prediction - target[batch]) ** 2).sum())
        count += len(batch)
    return total / max(count, 1)


@torch.no_grad()
def predict(model: Encoder, design, rows: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    model.eval()
    outputs, embeddings = [], []
    for start in range(0, len(rows), 4096):
        batch = rows[start : start + 4096]
        x = torch.as_tensor(np.asarray(design[batch].todense(), dtype=np.float32), device=DEVICE)
        representation = model.represent(x)
        outputs.append(model.head(representation).squeeze(-1).cpu().numpy())
        embeddings.append(representation.cpu().numpy())
    return np.concatenate(outputs), np.concatenate(embeddings)


def context_dependence(embedding: np.ndarray, offset: np.ndarray) -> float:
    """R^2 of a ridge probe from representation to the row's context offset."""

    if len(embedding) < 32 or np.std(offset) < 1e-9:
        return float("nan")
    centred = embedding - embedding.mean(axis=0)
    y = offset - offset.mean()
    gram = centred.T @ centred + 1.0 * np.eye(centred.shape[1])
    weight = np.linalg.solve(gram, centred.T @ y)
    residual = y - centred @ weight
    return float(1.0 - residual.var() / max(y.var(), 1e-12))


def run(substrate_dir: Path, role: str, output: Path, records_path: Path, seed: int) -> dict[str, object]:
    started = time.time()
    rows = pd.read_parquet(substrate_dir / "rows.parquet")
    rows = rows.loc[rows.role == role].reset_index(drop=True)
    design = build_features(substrate_dir, rows)

    affinity = rows.affinity.to_numpy(dtype=np.float64)
    base = rows.base.to_numpy(dtype=np.float64)
    documents = rows.docs.astype(str).to_numpy()
    context = pd.factorize(documents)[0]
    frame = pd.DataFrame({"c": context, "r": affinity - base})
    grouped = frame.groupby("c")["r"]
    context_offset = grouped.transform("mean").to_numpy()
    counts = grouped.transform("size").to_numpy()
    context_offset = np.where(counts >= MIN_CONTEXT_ROWS, context_offset, 0.0)
    document_mean = pd.DataFrame({"c": context, "y": affinity}).groupby("c")["y"].transform("mean").to_numpy()

    meta_role = rows.component.map(assign_meta_role).to_numpy()
    train_mask = meta_role == "meta_train"
    val_mask = meta_role == "meta_val"
    print(f"{len(rows)} rows | train {train_mask.sum()} val {val_mask.sum()} "
          f"test {(meta_role=='meta_test').sum()} | contexts {context.max()+1} | {DEVICE}", flush=True)

    class Shim:
        pass

    shim = Shim()
    shim.rows = rows
    shim.affinity = affinity
    shim.base = base
    shim.residual = affinity - base
    splits = [s for s in build_splits(rows) if s.regime == "separated"]
    units = []
    for split in splits:
        component_role = assign_meta_role(split.component)
        if component_role != "meta_test":
            continue
        left, right = within_document_pairs(documents[split.evaluation], affinity[split.evaluation])
        if len(left) < 2 * MIN_HALF_PAIRS:
            continue
        split.pair_left, split.pair_right = left, right
        units.append(split)
    print(f"{len(units)} meta_test units / {len({u.component for u in units})} components", flush=True)

    results, training = [], {}
    for arm in ARMS:
        print(f"training {arm} ...", flush=True)
        target = (affinity - document_mean) if arm == "centred" else affinity
        model, info = train_arm(arm, design, target, context, context_offset,
                                train_mask, val_mask, seed)
        training[arm] = info
        print(f"  {arm}: val_mse={info['best_val_mse']:.4f} epochs={info['epochs_run']}", flush=True)

        for unit in units:
            prediction, embedding = predict(model, design, unit.evaluation)
            label = affinity[unit.evaluation]
            reference = base[unit.evaluation]
            # `centred` predicts a document-centred target, which is on a
            # different additive scale; within-document ranking is unaffected and
            # absolute metrics are restored by adding the base's document mean.
            absolute = prediction + (document_mean[unit.evaluation] if arm == "centred" else 0.0)
            unit.pair_left_arr = unit.pair_left
            record = {
                "arm": arm, "unit": unit.unit, "component": unit.component, "endpoint": unit.endpoint,
                "ci_within": concordance(prediction, label, unit.pair_left, unit.pair_right),
                "base_ci_within": concordance(reference, label, unit.pair_left, unit.pair_right),
                "rmse": float(np.sqrt(np.mean((absolute - label) ** 2))),
                "base_rmse": float(np.sqrt(np.mean((reference - label) ** 2))),
                "context_r2": context_dependence(embedding, context_offset[unit.evaluation]),
            }
            from scipy.stats import spearmanr
            if len(label) > 2 and np.std(prediction) > 1e-9:
                record["spearman"] = float(spearmanr(prediction, label)[0])
                record["base_spearman"] = float(spearmanr(reference, label)[0])
            record["gain_ci_within"] = record["ci_within"] - record["base_ci_within"]
            record["gain_rmse"] = record["rmse"] - record["base_rmse"]
            results.append(record)

    records = pd.DataFrame.from_records(results)
    summary: dict[str, object] = {}
    for arm, part in records.groupby("arm"):
        summary[str(arm)] = {
            "gain_ci_within": paired_bootstrap(part, "gain_ci_within"),
            "gain_rmse": paired_bootstrap(part, "gain_rmse"),
            "context_r2": paired_bootstrap(part, "context_r2"),
            "ci_within": float(part.ci_within.mean()),
        }
    erm = summary["erm"]["gain_ci_within"]["mean"]
    for arm in ARMS:
        part = records.loc[records.arm == arm, ["component", "unit", "gain_ci_within"]].copy()
        other = records.loc[records.arm == "erm", ["component", "unit", "gain_ci_within"]]
        merged = part.merge(other, on=["component", "unit"], suffixes=("", "_erm"))
        merged["versus_erm"] = merged.gain_ci_within - merged.gain_ci_within_erm
        summary[arm]["versus_erm"] = paired_bootstrap(merged, "versus_erm")

    winners = [a for a in ARMS if a != "erm" and summary[a]["versus_erm"]["lower95"] > 0.005]
    verdict = ("INVARIANCE_MECHANISM_IMPROVES_REPRESENTATION" if winners
               else "NO_INVARIANCE_MECHANISM_BEATS_PLAIN_MULTITASK_ERM")

    payload = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "runtime_seconds": round(time.time() - started, 1),
        "device": str(DEVICE),
        "firewall": {"evaluated_role": role, "validate_read": False, "confirm_read": False,
                     "reported_on": "meta_test components inside discover"},
        "protocol": {
            "seed": seed, "arms": list(ARMS), "hidden": list(HIDDEN), "embed": EMBED,
            "epochs": EPOCHS, "batch": BATCH, "features": "1024 Morgan bits + 10 standardised descriptors",
            "loss": "MSE (standard multitask DTA); arms differ only in context treatment",
            "metric": "within-document pair concordance vs the frozen cross-fitted linear base",
            "capacity_matched": True,
        },
        "counts": {"rows": int(len(rows)), "meta_test_units": len(units),
                   "meta_test_components": int(len({u.component for u in units})),
                   "contexts": int(context.max() + 1)},
        "training": training,
        "summary": summary,
        "winners": winners,
        "verdict": verdict,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    records.to_parquet(records_path, index=False)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser(description="R1 context-invariance gate")
    parser.add_argument("--substrate", type=Path, default=DEFAULT_SUBSTRATE)
    parser.add_argument("--role", default="discover")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--records", type=Path, default=DEFAULT_RECORDS)
    parser.add_argument("--seed", type=int, default=SEED)
    arguments = parser.parse_args()
    payload = run(arguments.substrate, arguments.role, arguments.output, arguments.records, arguments.seed)
    print(json.dumps({"verdict": payload["verdict"], "winners": payload["winners"]}, indent=2))
    print(f"\n{'arm':<15s}{'gain CI vs base':>24s}{'vs ERM':>24s}{'context R2':>18s}{'gain RMSE':>12s}")
    for arm in ARMS:
        cell = payload["summary"][arm]
        gain, versus, r2 = cell["gain_ci_within"], cell["versus_erm"], cell["context_r2"]
        print(f"{arm:<15s}{gain['mean']:+.4f} [{gain['lower95']:+.4f},{gain['upper95']:+.4f}]"
              f"{versus['mean']:+.4f} [{versus['lower95']:+.4f},{versus['upper95']:+.4f}]"
              f"{r2['mean']:>17.3f}{cell['gain_rmse']['mean']:>+12.3f}")


if __name__ == "__main__":
    main()

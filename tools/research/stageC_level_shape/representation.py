"""Stage 2 measurement: is there a ligand-varying interaction signal at all?

Stage B showed the within-target readout activations have pairwise cosine 0.997,
but cosine alone is weak evidence — two vectors can be nearly parallel and still
differ enough to order ligands. This measures the quantities that actually
decide it, on the frozen leak-free baseline:

* **protein-constant vs ligand-varying variance** — split each representation's
  total variance into the part that moves only between targets and the part that
  moves between the ligands of one target. The second is all a within-target
  ordering mechanism can ever use.
* **effective rank** of the within-target centered covariance — how many
  directions the ligand-varying part actually occupies.
* **Euclidean separation** between a target's ligands, in units of the
  representation's own scale.
* **a frozen linear probe on within-target affinity differences** — the direct
  test. Fitted on meta_train targets by ordinary SGD, weight decay chosen on
  component folds, read once on meta_val. If the ligand-varying subspace carried
  usable ordering information, this probe would find it.

No training of the model, no meta_test.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import QPSMPData, stable_seed                 # noqa: E402
from scripts.stageR0_retrieval_falsification import (                 # noqa: E402
    component_bootstrap,
)
from scripts.stageR6_compare_arms import load_arm                     # noqa: E402
from scripts.train_qpsmp import (                                     # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    compact_episode, normalized_episode, training_label_scale,
)
from tools.research.stageA_innerloop.train_meta import (              # noqa: E402
    align_atoms, encode_parts, episode_tensors,
)

SPLIT = ROOT / "dataset/processed/meta_fewshot/bindingdb_ki_double_cold_v1"
PANEL_SEED = 20260817
PANEL_MAX = 24
PANEL_MIN = 6
FOLDS = 5
FOLD_SEED = 20260818
PROBE_STEPS = 500
PROBE_LR = 3e-3
WEIGHT_DECAYS = (1e-4, 1e-3, 1e-2, 1e-1, 1.0)
BOOTSTRAP_DRAWS = 9999
BOOTSTRAP_SEED = 20260816
OUT = Path(__file__).resolve().parent / "REPRESENTATION.json"

REPRESENTATIONS = ("embed", "readout_hidden", "occupancy", "section")


def panels(data: QPSMPData, split: str):
    from scripts.qpsmp_data import EpisodeSpec
    component_of = {c["target_id"]: c["protein_group_40"] for c in data.cells}
    out = []
    for target in sorted(data.tasks[split]):
        rng = np.random.default_rng(stable_seed("stageC", PANEL_SEED, split, target))
        order = data._unique_ligand_order(data.tasks[split][target], rng)
        if len(order) < PANEL_MIN:
            continue
        cells = tuple(int(i) for i in order[:PANEL_MAX])
        out.append((component_of[target], target, EpisodeSpec(
            split, component_of[target], target, (), cells, target)))
    return out


def capture(model, parts) -> dict[str, np.ndarray]:
    task = encode_parts(model, parts)
    body = model.interaction_head
    hidden = body[1](F.linear(task.query_hidden, body[0].weight, body[0].bias))
    embed_width = task.query_embed.shape[-1]
    return {
        "embed": task.query_embed.squeeze(0).cpu().numpy(),
        "readout_hidden": hidden.squeeze(0).cpu().numpy(),
        "occupancy": task.query_occupancy.squeeze(0).cpu().numpy(),
        "section": task.query_hidden.squeeze(0).cpu().numpy()[:, embed_width:],
    }


def effective_rank(matrix: np.ndarray) -> float:
    """exp(entropy of the normalized eigenvalue spectrum)."""
    if matrix.shape[0] < 2:
        return float("nan")
    centered = matrix - matrix.mean(0, keepdims=True)
    values = np.linalg.svd(centered, compute_uv=False) ** 2
    total = values.sum()
    if total <= 1e-20:
        return 0.0
    p = values / total
    p = p[p > 1e-12]
    return float(np.exp(-(p * np.log(p)).sum()))


def train_probe(blocks, width, decay, seed, steps=PROBE_STEPS):
    torch.manual_seed(seed)
    probe = nn.Linear(width, 1, bias=False)
    nn.init.zeros_(probe.weight)
    optimizer = torch.optim.AdamW(probe.parameters(), lr=PROBE_LR,
                                  weight_decay=decay)
    tensors = [(torch.as_tensor(x, dtype=torch.float32),
                torch.as_tensor(y, dtype=torch.float32)) for x, y in blocks]
    for _ in range(steps):
        loss = torch.zeros(())
        for x, y in tensors:
            p = probe(x).squeeze(-1)
            loss = loss + F.mse_loss(p - p.mean(), y - y.mean())
        (loss / len(tensors)).backward()
        optimizer.step()
        optimizer.zero_grad(set_to_none=True)
    return probe


def within_r(probe, blocks) -> list[float]:
    out = []
    with torch.no_grad():
        for x, y in blocks:
            p = probe(torch.as_tensor(x, dtype=torch.float32)).squeeze(-1).numpy()
            pc, yc = p - p.mean(), y - y.mean()
            d = float(np.sqrt((pc ** 2).sum()) * np.sqrt((yc ** 2).sum()))
            out.append(float((pc * yc).sum() / d) if d > 1e-12 else 0.0)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    arguments = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    label_scale = training_label_scale(data)
    model, _, _ = load_arm(arguments.checkpoint, data, arguments.device)
    model.eval()

    store: dict[str, dict] = {}
    for split in ("meta_train", "meta_val"):
        rows = []
        for component, target, spec in panels(data, split):
            episode = compact_episode(normalized_episode(
                data.materialize(spec), label_scale))
            parts = align_atoms(episode_tensors(model, episode,
                                                arguments.device, torch.float32))
            with torch.no_grad():
                features = capture(model, parts)
            labels = (parts["query_y"].squeeze(0).cpu().numpy()
                      * label_scale.scale + label_scale.mean)
            rows.append({"component": component, "target": target,
                         "features": features, "labels": labels})
        store[split] = rows
        print(f"{split}: {len(rows)} panels")

    payload = {"schema": "MetaSieve.StageC.Representation.v1",
               "date": "2026-08-17",
               "checkpoint": str(arguments.checkpoint.resolve().relative_to(ROOT)),
               "panels": {s: len(store[s]) for s in store},
               "representations": {}, "meta_test": data.seal_record()}

    for name in REPRESENTATIONS:
        rows = store["meta_val"]
        widths = rows[0]["features"][name].shape[1]
        target_means = np.stack([r["features"][name].mean(0) for r in rows])
        between = float(target_means.var(0).sum())
        within = float(np.mean([
            r["features"][name].var(0).sum() for r in rows]))
        ranks = [effective_rank(r["features"][name]) for r in rows]
        separations = []
        for r in rows:
            f = r["features"][name]
            centre = f.mean(0, keepdims=True)
            spread = np.linalg.norm(f - centre, axis=1).mean()
            separations.append(float(spread / max(np.linalg.norm(centre), 1e-9)))
        block = {
            "width": int(widths),
            "protein_constant_variance": between,
            "ligand_varying_variance": within,
            "ligand_varying_share": within / max(between + within, 1e-20),
            "within_target_effective_rank_mean": float(np.mean(ranks)),
            "within_target_effective_rank_max_possible": float(
                min(PANEL_MAX, widths)),
            "relative_euclidean_separation": float(np.mean(separations)),
        }
        payload["representations"][name] = block
        print(f"{name:<16} width {widths:>4}  between {between:>10.4f}  "
              f"within {within:>10.4f}  ligand-share {block['ligand_varying_share']:.5f}  "
              f"eff.rank {np.mean(ranks):>6.2f}  sep {np.mean(separations):.5f}")

    # --- the direct test: a frozen probe on within-target differences -------
    print("\nfrozen linear probe on within-target centered affinity:")
    components = [r["component"] for r in store["meta_train"]]
    unique = sorted(set(components))
    order = np.random.default_rng(FOLD_SEED).permutation(len(unique))
    fold_of = {unique[int(i)]: rank % FOLDS for rank, i in enumerate(order)}
    for name in REPRESENTATIONS:
        train_blocks = [(r["features"][name], r["labels"]) for r in store["meta_train"]]
        folds = np.asarray([fold_of[r["component"]] for r in store["meta_train"]])
        width = train_blocks[0][0].shape[1]
        table = {}
        for decay in WEIGHT_DECAYS:
            scores = []
            for fold in range(FOLDS):
                fit = [b for b, f in zip(train_blocks, folds) if f != fold]
                held = [b for b, f in zip(train_blocks, folds) if f == fold]
                if not fit or not held:
                    continue
                probe = train_probe(fit, width, decay, seed=fold)
                scores.extend(within_r(probe, held))
            table[f"{decay:g}"] = float(np.mean(scores)) if scores else float("nan")
        best = max(table, key=lambda key: table[key])
        probe = train_probe(train_blocks, width, float(best), seed=0)
        val_blocks = [(r["features"][name], r["labels"]) for r in store["meta_val"]]
        scores = within_r(probe, val_blocks)
        pairs = [(r["component"], r["target"], s)
                 for r, s in zip(store["meta_val"], scores)]
        interval = component_bootstrap(pairs, BOOTSTRAP_DRAWS, BOOTSTRAP_SEED)
        payload["representations"][name]["within_target_probe"] = {
            "fold_table": table, "selected_weight_decay": float(best),
            "meta_val_within_target_r": interval,
        }
        print(f"  {name:<16} train-fold r {table[best]:+.4f} (decay {best})  "
              f"meta_val r {interval['mean']:+.4f} "
              f"[{interval['lo']:+.4f}, {interval['hi']:+.4f}]")

    OUT.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

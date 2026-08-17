"""D0: does the occupancy ordering signal survive stratification?

Stage C found a resolved within-target ordering signal in the 24-dim occupancy
representation (meta_val r +0.218 [+0.075, +0.367], component bootstrap). The
governing task asks whether it holds:

* across different components (per-component r on meta_val);
* under scaffold novelty (panels whose ligands share scaffolds with the fit
  folds vs novel scaffolds);
* under low ligand recall (held-out-fold panels whose query ligands appear in
  the fit folds vs never-seen ligands).

The double-cold split makes every meta_val ligand and scaffold novel by
construction, so scaffold/recall strata are measured on meta_train held-out
component folds. Probe: frozen linear map on within-target centered affinity,
weight decay 1.0 (Stage C's selected value), fitted per fold on the other
components. No model training. meta_test never constructed.
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

from scripts.qpsmp_data import QPSMPData, stable_seed                 # noqa: E402
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
DECAY = 1.0
OUT = Path(__file__).resolve().parent / "D0_OCCUPANCY_STRATA.json"


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


def capture(model, parts) -> np.ndarray:
    task = encode_parts(model, parts)
    return task.query_occupancy.squeeze(0).cpu().numpy()


def train_probe(blocks, width, decay, seed):
    torch.manual_seed(seed)
    probe = nn.Linear(width, 1, bias=False)
    nn.init.zeros_(probe.weight)
    optimizer = torch.optim.AdamW(probe.parameters(), lr=PROBE_LR,
                                  weight_decay=decay)
    tensors = [(torch.as_tensor(x, dtype=torch.float32),
                torch.as_tensor(y, dtype=torch.float32)) for x, y in blocks]
    for _ in range(PROBE_STEPS):
        loss = torch.zeros(())
        for x, y in tensors:
            p = probe(x).squeeze(-1)
            loss = loss + nn.functional.mse_loss(p - p.mean(), y - y.mean())
        (loss / len(tensors)).backward()
        optimizer.step()
        optimizer.zero_grad(set_to_none=True)
    return probe


def within_r(probe, block) -> float:
    x, y = block
    with torch.no_grad():
        p = probe(torch.as_tensor(x, dtype=torch.float32)).squeeze(-1).numpy()
    pc, yc = p - p.mean(), y - y.mean()
    d = float(np.sqrt((pc ** 2).sum()) * np.sqrt((yc ** 2).sum()))
    return float((pc * yc).sum() / d) if d > 1e-12 else 0.0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--checkpoint", type=Path,
                        default=ROOT / "report/meta_fewshot/stageB_complementary_20260817/T/checkpoint.pt")
    parser.add_argument("--device",
                        default="cuda" if torch.cuda.is_available() else "cpu")
    args = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=SPLIT)
    label_scale = training_label_scale(data)
    model, _, _ = load_arm(args.checkpoint, data, args.device)
    model.eval()

    store = {}
    for split in ("meta_train", "meta_val"):
        rows = []
        for component, target, spec in panels(data, split):
            episode = compact_episode(normalized_episode(
                data.materialize(spec), label_scale))
            parts = align_atoms(episode_tensors(model, episode,
                                                args.device, torch.float32))
            with torch.no_grad():
                features = capture(model, parts)
            labels = (parts["query_y"].squeeze(0).cpu().numpy()
                      * label_scale.scale + label_scale.mean)
            rows.append({"component": component, "target": target,
                         "features": features, "labels": labels,
                         "split": split, "spec": spec})
        store[split] = rows

    train_rows = store["meta_train"]
    val_rows = store["meta_val"]
    width = train_rows[0]["features"].shape[1]
    components = sorted({r["component"] for r in train_rows})
    order = np.random.default_rng(FOLD_SEED).permutation(len(components))
    fold_of = {components[int(i)]: rank % FOLDS for rank, i in enumerate(order)}

    scaffold_freq: dict[str, int] = {}
    for row in data._read_jsonl(Path(CORPUS) / "ligands.jsonl"):
        scaffold_freq[str(row["scaffold"])] = scaffold_freq.get(str(row["scaffold"]), 0) + 1

    def panel_scaffold_mode(row) -> str:
        from collections import Counter
        counts = Counter()
        for idx in row["spec_query"]:
            ligand = data.cells[int(idx)]["ligand_id"]
            for entry in data._read_jsonl(Path(CORPUS) / "ligands.jsonl"):
                if entry["drug_key"] == ligand:
                    counts[str(entry["scaffold"])] += 1
        return counts.most_common(1)[0][0] if counts else ""

    # index ligand->scaffold once
    ligand_scaffold = {}
    for entry in data._read_jsonl(Path(CORPUS) / "ligands.jsonl"):
        ligand_scaffold[str(entry["drug_key"])] = str(entry.get("scaffold", ""))

    train_ligand_set = set()
    for split in ("meta_train",):
        for target, indices in data.tasks[split].items():
            for idx in indices:
                train_ligand_set.add(data.cells[int(idx)]["ligand_id"])

    # meta_train held-out-fold strata
    strat = {"scaffold_shared": [], "scaffold_novel": [],
             "recall_high": [], "recall_low": []}
    per_fold_scores = []
    for fold in range(FOLDS):
        fit_rows = [r for r in train_rows if fold_of[r["component"]] != fold]
        held_rows = [r for r in train_rows if fold_of[r["component"]] == fold]
        if not fit_rows or not held_rows:
            continue
        probe = train_probe(
            [(r["features"], r["labels"]) for r in fit_rows], width, DECAY,
            seed=fold)
        fit_ligands = set()
        for r in fit_rows:
            target = r["target"]
            for idx in data.tasks["meta_train"][target]:
                fit_ligands.add(data.cells[int(idx)]["ligand_id"])
        median_freq = float(np.median([scaffold_freq.get(
            ligand_scaffold.get(data.cells[int(idx)]["ligand_id"], ""), 0)
            for r in fit_rows for idx in r["spec"].query]))
        for r in held_rows:
            score = within_r(probe, (r["features"], r["labels"]))
            per_fold_scores.append((r["component"], r["target"], score))
            freqs = [scaffold_freq.get(ligand_scaffold.get(
                data.cells[int(idx)]["ligand_id"], ""), 0)
                for idx in r["spec"].query]
            if float(np.mean(freqs)) >= median_freq:
                strat["scaffold_shared"].append(score)
            else:
                strat["scaffold_novel"].append(score)
            recall = float(np.mean([
                data.cells[int(idx)]["ligand_id"] in fit_ligands
                for idx in r["spec"].query]))
            if recall >= 0.5:
                strat["recall_high"].append(score)
            else:
                strat["recall_low"].append(score)

    # meta_val per-component scores with the all-meta_train probe
    probe_all = train_probe(
        [(r["features"], r["labels"]) for r in train_rows], width, DECAY, seed=0)
    per_component = {}
    val_scores = []
    for r in val_rows:
        score = within_r(probe_all, (r["features"], r["labels"]))
        val_scores.append(score)
        per_component.setdefault(r["component"], []).append(score)

    payload = {
        "schema": "MetaSieve.StageD.OccupancyStrata.v1",
        "date": "2026-08-17",
        "checkpoint": str(args.checkpoint.resolve().relative_to(ROOT)),
        "meta_train_fold_strata": {
            key: {"n": len(values), "mean_r": float(np.mean(values))
                  if values else float("nan")}
            for key, values in strat.items()},
        "meta_val_per_component": {
            comp: {"n": len(values), "mean_r": float(np.mean(values))}
            for comp, values in sorted(per_component.items())},
        "meta_val_overall_r": float(np.mean(val_scores)),
        "meta_test": data.seal_record(),
    }
    print(json.dumps(payload, indent=1))
    OUT.write_text(json.dumps(payload, indent=1, sort_keys=True), encoding="utf-8")
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""P0 leakage-safe AnchorDelta diagnostic.

The encoder is trained only on the fit components carved from TRAIN.  An
uncertified checkpoint is deliberately rejected because a checkpoint trained
before the holdout split can leak target-specific information into the gate.
The strict development roster remains untouched.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch.nn import functional as F

from model.anchordelta import AnchorDelta, EncodedAnchorDelta, anchorabsolute
from model.reorder import ReorderingModel
from research.shared.priorgate import limitrows, splitcomponents, wrongtargets
from scripts.metric import evaluateprotocol, pairedcomponents
from scripts.preprocess import normalizeligands
from scripts.train import fitbase, loadlabels, loadproteins, maketrainroster


REPORT = Path("reports/active/anchordelta_p0_smoke.json")


def samplepairs(
    frame: pd.DataFrame,
    pairspertarget: int,
    seed: int,
) -> list[tuple[str, int, int]]:
    """Sample a fixed number of unordered within-target pairs per target."""

    if pairspertarget <= 0:
        raise ValueError("pairspertarget must be positive")
    generator = np.random.default_rng(seed)
    pairs: list[tuple[str, int, int]] = []
    for target, group in frame.groupby("target", sort=True):
        indices = group.index.to_numpy(dtype=np.int64)
        if len(indices) < 2:
            continue
        count = min(pairspertarget, len(indices) * (len(indices) - 1) // 2)
        chosen: set[tuple[int, int]] = set()
        while len(chosen) < count:
            left, right = generator.choice(indices, size=2, replace=False)
            pair = tuple(sorted((int(left), int(right))))
            chosen.add(pair)
        pairs.extend((str(target), left, right) for left, right in sorted(chosen))
    if not pairs:
        raise ValueError("no target has at least two labeled rows")
    return pairs


def buildinteraction(checkpoint: Path | None = None) -> ReorderingModel:
    """Build an encoder, loading only an explicitly provenance-checked artifact."""

    payload = None if checkpoint is None else torch.load(checkpoint, map_location="cuda")
    backbone = (payload or {}).get("backbone") or "hybrid"
    model = ReorderingModel(
        proteindim=1280,
        liganddim=1034,
        backbone=str(backbone),
        proteinconditioned=True,
        interactiononly=False,
    ).cuda()
    if payload is not None:
        state = payload.get("state", payload)
        model.load_state_dict(state, strict=True)
    model.eval()
    return model


@torch.no_grad()
def encodeframe(
    model: EncodedAnchorDelta,
    frame: pd.DataFrame,
    proteins: dict[str, torch.Tensor],
    feature: np.ndarray,
) -> dict[str, tuple[torch.Tensor, torch.Tensor, np.ndarray]]:
    encoded: dict[str, tuple[torch.Tensor, torch.Tensor, np.ndarray]] = {}
    for target, group in frame.groupby("target", sort=True):
        indices = group.index.to_numpy(dtype=np.int64)
        sources = group.source_row.to_numpy(dtype=np.int64)
        ligands = torch.as_tensor(feature[sources], device="cuda", dtype=torch.float32)
        values = torch.as_tensor(group.affinity.to_numpy(dtype=np.float32), device="cuda")
        encoded[str(target)] = (model.encode(proteins[str(target)], ligands), values, indices)
    return encoded


def fithead(
    model: AnchorDelta,
    encoded: dict[str, tuple[torch.Tensor, torch.Tensor, np.ndarray]],
    pairs: list[tuple[str, int, int]],
    epochs: int,
    seed: int,
) -> list[dict[str, float]]:
    optimizer = torch.optim.AdamW(model.head.parameters(), lr=2e-3, weight_decay=1e-4)
    lookup = {
        target: {int(index): position for position, index in enumerate(indices)}
        for target, (_, _, indices) in encoded.items()
    }
    history: list[dict[str, float]] = []
    generator = np.random.default_rng(seed)
    for epoch in range(epochs):
        order = generator.permutation(len(pairs))
        losses: list[float] = []
        model.train()
        for position in order:
            target, left, right = pairs[int(position)]
            features, labels, _ = encoded[target]
            leftpos, rightpos = lookup[target][left], lookup[target][right]
            # Random orientation prevents an ordered-row shortcut.
            if generator.random() < 0.5:
                leftpos, rightpos = rightpos, leftpos
            prediction = model.delta(features[rightpos : rightpos + 1], features[leftpos : leftpos + 1])
            targetdelta = labels[rightpos : rightpos + 1] - labels[leftpos : leftpos + 1]
            loss = F.huber_loss(prediction, targetdelta)
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.head.parameters(), 1.0)
            optimizer.step()
            losses.append(float(loss.detach().cpu()))
        history.append({"epoch": float(epoch + 1), "difference": float(np.mean(losses))})
    return history


def fitmodel(
    model: EncodedAnchorDelta,
    frame: pd.DataFrame,
    proteins: dict[str, torch.Tensor],
    feature: np.ndarray,
    pairs: list[tuple[str, int, int]],
    epochs: int,
    seed: int,
) -> list[dict[str, float]]:
    """Retrain interaction encoder and comparator on target-balanced pairs."""

    trainable = [parameter for parameter in model.parameters() if parameter.requires_grad]
    optimizer = torch.optim.AdamW(trainable, lr=2e-4, weight_decay=1e-4)
    bytarget: dict[str, list[tuple[int, int]]] = {}
    for target, left, right in pairs:
        bytarget.setdefault(target, []).append((left, right))
    generator = np.random.default_rng(seed)
    history: list[dict[str, float]] = []
    for epoch in range(epochs):
        model.train()
        losses: list[float] = []
        targets = list(bytarget)
        generator.shuffle(targets)
        for target in targets:
            targetpairs = bytarget[target]
            used = sorted({index for pair in targetpairs for index in pair})
            positions = {index: position for position, index in enumerate(used)}
            sources = frame.source_row.iloc[used].to_numpy(dtype=np.int64)
            ligands = torch.as_tensor(feature[sources], device="cuda", dtype=torch.float32)
            encoded = model.encode(proteins[target], ligands)
            order = generator.permutation(len(targetpairs))
            leftpositions = []
            rightpositions = []
            targetdeltas = []
            for position in order:
                left, right = targetpairs[int(position)]
                if generator.random() < 0.5:
                    left, right = right, left
                leftpositions.append(positions[left])
                rightpositions.append(positions[right])
                targetdeltas.append(
                    float(frame.affinity.iat[right] - frame.affinity.iat[left])
                )
            leftindex = torch.as_tensor(leftpositions, device="cuda", dtype=torch.long)
            rightindex = torch.as_tensor(rightpositions, device="cuda", dtype=torch.long)
            prediction = model.delta(encoded[rightindex], encoded[leftindex])
            targetdelta = torch.as_tensor(targetdeltas, device="cuda", dtype=torch.float32)
            loss = F.huber_loss(prediction, targetdelta)
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(trainable, 1.0)
            optimizer.step()
            losses.append(float(loss.detach().cpu()))
        history.append({"epoch": float(epoch + 1), "difference": float(np.mean(losses))})
    return history


@torch.no_grad()
def evaluate(
    model: EncodedAnchorDelta,
    frame: pd.DataFrame,
    episodes: list,
    proteins: dict[str, torch.Tensor],
    feature: np.ndarray,
    base: np.ndarray,
    wrong: dict[str, str],
    wronglabels: dict[str, torch.Tensor],
) -> dict[str, object]:
    model.eval()
    predictions = {
        name: []
        for name in (
            "anchordelta",
            "calibration",
            "permutedlabels",
            "wronglabels",
            "wrongprotein",
        )
    }
    labels: list[float] = []
    indices: list[int] = []
    for episode in episodes:
        support = list(episode.support_indices)
        query = list(episode.query_indices)
        supportsource = frame.source_row.iloc[support].to_numpy(dtype=np.int64)
        querysource = frame.source_row.iloc[query].to_numpy(dtype=np.int64)
        sx = torch.as_tensor(feature[supportsource], device="cuda")
        qx = torch.as_tensor(feature[querysource], device="cuda")
        sy = torch.as_tensor(frame.affinity.iloc[support].to_numpy(dtype=np.float32), device="cuda")
        qy = torch.as_tensor(frame.affinity.iloc[query].to_numpy(dtype=np.float32), device="cuda")
        # ``fitbase`` returns predictions aligned to the gate frame, not the
        # global ligand feature array.
        sb = torch.as_tensor(base[support], device="cuda")
        qb = torch.as_tensor(base[query], device="cuda")
        delta = model(proteins[episode.target_key], qx, sx)
        anchor, _ = anchorabsolute(sy, delta)
        permuted, _ = anchorabsolute(sy.roll(1), delta)
        wronglabelvalue, _ = anchorabsolute(wronglabels[episode.target_key], delta)
        wrongdelta = model(proteins[wrong[episode.target_key]], qx, sx)
        wrongvalue, _ = anchorabsolute(sy, wrongdelta)
        calibration = qb + (sy - sb).mean()
        predictions["anchordelta"].extend(anchor.cpu().tolist())
        predictions["calibration"].extend(calibration.cpu().tolist())
        predictions["permutedlabels"].extend(permuted.cpu().tolist())
        predictions["wronglabels"].extend(wronglabelvalue.cpu().tolist())
        predictions["wrongprotein"].extend(wrongvalue.cpu().tolist())
        labels.extend(qy.cpu().tolist())
        indices.extend(query)
    components = {episode.target_key: episode.homology_component for episode in episodes}
    metrics = {
        name: evaluateprotocol(
            predictions=values,
            labels=labels,
            episodes=episodes,
            prediction_indices=indices,
            component_by_target=components,
        )
        for name, values in predictions.items()
    }
    paired = pairedcomponents(
        predictions=predictions,
        labels=labels,
        episodes=episodes,
        prediction_indices=indices,
        component_by_target=components,
        reference="anchordelta",
    )
    return {"metrics": metrics, "paired_component_bootstrap": paired, "episodes": len(episodes)}


def makewronglabels(frame: pd.DataFrame, episodes: list, seed: int) -> dict[str, torch.Tensor]:
    """Rotate complete support-label vectors between targets, preserving k."""

    targets = sorted(episode.target_key for episode in episodes)
    if len(targets) < 2:
        raise ValueError("wrong-support control requires at least two targets")
    generator = np.random.default_rng(seed)
    order = np.asarray(targets, dtype=object)
    generator.shuffle(order)
    values = {
        target: torch.as_tensor(
            frame.affinity.iloc[list(next(item for item in episodes if item.target_key == source).support_indices)]
            .to_numpy(dtype=np.float32),
            device="cuda",
        )
        for target, source in zip(targets, np.roll(order, -1))
    }
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("dataset/public/chembl_37/processed/dualcold"))
    parser.add_argument(
        "--checkpoint",
        type=str,
        default="none",
        help="fit-only checkpoint with a recorded target/component manifest; default is a fresh encoder",
    )
    parser.add_argument("--fit-targets", type=int, default=32)
    parser.add_argument("--gate-targets", type=int, default=16)
    parser.add_argument("--pairs-per-target", type=int, default=32)
    parser.add_argument("--epochs", type=int, default=1)
    parser.add_argument("--train-encoder", action="store_true")
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--out", type=Path, default=REPORT)
    args = parser.parse_args()
    if not torch.cuda.is_available():
        raise RuntimeError("AnchorDelta P0 requires CUDA")
    if args.gate_targets < 2:
        raise ValueError("AnchorDelta wrong-support control requires at least two gate targets")
    if args.fit_targets < 1:
        raise ValueError("AnchorDelta requires at least one fit target")
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    full = loadlabels(args.root, "pKi")
    train = full.loc[full.dual_cold_split == "train"].reset_index(drop=True)
    fit, gate = splitcomponents(train, 0.2, args.seed)
    fit = limitrows(fit, args.fit_targets, 0, args.seed + 1)
    gate = limitrows(gate, args.gate_targets, 0, args.seed + 2)
    fit = fit.reset_index(drop=True)
    gate = gate.reset_index(drop=True)
    rawfeature = np.load(args.root / "ligand_features.npz", allow_pickle=False)["feat"]
    feature, normalization = normalizeligands(rawfeature, fit.source_row.to_numpy(), descriptors=10)
    proteins = loadproteins(args.root)
    _, basegate, _ = fitbase(fit, gate, feature)
    checkpoint = None if args.checkpoint.lower() == "none" else Path(args.checkpoint)
    if checkpoint is not None:
        raise ValueError(
            "uncertified checkpoints are disabled; provide a fit-only checkpoint manifest or use --checkpoint none"
        )
    if not args.train_encoder:
        raise ValueError("a fresh encoder must be trained on fit components; pass --train-encoder")
    interaction = buildinteraction(checkpoint).interaction
    model = EncodedAnchorDelta(
        interaction,
        feature_dim=8,
        hidden_dim=64,
        freeze_encoder=not args.train_encoder,
    ).cuda()
    pairs = samplepairs(fit, args.pairs_per_target, args.seed)
    if args.train_encoder:
        history = fitmodel(model, fit, proteins, feature, pairs, args.epochs, args.seed)
    else:
        fitencoded = encodeframe(model, fit, proteins, feature)
        history = fithead(model, fitencoded, pairs, args.epochs, args.seed)
    episodes = maketrainroster(gate, proteins, feature, targets=args.gate_targets, querycap=64, support=5)
    result = evaluate(
        model,
        gate,
        episodes,
        proteins,
        feature,
        basegate,
        wrongtargets(gate, args.seed + 9),
        makewronglabels(gate, episodes, args.seed + 10),
    )
    result.update({
        "protocol": "P0 leakage-safe fit-component antisymmetric AnchorDelta",
        "seed": args.seed,
        "fit_targets": len(fit.target.unique()),
        "gate_targets": len(gate.target.unique()),
        "pairs": len(pairs),
        "history": history,
        "normalization": normalization,
        "checkpoint": None if checkpoint is None else str(checkpoint),
        "train_component_holdout": True,
        "trainable_encoder": args.train_encoder,
        "checkpoint_provenance": "fresh random initialization trained only on fit components",
    })
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=True), encoding="utf-8")
    torch.save(model.state_dict(), args.out.with_suffix(".pt"))
    print(json.dumps({"metrics": result["metrics"], "out": str(args.out)}, indent=2))


if __name__ == "__main__":
    main()

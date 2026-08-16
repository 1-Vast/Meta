"""Run a small, auditable AdaMBind meta-learning smoke on downloaded data."""

from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import torch


def _metrics(labels, preds):
    labels = np.asarray(labels, dtype=float)
    preds = np.asarray(preds, dtype=float)
    finite = np.isfinite(labels) & np.isfinite(preds)
    labels, preds = labels[finite], preds[finite]
    if len(labels) == 0:
        return {"n": 0, "mse": None, "rmse": None, "mae": None}
    return {
        "n": int(len(labels)),
        "mse": float(np.mean((labels - preds) ** 2)),
        "rmse": float(np.sqrt(np.mean((labels - preds) ** 2))),
        "mae": float(np.mean(np.abs(labels - preds))),
    }


def _task_split(frame, seed, train_tasks, val_tasks, test_tasks, min_rows):
    counts = frame.groupby("target_sequence").size()
    targets = sorted(str(t) for t, n in counts.items() if n >= min_rows)
    rng = np.random.default_rng(seed)
    targets = list(np.asarray(targets, dtype=object)[rng.permutation(len(targets))])
    required = train_tasks + val_tasks + test_tasks
    if len(targets) < required:
        raise ValueError(f"only {len(targets)} targets have >= {min_rows} rows; need {required}")
    return (
        [str(x) for x in targets[:train_tasks]],
        [str(x) for x in targets[train_tasks : train_tasks + val_tasks]],
        [str(x) for x in targets[train_tasks + val_tasks : required]],
    )


def run(args):
    source = Path(args.source).resolve()
    data_root = Path(args.data_root).resolve()
    sys.path.insert(0, str(source))
    from model.Trainer import Trainer
    from model.gat_gcn import GAT_GCN
    from model.scheduler import Scheduler
    from utils.TestbedDataset import TestbedDataset

    csv_path = data_root / f"{args.dataset}-full-data.csv"
    processed_root = data_root
    if not csv_path.exists():
        raise FileNotFoundError(csv_path)
    frame = pd.read_csv(csv_path)
    for col in ("compound_iso_smiles", "target_sequence", "affinity"):
        if col not in frame.columns:
            raise ValueError(f"missing column {col}")
    if not np.isfinite(pd.to_numeric(frame["affinity"], errors="coerce")).all():
        raise ValueError("non-finite affinity values cannot be used in this smoke")

    train_targets, val_targets, test_targets = _task_split(
        frame, args.seed, args.train_tasks, args.val_tasks, args.test_tasks, args.support + args.query
    )
    dataset = TestbedDataset(root=str(processed_root), dataset=f"{args.dataset}-full-data")
    if len(dataset) != len(frame):
        raise ValueError(f"processed rows {len(dataset)} != CSV rows {len(frame)}")

    groups = {}
    for target, group in frame.groupby("target_sequence", sort=True):
        indices = list(group.index.astype(int))
        if len(indices) >= args.support + args.query:
            groups[str(target)] = [dataset[i] for i in indices]

    def build_fdata(targets):
        out = {}
        for target in targets:
            rows = groups[target]
            # Deterministic within-target support/query selection mirrors
            # DataSplit.py but avoids a second random stream in the smoke.
            out[target] = [rows[: args.support], rows[args.support : args.support + args.query]]
        return out

    fdata = build_fdata(train_targets + val_targets + test_targets)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    if device.type == "cuda":
        torch.cuda.manual_seed_all(args.seed)
        torch.cuda.reset_peak_memory_stats()

    net = GAT_GCN().to(device)
    n_params = len(list(net.parameters()))
    scheduler = Scheduler(n_params, args.meta_batch_size, list(range(n_params))).to(device)
    scheduler_optimizer = torch.optim.Adam(scheduler.parameters(), lr=1e-4)
    trainer = Trainer(net)
    train_idx = np.asarray(train_targets, dtype=object)
    val_idx = np.asarray(val_targets, dtype=object)
    task_args = SimpleNamespace(
        batch_size=args.batch_size,
        noise=0,
        noise_val=0.0,
        reg_lr=args.reg_lr,
        meta_lr=args.meta_lr,
        update_step_train=args.update_step_train,
        update_step_test=args.update_step_test,
    )

    started = time.perf_counter()
    fast_weights = trainer.train(net, task_args, 0, train_idx, fdata, update=0)
    losses = trainer.get_rloss()
    grad1, grad2 = trainer.get_grads()
    _, _, scores = scheduler.ats(losses, grad1, grad2, 0)
    probabilities = torch.softmax(scores.reshape(-1), dim=-1)
    selected = scheduler.sample_task(probabilities, args.meta_batch_size)
    selected_targets = train_idx[selected].tolist()
    fast_weights = trainer.train(net, task_args, 0, np.asarray(selected_targets), fdata, update=0)
    preds, labels, val_loss, val_ci, val_r2, val_spear, val_pear = trainer.predict(
        net, task_args, val_idx, fdata, fast_weights
    )

    # Exercise the policy-gradient update used by AdaMBind.  The returned
    # score is diagnostic only; no checkpoint is retained from this smoke.
    policy_loss = 0.0
    for index in selected:
        policy_loss = policy_loss - scheduler.m.log_prob(torch.tensor(index, device=device))
    policy_loss = policy_loss * (-float(val_loss))
    scheduler_optimizer.zero_grad()
    policy_loss.backward()
    scheduler_optimizer.step()

    result = {
        "status": "ok",
        "dataset": args.dataset,
        "csv": str(csv_path),
        "processed": str(processed_root / "processed" / f"{args.dataset}-full-data.pt"),
        "seed": args.seed,
        "device": str(device),
        "cuda_available": bool(torch.cuda.is_available()),
        "cuda_device": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
        "rows": int(len(frame)),
        "targets": {"train": train_targets, "val": val_targets, "test": test_targets},
        "support": args.support,
        "query": args.query,
        "selected_targets": selected_targets,
        "selected_task_probabilities": [float(probabilities[i].detach().cpu()) for i in selected],
        "trainer": {
            "n_parameter_tensors": n_params,
            "update_step_train": args.update_step_train,
            "update_step_test": args.update_step_test,
            "train_query_losses": [float(x.detach().cpu()) for x in losses],
        },
        "validation": {
            "mse": float(val_loss),
            "ci": float(val_ci),
            "r2": float(val_r2),
            "spearman": float(val_spear),
            "pearson": float(val_pear),
            "finite_metrics": _metrics(labels, preds),
        },
        "wall_seconds": float(time.perf_counter() - started),
        "peak_cuda_allocated_bytes": int(torch.cuda.max_memory_allocated())
        if device.type == "cuda"
        else 0,
        "scientific_use": "mechanism smoke only; not a strict FORT cold-split performance claim",
    }
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--dataset", default="davis")
    parser.add_argument("--seed", type=int, default=168)
    parser.add_argument("--train-tasks", type=int, default=2)
    parser.add_argument("--val-tasks", type=int, default=1)
    parser.add_argument("--test-tasks", type=int, default=1)
    parser.add_argument("--support", type=int, default=3)
    parser.add_argument("--query", type=int, default=8)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--update-step-train", type=int, default=1)
    parser.add_argument("--update-step-test", type=int, default=1)
    parser.add_argument("--reg-lr", type=float, default=1e-4)
    parser.add_argument("--meta-lr", type=float, default=1e-5)
    parser.add_argument("--meta-batch-size", type=int, default=1)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    command = " ".join(sys.argv)
    try:
        result = run(args)
    except Exception as exc:
        result = {
            "status": "failed",
            "command": command,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(),
        }
        args.out.write_text(json.dumps(result, indent=2), encoding="utf-8")
        print(json.dumps(result, indent=2))
        raise
    result["command"] = command
    args.out.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

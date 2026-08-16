"""Controlled synthetic capacity probe for the zero-shot interaction trunk.

Both architectures are trained, at matched budget, on a synthetic task whose
label is a protein-by-ligand contact-type bilinear form.  The probe answers a
single question the real corpus cannot answer cheaply: can the trunk express a
protein-conditioned interaction at all, given enough gradient steps?

It is a falsification instrument, not a performance claim.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

import torch

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM
from scripts.train_qpsmp import resolve_architecture


PROTEIN_DIM = 32
RESIDUES = 8
ATOMS = 8
QUERIES = 6


def corpus(generator: torch.Generator, tasks: int, types: int,
           codebook: torch.Tensor, weights: torch.Tensor):
    atom_type = torch.randint(types, (tasks, QUERIES, ATOMS), generator=generator)
    residue_type = torch.randint(types, (tasks, RESIDUES), generator=generator)
    raw_atoms = torch.zeros(tasks, QUERIES, ATOMS, ATOM_FEAT_DIM)
    raw_atoms.scatter_(-1, atom_type.unsqueeze(-1), 1.0)
    tokens = codebook[residue_type]
    ligand_share = torch.nn.functional.one_hot(atom_type, types).float().mean(-2)
    pocket_share = torch.nn.functional.one_hot(residue_type, types).float().mean(-2)
    truth = (ligand_share * pocket_share[:, None, :] * weights).sum(-1)
    bonds = torch.zeros(tasks, QUERIES, ATOMS, ATOMS, BOND_FEAT_DIM)
    bonds[..., 0] = 1.0
    episode = {
        "protein_pooled": tokens.mean(1), "protein_tokens": tokens,
        "protein_mask": torch.ones(tasks, RESIDUES),
        "protein_chemistry": torch.zeros(tasks, RESIDUES, 4),
        "support_atoms": raw_atoms[:, :0], "support_bonds": bonds[:, :0],
        "support_mask": torch.ones(tasks, 0, ATOMS),
        "support_y": torch.zeros(tasks, 0),
        "query_atoms": raw_atoms, "query_bonds": bonds,
        "query_mask": torch.ones(tasks, QUERIES, ATOMS),
    }
    return episode, truth


def run(model, episode, adapt=False):
    return model(
        episode["protein_pooled"], episode["protein_tokens"],
        episode["protein_mask"], episode["support_atoms"],
        episode["support_bonds"], episode["support_mask"], episode["support_y"],
        episode["query_atoms"], episode["query_bonds"], episode["query_mask"],
        adapt=adapt, protein_chemistry=episode["protein_chemistry"])


def probe(arch: str, steps: int, seed: int, device: str,
          learning_rate: float = 3e-3) -> dict:
    torch.manual_seed(seed)
    generator = torch.Generator().manual_seed(seed + 1)
    types = 6
    codebook = torch.randn(types, PROTEIN_DIM, generator=generator)
    weights = torch.randn(types, generator=generator)
    model = resolve_architecture(arch)(
        protein_dim=PROTEIN_DIM, hidden_dim=64, task_dim=16, ligand_layers=3,
        pair_dim=48, pair_blocks=3, pair_latents=12, pair_heads=4,
        pair_chunk_size=32).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=learning_rate)
    started = time.monotonic()
    trace = []
    for step in range(steps):
        episode, truth = corpus(generator, 16, types, codebook, weights)
        episode = {k: v.to(device) for k, v in episode.items()}
        out = run(model, {k: v for k, v in episode.items()})
        loss = (out.prediction - truth.to(device)).square().mean()
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        if (step + 1) % max(steps // 8, 1) == 0:
            trace.append({"step": step + 1, "loss": float(loss.detach())})
    heldout, truth = corpus(
        torch.Generator().manual_seed(73201), 256, types, codebook, weights)
    heldout = {k: v.to(device) for k, v in heldout.items()}
    truth = truth.to(device)
    with torch.no_grad():
        out = run(model, heldout)
        order = torch.randperm(truth.shape[0], device=device)
        shuffled = dict(heldout)
        for key in ("protein_pooled", "protein_tokens", "protein_mask",
                    "protein_chemistry"):
            shuffled[key] = heldout[key].index_select(0, order)
        wrong = run(model, shuffled)
    variance = float(truth.var())
    return {
        "arch": arch,
        "learning_rate": learning_rate,
        "trainable_parameters": int(sum(
            p.numel() for p in model.parameters() if p.requires_grad)),
        "steps": steps,
        "seconds": time.monotonic() - started,
        "label_variance": variance,
        "heldout_mse": float((out.prediction - truth).square().mean()),
        "heldout_mse_relative": float(
            (out.prediction - truth).square().mean()) / variance,
        "shuffled_protein_mse": float((wrong.prediction - truth).square().mean()),
        "prediction_spread_within_task": float(out.prediction.std(-1).mean()),
        "truth_spread_within_task": float(truth.std(-1).mean()),
        "loss_trace": trace,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--steps", type=int, default=1500)
    parser.add_argument("--seed", type=int, default=9101)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--archs", nargs="+", default=["bpsf", "grammar"])
    parser.add_argument("--learning-rates", type=float, nargs="+",
                        default=[3e-3])
    args = parser.parse_args()
    rows = [probe(arch, args.steps, args.seed, args.device, lr)
            for arch in args.archs for lr in args.learning_rates]
    for row in rows:
        print(json.dumps({k: v for k, v in row.items() if k != "loss_trace"},
                         indent=2, sort_keys=True))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(
        {"schema": "MetaSieve.QPSMPTrunkCapacityProbe.v1",
         "task": "protein-by-ligand contact-type bilinear form, zero-shot only",
         "device": args.device, "probes": rows},
        indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

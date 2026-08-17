"""Where does A0's protein-conditioning of the ordering disappear?

`branch_ordering_probe.py` establishes that the trained interaction branch
orders ligands almost identically under a wrong protein (centered shift
0.0007 pK) while a randomly initialised copy of the *same architecture* shifts
by 0.0766 pK. Training removes the conditioning; the architecture can express
it.

This probe localizes the collapse along the one protein-conditioned path in
the trunk (`model/interaction_grammar.py`):

    residues  --ResidueEncoder-->  residue slots
    atoms     --LigandEncoder--->  atom states
              --ContactGrammar--->  attention weights  ->  (occupancy,
                                                            mean_state,
                                                            max_state)
              --embed/section---->  interaction readout

Two candidate loci, with different remedies:

* **attention**: the atom-to-residue attention weights themselves become
  protein-invariant. The protein stops being read at all, and a training
  signal has to restore the conditioning upstream.
* **readout**: attention still responds to the protein, but the downstream
  `embed`/`interaction_head` washes the difference out of the *ligand-
  differential* while keeping it in the level. Then the readout, or a
  normalization inside it, is the structural culprit.

Measured per target, for the correct protein and a matched wrong protein:

* Jensen-Shannon divergence between the two attention weight distributions,
  averaged over atoms and heads;
* the same divergence for the *ligand-differential* — how much the attention
  difference varies across the queries of one target, which is the part that
  can carry ordering;
* the centered spread of `occupancy`, `mean_state` and `max_state` shifts.

No training, no gradients. `meta_test` excluded logically after parsing.
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

from scripts.qpsmp_data import QPSMPData                          # noqa: E402
from scripts.stageR0_retrieval_falsification import (             # noqa: E402
    component_target_mean,
)
from scripts.stageR6_compare_arms import SUPPORT_SIZES            # noqa: E402
from tools.research.a2_readiness._arms import build_arm           # noqa: E402
from scripts.train_level_shape import matched_donors, normalized  # noqa: E402
from scripts.train_qpsmp import (                                 # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    compact_episode, training_label_scale,
)


def capture(model, data, target, episode, device, dtype) -> dict:
    """Run the trunk for one protein, capturing attention and grammar outputs."""
    store: dict = {}

    def hook(module, inputs, output):
        # ContactGrammar.forward returns (occupancy, mean_state, max_state).
        store["occupancy"], store["mean_state"], store["max_state"] = output

    handle = model.grammar.register_forward_hook(hook)

    # Recompute the attention weights explicitly: the module does not expose
    # them, and re-deriving them here keeps the probe read-only.
    original = model.grammar.forward

    def instrumented(atoms, atom_mask, atom_chemistry, residues, residue_mask):
        grammar = model.grammar
        biased = atoms + grammar.atom_chemistry(atom_chemistry)
        pairs, atom_count, hidden = biased.shape
        residue_count = residues.shape[1]
        query = grammar.atom_query(biased).reshape(
            pairs, atom_count, grammar.heads, grammar.head_dim).transpose(1, 2)
        key = grammar.residue_key(residues).reshape(
            pairs, residue_count, grammar.heads, grammar.head_dim).transpose(1, 2)
        score = query @ key.transpose(-1, -2) / grammar.head_dim ** 0.5
        score = score.masked_fill(
            ~residue_mask[:, None, None, :].bool(), torch.finfo(score.dtype).min)
        store["attention"] = torch.softmax(score, dim=-1).detach()
        store["atom_mask"] = atom_mask.detach()
        return original(atoms, atom_mask, atom_chemistry, residues, residue_mask)

    model.grammar.forward = instrumented
    try:
        pooled, tokens, mask = data.protein_for_target(target)
        chemistry = data.protein_chemistry_for_target(target)
        model(pooled.to(device, dtype).unsqueeze(0),
              tokens.to(device, dtype).unsqueeze(0),
              mask.to(device, dtype).unsqueeze(0),
              episode.query_atoms[:0].to(device, dtype).unsqueeze(0),
              episode.query_bonds[:0].to(device, dtype).unsqueeze(0),
              episode.query_mask[:0].to(device, dtype).unsqueeze(0),
              episode.query_y[:0].to(device, dtype).unsqueeze(0),
              episode.query_atoms.to(device, dtype).unsqueeze(0),
              episode.query_bonds.to(device, dtype).unsqueeze(0),
              episode.query_mask.to(device, dtype).unsqueeze(0),
              adapt=False,
              protein_chemistry=chemistry.to(device, dtype).unsqueeze(0),
              support_fingerprint=episode.query_fingerprint[:0].to(device, dtype).unsqueeze(0),
              query_fingerprint=episode.query_fingerprint.to(device, dtype).unsqueeze(0))
    finally:
        model.grammar.forward = original
        handle.remove()
    return store


def jensen_shannon(left: torch.Tensor, right: torch.Tensor,
                   mask: torch.Tensor) -> float:
    """Mean JS divergence over valid atoms and heads, in nats."""
    middle = 0.5 * (left + right)

    def kl(p, q):
        return (p * ((p.clamp_min(1e-12)).log() - q.clamp_min(1e-12).log())).sum(-1)

    divergence = 0.5 * kl(left, middle) + 0.5 * kl(right, middle)
    weight = mask[:, None, :].expand_as(divergence)
    total = weight.sum().clamp_min(1.0)
    return float((divergence * weight).sum() / total)


def centered_spread(correct: torch.Tensor, wrong: torch.Tensor) -> float:
    """Spread across queries of the correct-minus-wrong difference."""
    difference = (correct - wrong).float()
    difference = difference - difference.mean(0, keepdim=True)
    return float(difference.std())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", action="append", required=True)
    parser.add_argument("--split-directory", type=Path, required=True)
    parser.add_argument("--evaluation-seed", type=int, default=73101)
    parser.add_argument("--query-size", type=int, default=16)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=arguments.split_directory,
                     include_meta_test=False)
    scale = training_label_scale(data)
    donors = matched_donors(data, "meta_val", donor_pool="meta_val",
                            whitening_pool="meta_train")
    specs = data.fixed_nested_episode_banks(
        "meta_val", SUPPORT_SIZES, arguments.query_size, 1,
        arguments.evaluation_seed, None)[0]

    rows: list[dict] = []
    for item in arguments.arm:
        name, _, path = item.partition("=")
        model, _, seed = build_arm(path, name, data, arguments.device)
        dtype = next(model.parameters()).dtype
        with torch.no_grad():
            for spec in specs:
                episode = compact_episode(
                    normalized(data.materialize(spec), scale))
                correct = capture(model, data, spec.target, episode,
                                  arguments.device, dtype)
                wrong = capture(model, data, donors[spec.target], episode,
                                arguments.device, dtype)
                rows.append({
                    "arm": name, "seed": seed, "target": spec.target,
                    "component": spec.component,
                    "attention_js": jensen_shannon(
                        correct["attention"], wrong["attention"],
                        correct["atom_mask"]),
                    "occupancy_centered_shift": centered_spread(
                        correct["occupancy"], wrong["occupancy"]),
                    "mean_state_centered_shift": centered_spread(
                        correct["mean_state"], wrong["mean_state"]),
                    "max_state_centered_shift": centered_spread(
                        correct["max_state"], wrong["max_state"]),
                })
        del model
        if arguments.device.startswith("cuda"):
            torch.cuda.empty_cache()

    fields = ("attention_js", "occupancy_centered_shift",
              "mean_state_centered_shift", "max_state_centered_shift")
    report = {
        "schema": "MetaSieve.A2ReadinessAttentionLocus.v1",
        "split": "meta_val", "k": 0,
        "split_assignment_sha256": data.split_manifest["assignment_sha256"],
        "meta_test": {"evaluated": False, "included": False,
                      "seal": "physical: QPSMPData(include_meta_test=False)"},
        "arms": {
            name: {field: component_target_mean(
                (r["component"], r["target"], r[field])
                for r in rows if r["arm"] == name) for field in fields}
            for name in sorted({r["arm"] for r in rows})},
    }
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(report, indent=1), encoding="utf-8")

    print(f"{'arm':<10}{'attn JS':>10}{'occupancy':>12}{'mean_state':>12}{'max_state':>12}")
    for name, cell in report["arms"].items():
        print(f"{name:<10}{cell['attention_js']:>10.5f}"
              f"{cell['occupancy_centered_shift']:>12.5f}"
              f"{cell['mean_state_centered_shift']:>12.5f}"
              f"{cell['max_state_centered_shift']:>12.5f}")
    print(f"\nwrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

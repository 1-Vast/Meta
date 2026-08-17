"""Phase 3: where in the trunk is the protein's ligand-differential lost?

Attention-weight divergence is not evidence. A randomly initialised model's
attention also changes under input substitution, and the v1 probe's headline —
"trained attention is as protein-sensitive as random (JS 0.241 vs 0.218), so
the collapse is downstream in the readout" — is an inference from a
correlational statistic, not a causal measurement. Phase 2 already contradicts
it: the ligand-differential is protein-invariant at `occupancy` and
`mean_state` (cosine 1.0000), which are the *immediate outputs* of the
attention block, with no readout in between.

This audit settles it by intervention. The protein reaches the ligand-varying
path through exactly two channels inside `ContactGrammar`:

    weight  = softmax( atom_query(atoms) @ residue_key(residues)^T / sqrt(d) )
    context = weight @ residue_value(residues)

* the **routing** channel — `residue_key`, which sets *where* each atom looks;
* the **content** channel — `residue_value`, which sets *what it retrieves*.

Both are then fused with the ligand's own atom states:

    state = atom_context( [atoms, context, atoms * context] )

so a third possibility is that both channels deliver protein information and
the *fusion* discards it.

The audit runs each channel separately, holding the other at the correct
protein, and measures the change in the ligand-differential (the centered part,
which is all that can order ligands) at every stage. A channel that carries no
protein-differential cannot be repaired by changing the readout; a channel that
carries one which then vanishes downstream localises the loss to the fusion.

Gradient magnitudes complete the picture: the Jacobian of the *level* and of
the *ligand-differential* with respect to the protein tokens, measured
separately. A trunk whose level-gradient is large and differential-gradient is
small is not failing to see the protein — it is routing it entirely into the
level.

No training. No claim here describes attention as pocket-aware, contact-
resolved or biologically localized: the protein path is provably invariant to
residue-slot order (`tests/test_probe_structure.py`), so the attention cannot
be reading ordered structure of any kind.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import torch
import torch.nn.functional as F

ROOT = Path(__file__).resolve().parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.qpsmp_data import QPSMPData                              # noqa: E402
from scripts.stageR0_retrieval_falsification import (                 # noqa: E402
    component_bootstrap, component_target_mean,
)
from scripts.train_level_shape import normalized                      # noqa: E402
from scripts.train_qpsmp import (                                     # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    compact_episode, training_label_scale,
)
from tools.research.a2_readiness_v2 import _frozen                    # noqa: E402
from tools.research.a2_readiness_v2._arms import random_arm, trained_arm  # noqa: E402
from tools.research.a2_readiness_v2._donors import stratified_donors  # noqa: E402

CHANNELS = ("both", "routing_only", "content_only", "none")


def protein_inputs(data, target: str, device: str, dtype):
    pooled, tokens, mask = data.protein_for_target(target)
    chemistry = data.protein_chemistry_for_target(target)
    return [pooled.to(device, dtype).unsqueeze(0),
            tokens.to(device, dtype).unsqueeze(0),
            mask.to(device, dtype).unsqueeze(0),
            chemistry.to(device, dtype).unsqueeze(0)]


def encode_residues(model, parts: list):
    """`ResidueEncoder` output for one protein: residue slots and the summary."""
    pooled, tokens, mask, chemistry = parts
    residues, summary = model.protein_encoder(pooled, tokens, mask, chemistry)
    return model.refine_slots(residues, mask), summary, mask


def grammar_split(model, atom_states, atom_mask, atom_chemistry,
                  key_residues, value_residues, residue_mask):
    """`ContactGrammar.forward` with the two protein channels separable.

    Reproduces the module exactly (verified by
    `tests/test_attention_intervention.py`) except that the residues feeding
    `residue_key` and `residue_value` may come from different proteins.
    """
    grammar = model.grammar
    pairs, atom_count, hidden = atom_states.shape
    residue_count = key_residues.shape[1]
    atoms = atom_states + grammar.atom_chemistry(atom_chemistry)
    query = grammar.atom_query(atoms).reshape(
        pairs, atom_count, grammar.heads, grammar.head_dim).transpose(1, 2)
    key = grammar.residue_key(key_residues).reshape(
        pairs, residue_count, grammar.heads, grammar.head_dim).transpose(1, 2)
    value = grammar.residue_value(value_residues).reshape(
        pairs, residue_count, grammar.heads, grammar.head_dim).transpose(1, 2)
    score = query @ key.transpose(-1, -2) / grammar.head_dim ** 0.5
    score = score.masked_fill(
        ~residue_mask[:, None, None, :].bool(), torch.finfo(score.dtype).min)
    weight = torch.softmax(score, dim=-1)
    context = (weight @ value).transpose(1, 2).reshape(pairs, atom_count, hidden)
    state = grammar.atom_context(
        torch.cat((atoms, context, atoms * context), dim=-1))
    state = state * atom_mask.unsqueeze(-1)
    occupancy = (torch.sigmoid(grammar.type_logit(state))
                 * F.softplus(grammar.type_strength(state)))
    occupancy = (occupancy * atom_mask.unsqueeze(-1)).sum(1)
    occupancy = occupancy / (1.0 + occupancy)
    denominator = atom_mask.sum(1, keepdim=True).clamp_min(1.0)
    mean_state = (state * atom_mask.unsqueeze(-1)).sum(1) / denominator
    max_state = state.masked_fill(
        atom_mask.unsqueeze(-1) == 0, torch.finfo(state.dtype).min).amax(1)
    max_state = torch.nan_to_num(max_state, neginf=0.0)
    return {"weight": weight, "context": context, "occupancy": occupancy,
            "mean_state": mean_state, "max_state": max_state}


def downstream(model, ligand, parts, summary, occupancy, mean_state, max_state):
    """`encode`'s tail: embed -> section -> the three endpoint branches."""
    count = mean_state.shape[0]
    wide_summary = summary.expand(count, -1)
    embed = model.embed_norm(model.embed(torch.cat(
        (ligand, mean_state, max_state, wide_summary, occupancy), -1)))
    section = model.section_norm(model.section(embed))
    interaction = (model.interaction_head(torch.cat((embed, section), -1)).squeeze(-1)
                   + model.contact_weight(occupancy).squeeze(-1))
    return {"embed": embed, "section": section,
            "interaction": interaction[:, None]}


def run(model, data, spec, episode, key_target, value_target, device, dtype):
    """One forward with independently chosen routing and content proteins."""
    atoms = episode.query_atoms.to(device, dtype)
    bonds = episode.query_bonds.to(device, dtype)
    mask = episode.query_mask.to(device, dtype)
    ligand, atom_states = model.ligand_encoder(atoms, bonds, mask)
    key_parts = protein_inputs(data, key_target, device, dtype)
    key_residues, key_summary, residue_mask = encode_residues(model, key_parts)
    if value_target == key_target:
        value_residues = key_residues
    else:
        value_residues, _, _ = encode_residues(
            model, protein_inputs(data, value_target, device, dtype))
    count = atoms.shape[0]
    wide_key = key_residues.expand(count, -1, -1)
    wide_value = value_residues.expand(count, -1, -1)
    wide_mask = residue_mask.expand(count, -1)
    stages = grammar_split(
        model, atom_states, mask, model.atom_chemistry(atoms),
        wide_key, wide_value, wide_mask)
    stages.update(downstream(
        model, ligand, key_parts, key_summary,
        stages["occupancy"], stages["mean_state"], stages["max_state"]))
    return stages


def differential(values: torch.Tensor) -> np.ndarray:
    """The within-target centered part: everything that can order ligands."""
    flat = values.reshape(values.shape[0], -1).detach().float().cpu().numpy()
    return flat - flat.mean(0, keepdims=True)


def relative_change(correct: torch.Tensor, other: torch.Tensor) -> float:
    left, right = differential(correct), differential(other)
    norm = float(np.linalg.norm(left))
    return float(np.linalg.norm(left - right) / norm) if norm > 1e-12 else float("nan")


def level_change(correct: torch.Tensor, other: torch.Tensor) -> float:
    left = correct.reshape(correct.shape[0], -1).detach().float().cpu().numpy()
    right = other.reshape(other.shape[0], -1).detach().float().cpu().numpy()
    norm = float(np.linalg.norm(left.mean(0)))
    return (float(np.linalg.norm(left.mean(0) - right.mean(0)) / norm)
            if norm > 1e-12 else float("nan"))


def jacobians(model, data, spec, episode, device, dtype) -> dict:
    """Separate gradient magnitudes for the level and the ligand-differential.

    `d level / d P` and `d differential / d P` answer different questions. The
    first asks whether the protein reaches the output at all; the second asks
    whether it reaches the only part of the output that can order ligands.
    Their ratio is the quantity the whole cycle turns on.
    """
    parts = protein_inputs(data, spec.target, device, dtype)
    parts[1] = parts[1].clone().requires_grad_(True)
    tokens = parts[1]

    atoms = episode.query_atoms.to(device, dtype)
    bonds = episode.query_bonds.to(device, dtype)
    mask = episode.query_mask.to(device, dtype)
    ligand, atom_states = model.ligand_encoder(atoms, bonds, mask)
    residues, summary, residue_mask = encode_residues(model, parts)
    count = atoms.shape[0]
    stages = grammar_split(
        model, atom_states, mask, model.atom_chemistry(atoms),
        residues.expand(count, -1, -1), residues.expand(count, -1, -1),
        residue_mask.expand(count, -1))
    stages.update(downstream(model, ligand, parts, summary,
                             stages["occupancy"], stages["mean_state"],
                             stages["max_state"]))
    stages["protein_value"] = model.protein_head(
        summary.expand(count, -1)).squeeze(-1)[:, None]

    out: dict = {}
    for name in ("weight", "context", "occupancy", "mean_state", "max_state",
                 "embed", "section", "interaction", "protein_value"):
        value = stages[name].reshape(stages[name].shape[0], -1)
        for part, tensor in (("level", value.mean(0)),
                             ("differential", value - value.mean(0, keepdim=True))):
            scalar = tensor.square().sum()
            grad = torch.autograd.grad(scalar, tokens, retain_graph=True,
                                       allow_unused=True)[0]
            magnitude = 0.0 if grad is None else float(grad.norm())
            # Normalise by the quantity's own scale so widths are comparable.
            denominator = float(tensor.norm()) + 1e-12
            out[f"{name}__{part}"] = magnitude / denominator
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=_frozen.SPLIT_DIRECTORY)
    scale = training_label_scale(data)
    donors = stratified_donors(data, "meta_val", _frozen.DONOR_POOL,
                               _frozen.WHITENING_POOL)
    specs = data.fixed_nested_episode_banks(
        "meta_val", (0,), _frozen.QUERY_SIZE, 1, _frozen.EVALUATION_SEED,
        None)[0]
    episodes = {spec.target: compact_episode(normalized(data.materialize(spec), scale))
                for spec in specs}

    stage_names = ("weight", "context", "occupancy", "mean_state", "max_state",
                   "embed", "section", "interaction")
    rows: list[dict] = []
    for arm, seed in (("A0", None), ("randinit", _frozen.RANDOM_INIT_SEEDS[0])):
        path = _frozen.A0_CHECKPOINTS[0]
        model, _, _ = (trained_arm(path, data, arguments.device) if seed is None
                       else random_arm(path, data, arguments.device, seed))
        model.eval()
        dtype = next(model.parameters()).dtype
        for spec in specs:
            episode = episodes[spec.target]
            donor = donors[spec.target]["nearest"][0]
            with torch.no_grad():
                correct = run(model, data, spec, episode, spec.target,
                              spec.target, arguments.device, dtype)
                variants = {
                    "both": run(model, data, spec, episode, donor, donor,
                                arguments.device, dtype),
                    "routing_only": run(model, data, spec, episode, donor,
                                        spec.target, arguments.device, dtype),
                    "content_only": run(model, data, spec, episode, spec.target,
                                        donor, arguments.device, dtype),
                }
            row = {"arm": arm, "target": spec.target, "component": spec.component}
            for channel, values in variants.items():
                for name in stage_names:
                    row[f"{name}__{channel}"] = relative_change(
                        correct[name], values[name])
                    row[f"level__{name}__{channel}"] = level_change(
                        correct[name], values[name])
            row.update({f"jac__{k}": v for k, v in jacobians(
                model, data, spec, episode, arguments.device, dtype).items()})
            rows.append(row)
        del model
        if arguments.device.startswith("cuda"):
            torch.cuda.empty_cache()
        print(f"  measured {arm}")

    def weighted(subset, field):
        return component_target_mean(
            (r["component"], r["target"], r.get(field)) for r in subset)

    payload: dict = {
        "schema": "MetaSieve.A2ReadinessV2.AttentionCausalAudit.v1",
        "split": "meta_val", "k": 0,
        "frozen_design": _frozen.frozen_manifest(),
        "meta_test": data.seal_record(),
        "interventions": {
            "both": "routing and content from the donor protein",
            "routing_only": "residue_key from the donor, residue_value correct",
            "content_only": "residue_value from the donor, residue_key correct",
        },
        "arms": {},
    }
    for arm in ("A0", "randinit"):
        subset = [r for r in rows if r["arm"] == arm]
        cell: dict = {"stages": {}, "jacobian": {}}
        for name in stage_names:
            cell["stages"][name] = {
                channel: {
                    "differential_relative_change": weighted(
                        subset, f"{name}__{channel}"),
                    "level_relative_change": weighted(
                        subset, f"level__{name}__{channel}"),
                    "differential_ci": component_bootstrap(
                        [(r["component"], r["target"], r[f"{name}__{channel}"])
                         for r in subset],
                        _frozen.BOOTSTRAP_DRAWS, _frozen.BOOTSTRAP_SEED),
                }
                for channel in ("both", "routing_only", "content_only")}
        for key in sorted(k for k in subset[0] if k.startswith("jac__")):
            cell["jacobian"][key[len("jac__"):]] = weighted(subset, key)
        payload["arms"][arm] = cell

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    with arguments.output.with_suffix(".rows.jsonl").open(
            "w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row) + "\n")
    report(payload, stage_names)
    print(f"\nwrote {arguments.output}")
    return 0


def report(payload: dict, stage_names) -> None:
    for arm, cell in payload["arms"].items():
        print(f"\n=== {arm}: relative change in the LIGAND-DIFFERENTIAL "
              f"(level in brackets)")
        print(f"{'stage':<13}{'both':>20}{'routing only':>20}{'content only':>20}")
        for name in stage_names:
            block = cell["stages"][name]
            cells = []
            for channel in ("both", "routing_only", "content_only"):
                cells.append(
                    f"{block[channel]['differential_relative_change']:.4f}"
                    f" [{block[channel]['level_relative_change']:.4f}]")
            print(f"{name:<13}{cells[0]:>20}{cells[1]:>20}{cells[2]:>20}")
        print(f"\n  Jacobian wrt protein tokens, normalised "
              f"(level vs differential):")
        print(f"    {'stage':<16}{'level':>12}{'differential':>14}{'ratio':>10}")
        for name in ("weight", "context", "occupancy", "mean_state", "max_state",
                     "embed", "section", "interaction", "protein_value"):
            level = cell["jacobian"].get(f"{name}__level", float("nan"))
            diff = cell["jacobian"].get(f"{name}__differential", float("nan"))
            ratio = level / diff if diff and np.isfinite(diff) and diff > 0 else float("nan")
            print(f"    {name:<16}{level:>12.4e}{diff:>14.4e}{ratio:>10.2f}")


if __name__ == "__main__":
    raise SystemExit(main())

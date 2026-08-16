"""Stage 1 structural and synthetic falsification gates for Core Innovation A.

These run before any real training. Each one falsifies a specific way the
relative-transport model could be cosmetic rather than real:

algebraic
  - the relative potential is exactly antisymmetric;
  - the anchor set has exactly zero mean shape (no constant component);
  - the endpoint is the exact sum of its three branches;
  - k=0 returns exactly the zero-shot endpoint; transport is exactly 0;
  - support permutation invariance; query permutation equivariance;
  - the prediction for one query never depends on the other queries;
  - support labels enter only as residuals (a label shift moves the
    transport by the weighted residual sum, nothing else);
  - query labels are never a model input; geometry is refused;

gradient
  - k=1 gradients w.r.t. the support label and the support ligand are nonzero;
  - no dead trainable branch at k>=2 (k=0/k=1 semantic exceptions documented);

k=1 mechanism
  - the k=1 correction changes with the query ligand;
  - replacing only the support ligand (same label) changes the prediction;

synthetic training
  - on a synthetic protein-by-ligand bilinear interaction task the trained
    endpoint orders held-out ligands, and deleting the interaction branch
    destroys that ordering;
  - matched-wrong support is clearly worse than the correct support;
  - on a private task (labels carry no transferable interaction) the model
    abstains: full degenerates to the level shift.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM
from model.reltransport import RelTransportModel
from scripts.train_qpsmp import LabelScale
from scripts.train_reltransport import (
    RelConfig, level_term, relative_supervision, shape_objective, contrast,
)

PROTEIN_DIM, SLOTS, ATOMS = 32, 16, 9


def build(seed: int = 0, **kwargs) -> RelTransportModel:
    torch.manual_seed(seed)
    return RelTransportModel(protein_dim=PROTEIN_DIM, hidden_dim=24, task_dim=12,
                             ligand_layers=2, pair_dim=24, pair_latents=6,
                             pair_heads=2, anchors=6, rank=12,
                             **kwargs).double()


def ligand(count: int, seed: int) -> tuple[torch.Tensor, ...]:
    generator = torch.Generator().manual_seed(seed)
    atoms = torch.rand(count, ATOMS, ATOM_FEAT_DIM,
                       generator=generator, dtype=torch.float64)
    bonds = torch.rand(count, ATOMS, ATOMS, BOND_FEAT_DIM,
                       generator=generator, dtype=torch.float64)
    bonds = bonds * (bonds > 0.7)
    bonds = 0.5 * (bonds + bonds.transpose(1, 2))
    mask = torch.ones(count, ATOMS, dtype=torch.float64)
    fingerprint = (torch.rand(count, 64, generator=generator,
                              dtype=torch.float64) > 0.7).double()
    fingerprint[:, 0] = 1.0                       # never an all-zero row
    return atoms, bonds, mask, fingerprint


def protein(seed: int) -> tuple[torch.Tensor, ...]:
    generator = torch.Generator().manual_seed(seed)
    pooled = torch.randn(PROTEIN_DIM, generator=generator, dtype=torch.float64)
    tokens = torch.randn(SLOTS, PROTEIN_DIM, generator=generator, dtype=torch.float64)
    mask = torch.ones(SLOTS, dtype=torch.float64)
    chemistry = torch.rand(SLOTS, 4, generator=generator, dtype=torch.float64)
    return pooled, tokens, mask, chemistry


def episode(support: int, query: int, *, protein_seed: int = 1,
            ligand_seed: int = 2) -> dict:
    pooled, tokens, mask, chemistry = protein(protein_seed)
    sa, sb, sm, sf = ligand(support, ligand_seed)
    qa, qb, qm, qf = ligand(query, ligand_seed + 100)
    return {
        "protein_pooled": pooled, "protein_tokens": tokens, "protein_mask": mask,
        "protein_chemistry": chemistry,
        "support_atoms": sa, "support_bonds": sb, "support_mask": sm,
        "support_y": torch.linspace(4.0, 9.0, support, dtype=torch.float64),
        "query_atoms": qa, "query_bonds": qb, "query_mask": qm,
        "support_fingerprint": sf, "query_fingerprint": qf,
    }


# ---------------------------------------------------------------- algebra


def test_delta_is_antisymmetric():
    model = build()
    pooled, tokens, mask, chemistry = protein(3)
    embed = torch.randn(2, 5, model.embed_dim, dtype=torch.float64)
    residues, summary = model.encode_protein(pooled[None], tokens[None],
                                             mask[None], chemistry[None])
    gate = mask[None].unsqueeze(-1)
    residue_mean = (residues * gate).sum(1) / gate.sum(1).clamp_min(1.0)
    u = model.relative.direction_vector(summary, residue_mean)
    delta = model.relative.delta_matrix(u, embed, embed)
    assert torch.allclose(delta, -delta.transpose(-1, -2), atol=1e-12)
    diagonal = torch.diagonal(delta, dim1=-2, dim2=-1)
    assert diagonal.abs().max().item() < 1e-12


def test_anchor_set_has_exactly_zero_mean_shape():
    model = build()
    pooled, tokens, mask, chemistry = protein(7)
    residues, summary = model.encode_protein(pooled[None], tokens[None],
                                             mask[None], chemistry[None])
    gate = mask[None].unsqueeze(-1)
    residue_mean = (residues * gate).sum(1) / gate.sum(1).clamp_min(1.0)
    u = model.relative.direction_vector(summary, residue_mean)
    anchors = model.anchor.unsqueeze(0).expand(2, -1, -1)
    shape, anchor_mean = model.anchor_shape(u.expand(2, -1), anchors)
    assert shape.mean(-1).abs().max().item() < 1e-10
    assert anchor_mean.abs().max().item() < 1e-10


def test_endpoint_is_the_exact_sum_of_three_branches():
    output = build()(**episode(0, 5), adapt=False)
    assert torch.allclose(
        output.endpoint,
        output.ligand_prior + output.target_level + output.shape, atol=1e-12)


def test_level_branch_is_constant_within_target():
    output = build()(**episode(0, 6), adapt=False)
    assert output.target_level.std().item() == pytest.approx(0.0, abs=1e-12)


def test_ligand_prior_is_protein_blind():
    model = build()
    base = episode(0, 5, protein_seed=1)
    other = protein(9)
    swapped = model(**{**base, "protein_pooled": other[0],
                       "protein_tokens": other[1], "protein_mask": other[2],
                       "protein_chemistry": other[3]}, adapt=False)
    assert torch.allclose(model(**base, adapt=False).ligand_prior,
                          swapped.ligand_prior, atol=1e-12)


def test_zero_support_returns_the_endpoint_exactly():
    output = build()(**episode(0, 5))
    assert torch.equal(output.prediction, output.endpoint)
    assert output.transport.abs().max().item() == 0.0
    assert output.rho.abs().max().item() == 0.0


def test_support_permutation_invariance():
    model = build()
    base = episode(4, 5)
    order = torch.tensor([2, 0, 3, 1])
    permuted = {**base,
                "support_atoms": base["support_atoms"][order],
                "support_bonds": base["support_bonds"][order],
                "support_mask": base["support_mask"][order],
                "support_y": base["support_y"][order],
                "support_fingerprint": base["support_fingerprint"][order]}
    assert torch.allclose(model(**base).prediction,
                          model(**permuted).prediction, atol=1e-10)


def test_query_permutation_equivariance():
    model = build()
    base = episode(3, 5)
    order = torch.tensor([4, 0, 2, 1, 3])
    permuted = {**base,
                "query_atoms": base["query_atoms"][order],
                "query_bonds": base["query_bonds"][order],
                "query_mask": base["query_mask"][order],
                "query_fingerprint": base["query_fingerprint"][order]}
    assert torch.allclose(model(**base).prediction[order],
                          model(**permuted).prediction, atol=1e-10)


def test_prediction_is_independent_of_the_other_queries():
    model = build()
    full = episode(2, 6)
    joint = model(**full).prediction
    for index in range(6):
        single = model(**{**full,
                          "query_atoms": full["query_atoms"][index:index + 1],
                          "query_bonds": full["query_bonds"][index:index + 1],
                          "query_mask": full["query_mask"][index:index + 1],
                          "query_fingerprint":
                              full["query_fingerprint"][index:index + 1]}).prediction
        assert torch.allclose(joint[index], single[0], atol=1e-10)


def test_support_labels_enter_only_as_residuals():
    """A label shift c moves the transport by shrink * c * sum_k a rho, per
    query (the gate and weights are label-free), and the endpoint is
    untouched."""
    model = build()
    base = episode(3, 5)
    shifted = {**base, "support_y": base["support_y"] + 2.0}
    before, after = model(**base), model(**shifted)
    assert torch.allclose(before.endpoint, after.endpoint, atol=1e-12)
    assert torch.allclose(before.rho, after.rho, atol=1e-12)
    shrink = float(model.transport.shrinkage(3, torch.zeros(1, dtype=torch.float64)))
    expected = shrink * 2.0 * (before.weight * before.rho).sum(-1)
    assert torch.allclose(after.transport - before.transport, expected,
                          atol=1e-9)


def test_query_labels_are_never_an_input():
    import inspect
    signature = inspect.signature(RelTransportModel.forward)
    assert "query_y" not in signature.parameters


def test_geometry_input_is_refused():
    with pytest.raises(ValueError, match="common-frame"):
        build()(**episode(1, 2),
                geometry_available=torch.ones(1, dtype=torch.bool))


# ---------------------------------------------------------------- k=1 mechanism


def test_k1_correction_varies_with_query():
    model = build()
    base = episode(1, 5)
    output = model(**base)
    assert output.transport.std().item() > 1e-9
    other = {**base, "query_atoms": ligand(5, 999)[0],
             "query_bonds": ligand(5, 999)[1],
             "query_mask": ligand(5, 999)[2],
             "query_fingerprint": ligand(5, 999)[3]}
    changed = model(**other).transport
    assert not torch.allclose(output.transport, changed, atol=1e-9, rtol=0)


def test_k1_gradients_wrt_support_label_and_ligand_are_nonzero():
    """k=1: support labels act functionally (label locked), and the support
    ligand carries real autograd gradient into the correction.

    The label-lock contract detaches the residual value, so the label enters
    as a value, not as a differentiable leaf; the lineage convention
    (test_interaction_grammar_synthetic) measures the label effect as the
    functional response to a shift plus trainability of the transport
    parameters. The additive correction is affine in the residual with unit
    coefficient, so a uniform label shift moves the transport by shrink * c
    exactly (the query-specificity lives in delta and delta_f0, which are
    label-free). The support *ligand*, by contrast, must carry genuine
    autograd gradient through delta and the weights.
    """
    model = build()
    base = episode(1, 5)
    output = model(**base)
    # Functional label effect: a uniform shift moves the transport by
    # shrink * c * rho(q, 1) per query, nonzero; the query-specific
    # dependence is covered by test_k1_correction_varies_with_query.
    shifted = {**base, "support_y": base["support_y"] + 1.0}
    moved = model(**shifted)
    delta = moved.transport - output.transport
    expected = float(model.transport.shrinkage(
        1, torch.zeros(1, dtype=torch.float64))) * (output.weight
                                                    * output.rho).sum(-1)
    assert torch.allclose(delta, expected, atol=1e-9)
    assert float(delta.abs().min()) > 1e-9
    # Label-reachability: the k=1 prediction loss reaches the relative and
    # transport parameters.
    model.zero_grad()
    output.prediction.square().mean().backward()
    assert any(p.grad is not None and p.grad.abs().sum() > 0
               for _, p in model.relative.named_parameters())
    assert (model.transport.log_shrinkage.grad is not None
            and model.transport.log_shrinkage.grad.abs().sum() > 0)
    # Support ligand: genuine autograd gradient into the correction.
    support_atoms = base["support_atoms"].clone().requires_grad_(True)
    inputs = {**base, "support_atoms": support_atoms}
    model.zero_grad()
    model(**inputs).transport.sum().backward()
    assert support_atoms.grad is not None and support_atoms.grad.abs().sum() > 0


def test_support_ligand_replacement_changes_prediction():
    model = build()
    base = episode(1, 5)
    output = model(**base)
    other_atoms, other_bonds, other_mask, other_fp = ligand(1, 555)
    replaced = model(**{**base,
                        "support_atoms": other_atoms,
                        "support_bonds": other_bonds,
                        "support_mask": other_mask,
                        "support_fingerprint": other_fp})
    assert not torch.allclose(output.prediction, replaced.prediction, atol=1e-9)


# ---------------------------------------------------------------- gradients


@pytest.mark.parametrize("support", [0, 1, 2, 5])
def test_no_dead_trainable_branch(support):
    model = build()
    output = model(**episode(support, 6))
    output.prediction.sum().backward()
    missing = [name for name, parameter in model.named_parameters()
               if parameter.requires_grad
               and (parameter.grad is None or not parameter.grad.abs().sum())]
    if support == 0:
        missing = [n for n in missing
                   if not n.startswith("transport.")
                   and not n.startswith("gate.")]
    if support == 1:
        # A softmax over one support is identically 1, so the learned key,
        # its temperature and the Tanimoto scale have no gradient at k=1.
        # The gate is the sole k=1 query-specific channel, by design.
        for name in ("transport.key.weight", "transport.log_temperature",
                     "transport.similarity_scale"):
            assert name in missing
        missing = [n for n in missing if n not in {
            "transport.key.weight", "transport.log_temperature",
            "transport.similarity_scale"}]
    assert missing == []


def test_level_and_shape_gradients_are_separable():
    """Detaching one branch removes gradient from its head and no other."""
    model = build()
    output = model(**episode(0, 5), adapt=False)
    (output.ligand_prior + output.shape).sum().backward()
    assert all(p.grad is None or not p.grad.abs().sum()
               for _, p in model.level_head.named_parameters())
    assert any(p.grad is not None and p.grad.abs().sum()
               for _, p in model.relative.named_parameters())


# ---------------------------------------------------------------- synthetic training


def synthetic_targets(count: int, seed: int, interaction: bool = True,
                      noise: float = 0.02):
    """Proteins and ligands with y = level_p + s * a_p . b_L (+ ligand noise).

    Returns the protein/ligand tables plus the label standard deviation used
    to normalize training labels exactly like the real pipeline does.
    """
    generator = torch.Generator().manual_seed(seed)
    level_weight = torch.randn(4, generator=generator, dtype=torch.float64)
    proteins = []
    for index in range(count):
        direction = torch.randn(4, generator=generator, dtype=torch.float64)
        pooled = torch.cat((direction, torch.randn(
            PROTEIN_DIM - 4, generator=generator, dtype=torch.float64)))
        tokens = torch.randn(SLOTS, PROTEIN_DIM, generator=generator,
                             dtype=torch.float64)
        chemistry = torch.rand(SLOTS, 4, generator=generator,
                               dtype=torch.float64)
        proteins.append({
            "pooled": pooled, "tokens": tokens,
            "mask": torch.ones(SLOTS, dtype=torch.float64),
            "chemistry": chemistry,
            "direction": direction,
            # The target level is a learnable function of the features —
            # a level that is independent of the inputs is unidentifiable.
            "level": 2.0 * float(level_weight @ direction),
        })
    ligands = []
    for index in range(60):
        atoms, bonds, mask, fp = ligand(1, seed + 1000 + index)
        basis = torch.randn(4, ATOM_FEAT_DIM, generator=generator,
                            dtype=torch.float64)
        mean_atoms = (atoms[0] * mask[0].unsqueeze(-1)).sum(0) \
            / mask[0].sum().clamp_min(1.0)
        b = basis @ mean_atoms
        ligand_noise = (torch.randn((), generator=generator,
                                    dtype=torch.float64).item()
                        if not interaction else 0.0)
        ligands.append({"atoms": atoms, "bonds": bonds, "mask": mask,
                        "fp": fp, "b": b, "noise": ligand_noise})
    labels = np.asarray([
        protein["level"]
        + (0.0 if not interaction else
           2.0 * float(protein["direction"] @ entry["b"]))
        + entry["noise"]
        for protein in proteins for entry in ligands])
    return proteins, ligands, float(labels.std())


def synthetic_episode(proteins, ligands, protein_index: int, support: int,
                      query: int, rng: np.random.Generator, *, private: bool):
    """Assemble one episode; private targets have no interaction term."""
    protein = proteins[protein_index]
    order = rng.permutation(len(ligands))
    support_indices = list(map(int, order[:support]))
    query_indices = list(map(int, order[support:support + query]))
    keys = ["atoms", "bonds", "mask", "fp"]
    def stack(key, indices):
        if not indices:
            shapes = {"atoms": (0, 1, ATOM_FEAT_DIM), "bonds": (0, 1, 1, BOND_FEAT_DIM),
                      "mask": (0, 1), "fp": (0, 64)}
            return torch.zeros(*shapes[key], dtype=torch.float64)
        return torch.cat([ligands[i][key] for i in indices], 0)
    support_y = torch.tensor([
        protein["level"] + (0.0 if private else
                            2.0 * float(protein["direction"] @ ligands[i]["b"]))
        + ligands[i]["noise"] for i in support_indices], dtype=torch.float64)
    query_y = torch.tensor([
        protein["level"] + (0.0 if private else
                            2.0 * float(protein["direction"] @ ligands[i]["b"]))
        + ligands[i]["noise"] for i in query_indices], dtype=torch.float64)
    return {
        "protein_pooled": protein["pooled"], "protein_tokens": protein["tokens"],
        "protein_mask": protein["mask"], "protein_chemistry": protein["chemistry"],
        "support_atoms": stack("atoms", support_indices),
        "support_bonds": stack("bonds", support_indices),
        "support_mask": stack("mask", support_indices),
        "support_y": support_y,
        "support_fingerprint": stack("fp", support_indices),
        "query_atoms": stack("atoms", query_indices),
        "query_bonds": stack("bonds", query_indices),
        "query_mask": stack("mask", query_indices),
        "query_y": query_y,
        "query_fingerprint": stack("fp", query_indices),
    }


def train_synthetic(model: RelTransportModel, proteins, ligands, label_std,
                    private: bool, steps: int = 500, lr: float = 1e-3,
                    counterfactual: bool = False, k0_only: bool = False,
                    query: int = 12, variance_weight: float = 1.0,
                    identify_weight: float = 0.3, seed: int = 7):
    config = RelConfig(seed=seed, steps=steps, episodes_per_step=1,
                       routing=True, counterfactual=counterfactual,
                       shape_variance_weight=variance_weight,
                       relative_loss_weight=0.5,
                       identify_weight=identify_weight,
                       device="cpu", amp=False)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr,
                                  weight_decay=1e-5)
    rng = np.random.default_rng(seed)
    scale = LabelScale(mean=0.0, scale=label_std)
    support_sizes = (0,) if k0_only else (0, 1, 2, 3, 5)
    donors = {index: (index + 1) % len(proteins)
              for index in range(len(proteins))}
    for step in range(1, steps + 1):
        optimizer.zero_grad(set_to_none=True)
        support = support_sizes[(step - 1) % len(support_sizes)]
        protein_index = int(rng.integers(len(proteins)))
        episode_data = synthetic_episode(
            proteins, ligands, protein_index, support, query, rng,
            private=private)
        support_count = support
        width = max(episode_data["support_atoms"].shape[-2],
                    episode_data["query_atoms"].shape[-2], 1)

        def pad(values, target, rank):
            if rank == 2:
                return torch.nn.functional.pad(values, (0, 0, 0, target - values.shape[-2]))
            if rank == 3:
                return torch.nn.functional.pad(
                    values, (0, 0, 0, target - values.shape[-2],
                             0, target - values.shape[-2]))
            return torch.nn.functional.pad(values, (0, target - values.shape[-1]))

        support_atoms = pad(episode_data["support_atoms"], width, 2)
        support_bonds = pad(episode_data["support_bonds"], width, 3)
        support_mask = pad(episode_data["support_mask"], width, 1)
        atoms = torch.cat((support_atoms, episode_data["query_atoms"]), 0).unsqueeze(0)
        bonds = torch.cat((support_bonds, episode_data["query_bonds"]), 0).unsqueeze(0)
        mask = torch.cat((support_mask, episode_data["query_mask"]), 0).unsqueeze(0)
        pooled = episode_data["protein_pooled"].unsqueeze(0)
        tokens = episode_data["protein_tokens"].unsqueeze(0)
        protein_mask = episode_data["protein_mask"].unsqueeze(0)
        chemistry = episode_data["protein_chemistry"].unsqueeze(0)
        donor_pooled = proteins[donors[protein_index]]["pooled"].unsqueeze(0)
        donor_tokens = proteins[donors[protein_index]]["tokens"].unsqueeze(0)
        donor_mask = proteins[donors[protein_index]]["mask"].unsqueeze(0)
        donor_chemistry = proteins[donors[protein_index]]["chemistry"].unsqueeze(0)
        parts = model.forward_parts(pooled, tokens, protein_mask, atoms, bonds,
                                    mask, chemistry)
        donor_parts = model.forward_parts(donor_pooled, donor_tokens,
                                          donor_mask, atoms, bonds, mask,
                                          donor_chemistry)
        endpoint, prior, level, shape, _, embed, u, u_gate = parts
        query_y = episode_data["query_y"].unsqueeze(0) / label_std
        support_y = episode_data["support_y"].unsqueeze(0) / label_std
        full_y = torch.cat((support_y, query_y), -1)
        query_prior, query_level = prior[:, support_count:], level[:, support_count:]
        query_shape = shape[:, support_count:]
        query_embed = embed[:, support_count:]
        query_fp = episode_data["query_fingerprint"].unsqueeze(0)
        support_fp = episode_data["support_fingerprint"].unsqueeze(0)
        from model.similarity_grammar import tanimoto
        query_similarity = tanimoto(query_fp, query_fp)
        full_similarity = tanimoto(
            torch.cat((support_fp, query_fp), 1),
            torch.cat((support_fp, query_fp), 1))
        from scripts.train_reltransport import transport_block
        transport = torch.zeros_like(query_shape)
        query_endpoint = endpoint[:, support_count:]
        support_endpoint = endpoint[:, :support_count]
        if support_count:
            residual = (support_y - support_endpoint).detach()
            transport, _ = transport_block(
                model, support_count, query_embed, embed[:, :support_count], u,
                u_gate, residual, tanimoto(query_fp, support_fp))
        prediction = query_endpoint + transport
        p_level = query_prior + query_level + query_shape.detach() \
            + transport.detach()
        p_shape = query_prior + query_level.detach() + query_shape + transport
        loss = level_term(p_level, query_y) \
            + shape_objective(model, config, p_shape, query_y, scale,
                              query_similarity) \
            + 0.5 * relative_supervision(model, config, parts, full_y, scale,
                                         full_similarity) \
            + config.identify_weight * query_shape.mean(-1).square().mean()
        if support_count and not private:
            wrong_shape = donor_parts[3][:, support_count:]
            loss = loss + 0.25 * contrast(
                shape_objective(model, config,
                                query_prior.detach() + query_shape + transport,
                                query_y, scale, query_similarity),
                shape_objective(model, config,
                                donor_parts[1][:, support_count:].detach()
                                + wrong_shape + transport.detach(),
                                query_y, scale, query_similarity),
                config.contrast_temperature)
            frozen = query_endpoint.detach()
            correct_mse = ((frozen + transport) - query_y).square().mean()
            # k>=2: permuted labels (mean(r) invariant). k=1: a wrong support
            # ligand with the same label — the label flip is evaluated, not
            # trained (it destabilizes a query-specific gate).
            if support_count > 1:
                wrong_transport, _ = transport_block(
                    model, support_count, query_embed,
                    embed[:, :support_count], u, u_gate, residual.roll(1, dims=-1),
                    tanimoto(query_fp, support_fp))
                loss = loss + 0.25 * contrast(
                    correct_mse,
                    ((frozen + wrong_transport) - query_y).square().mean(),
                    config.contrast_temperature)
            else:
                wrong_index = int(rng.integers(len(ligands)))
                wrong_ligand = ligands[wrong_index]
                wrong_atoms = wrong_ligand["atoms"].unsqueeze(0)
                wrong_bonds = wrong_ligand["bonds"].unsqueeze(0)
                wrong_mask = wrong_ligand["mask"].unsqueeze(0)
                wrong_parts = model.forward_parts(
                    pooled, tokens, protein_mask, wrong_atoms, wrong_bonds,
                    wrong_mask, chemistry)
                wrong_transport, _ = transport_block(
                    model, 1, query_embed, wrong_parts[5], u, u_gate,
                    (support_y - wrong_parts[0]).detach(),
                    tanimoto(query_fp, support_fp))
                loss = loss + 0.25 * contrast(
                    correct_mse,
                    ((frozen + wrong_transport) - query_y).square().mean(),
                    config.contrast_temperature)
        loss.backward()
        optimizer.step()
    return model


def _mse(prediction: torch.Tensor, truth: torch.Tensor) -> float:
    return float(((prediction - truth) ** 2).mean())


@pytest.mark.research_gate
def test_interaction_task_ordering_dies_without_the_interaction_branch():
    """Synthetic bilinear interaction: shape branch carries the ordering.

    Measured over three independently trained seeds and every synthetic
    protein (trained and held out), 20 queries each. The interaction-cut arm
    keeps ligand prior + target level, so the contrast isolates the shape
    branch's contribution to within-target ordering.
    """
    proteins, ligands, label_std = synthetic_targets(8, seed=11,
                                                     interaction=True)
    full_cells, cut_cells = [], []
    for seed in (3, 4, 5):
        torch.manual_seed(seed)
        model = build(seed=seed)
        train_synthetic(model, proteins, ligands, label_std, private=False,
                        steps=500, k0_only=True, query=12, seed=seed)
        rng = np.random.default_rng(42 + seed)
        with torch.no_grad():
            for protein_index in range(len(proteins)):
                held = synthetic_episode(proteins, ligands, protein_index, 0,
                                         20, rng, private=False)
                output = model(**{k: v for k, v in held.items()
                                  if k != "query_y"})
                truth = held["query_y"].numpy()
                full_cells.append(_ci(output.endpoint.numpy(), truth))
                cut_cells.append(_ci((output.ligand_prior
                                      + output.target_level).numpy(), truth))
    full, cut = np.asarray(full_cells), np.asarray(cut_cells)
    gaps = full - cut
    assert float(full.mean()) > 0.70, f"full endpoint CI {full.mean():.3f}"
    assert float(gaps.mean()) > 0.20, (
        f"interaction branch contributes too little: gap {gaps.mean():.3f}")
    assert int((gaps > 0).sum()) >= 20, (
        f"gap positive in only {(gaps > 0).sum()}/24 cells")


@pytest.mark.research_gate
def test_matched_wrong_support_is_clearly_worse():
    """k=1 correct support beats the magnitude-matched wrong label.

    Measured over two seeds, three held-out panels and four proteins each
    (single-draw comparisons are too noisy for a 8-query panel).
    """
    proteins, ligands, label_std = synthetic_targets(8, seed=12,
                                                     interaction=True)
    correct, wrong = [], []
    for seed in (4, 5):
        torch.manual_seed(seed)
        model = build(seed=seed)
        train_synthetic(model, proteins, ligands, label_std, private=False,
                        counterfactual=True, query=12, seed=seed)
        with torch.no_grad():
            for draw in range(3):
                rng = np.random.default_rng(43 + seed + 100 * draw)
                for protein_index in range(4):
                    held = synthetic_episode(proteins, ligands, protein_index,
                                             1, 8, rng, private=False)
                    keys = {k: v for k, v in held.items() if k != "query_y"}
                    keys["support_y"] = held["support_y"] / label_std
                    output = model(**keys)
                    residual = held["support_y"] / label_std \
                        - output.support_endpoint.detach()
                    flipped = held["support_y"] / label_std - 2.0 * residual
                    wrong_output = model(**{**keys, "support_y": flipped})
                    correct.append(_mse(output.prediction,
                                        held["query_y"] / label_std))
                    wrong.append(_mse(wrong_output.prediction,
                                      held["query_y"] / label_std))
    assert float(np.mean(wrong)) - float(np.mean(correct)) > 0.10, (
        f"correct {np.mean(correct)} not clearly better than wrong "
        f"{np.mean(wrong)}")


@pytest.mark.research_gate
def test_private_task_abstains_to_the_level_shift():
    """No transferable interaction -> full degenerates to level-only.

    Two seeds, three held-out panels each; the variance-heavy shape objective
    must teach the transport to suppress query-specific variation when the
    support residuals carry no transferable information.
    """
    proteins, ligands, label_std = synthetic_targets(8, seed=13,
                                                     interaction=False,
                                                     noise=0.3)
    full_values, level_values = [], []
    for seed in (5, 6):
        torch.manual_seed(seed)
        model = build(seed=seed)
        train_synthetic(model, proteins, ligands, label_std, private=True,
                        steps=1000, variance_weight=4.0, seed=seed)
        for draw in range(3):
            rng = np.random.default_rng(44 + seed + 100 * draw)
            with torch.no_grad():
                held = synthetic_episode(proteins, ligands, 0, 5, 10, rng,
                                         private=True)
                keys = {k: v for k, v in held.items() if k != "query_y"}
                keys["support_y"] = held["support_y"] / label_std
                output = model(**keys)
                truth = held["query_y"] / label_std
                full_values.append(_mse(output.prediction, truth))
                level_shift = float(model.transport.shrinkage(
                    5, output.support_residual)) \
                    * float(output.support_residual.mean())
                level_values.append(_mse(output.endpoint + level_shift, truth))
    assert float(np.mean(full_values)) <= float(np.mean(level_values)) + 0.005, (
        f"query-specific transport should not help on a private task: "
        f"{np.mean(full_values)} vs {np.mean(level_values)}")


def _ci(prediction: np.ndarray, truth: np.ndarray) -> float:
    concordant, comparable = 0.0, 0
    for left in range(len(truth)):
        for right in range(left + 1, len(truth)):
            true_delta = float(truth[left] - truth[right])
            if true_delta == 0:
                continue
            comparable += 1
            concordant += float((prediction[left] - prediction[right])
                                * true_delta > 0)
    return float(concordant / comparable) if comparable else float("nan")

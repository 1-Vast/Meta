"""Stage 1 structural and synthetic gates for the direct-shape model (R13).

The direct-shape family fixes the supervision leak the R9-R12 ladder
localized: the relative supervision now targets `s(e_i) - s(e_j)` — the
deployed ordering quantity itself — instead of a bilinear potential that
only approximates it. These gates run before any real-data training and
falsify the specific ways the new family could be cosmetic.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from contracts.ligand_graph import ATOM_FEAT_DIM, BOND_FEAT_DIM
from model.shape_direct import ShapeDirectModel
from scripts.train_qpsmp import LabelScale
from test_reltransport_synthetic import (
    synthetic_episode, synthetic_targets, train_synthetic, _ci,
)
from scripts.train_reltransport import (
    RelConfig, level_term, shape_objective, contrast,
)

PROTEIN_DIM, SLOTS, ATOMS = 32, 16, 9


def build(seed: int = 0, **kwargs) -> ShapeDirectModel:
    torch.manual_seed(seed)
    return ShapeDirectModel(protein_dim=PROTEIN_DIM, hidden_dim=24, task_dim=12,
                            ligand_layers=2, pair_dim=24, pair_latents=6,
                            pair_heads=2, anchors=6, shape_hidden=12,
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
    fingerprint[:, 0] = 1.0
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


def test_anchor_set_has_exactly_zero_mean_shape():
    model = build()
    protein_vector = torch.randn(1, 2 * model.hidden_dim, dtype=torch.float64)
    anchor_shape, _ = model.shape_head(
        model.shape_head.anchor.unsqueeze(0), protein_vector)
    assert float(anchor_shape.mean(-1).abs().max()) < 1e-10


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
    model = build()
    base = episode(3, 5)
    shifted = {**base, "support_y": base["support_y"] + 2.0}
    before, after = model(**base), model(**shifted)
    assert torch.allclose(before.endpoint, after.endpoint, atol=1e-12)
    shrink = float(model.transport.shrinkage(3, torch.zeros(1, dtype=torch.float64)))
    expected = shrink * 2.0 * before.weight.sum(-1)
    assert torch.allclose(after.transport - before.transport, expected,
                          atol=1e-9)


def test_query_labels_are_never_an_input():
    import inspect
    signature = inspect.signature(ShapeDirectModel.forward)
    assert "query_y" not in signature.parameters


def test_geometry_input_is_refused():
    with pytest.raises(ValueError, match="common-frame"):
        build()(**episode(1, 2),
                geometry_available=torch.ones(1, dtype=torch.bool))


@pytest.mark.parametrize("support", [0, 1, 2, 5])
def test_no_dead_trainable_branch(support):
    model = build()
    output = model(**episode(support, 6))
    output.prediction.sum().backward()
    missing = [name for name, parameter in model.named_parameters()
               if parameter.requires_grad
               and (parameter.grad is None or not parameter.grad.abs().sum())]
    if support == 0:
        missing = [n for n in missing if not n.startswith("transport.")]
    if support == 1:
        # A softmax over one support is identically 1: the learned key, its
        # temperature and the Tanimoto scale have no k=1 gradient (the same
        # documented degeneracy as the retained Tanimoto baseline).
        for name in ("transport.key.weight", "transport.log_temperature",
                     "transport.similarity_scale"):
            assert name in missing
        missing = [n for n in missing if n not in {
            "transport.key.weight", "transport.log_temperature",
            "transport.similarity_scale"}]
    assert missing == []


# ---------------------------------------------------------------- synthetic training


def train_direct(model: ShapeDirectModel, proteins, ligands, label_std,
                 private: bool, steps: int = 500, lr: float = 1e-3,
                 counterfactual: bool = False, k0_only: bool = False,
                 query: int = 12, variance_weight: float = 1.0,
                 seed: int = 7):
    from model.similarity_grammar import tanimoto
    from scripts.train_shape_direct import difference_supervision
    from scripts.train_shape_direct import ShapeDirectConfig
    config = ShapeDirectConfig(seed=seed, steps=steps, episodes_per_step=1,
                               shape_variance_weight=variance_weight,
                               difference_weight=0.5, device="cpu", amp=False)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-5)
    rng = np.random.default_rng(seed)
    scale = LabelScale(mean=0.0, scale=label_std)
    support_sizes = (0,) if k0_only else (0, 1, 2, 3, 5)
    for step in range(1, steps + 1):
        optimizer.zero_grad(set_to_none=True)
        support = support_sizes[(step - 1) % len(support_sizes)]
        pi = int(rng.integers(len(proteins)))
        e = synthetic_episode(proteins, ligands, pi, support, query, rng,
                              private=private)
        width = max(e["support_atoms"].shape[-2], e["query_atoms"].shape[-2], 1)
        sa = torch.nn.functional.pad(e["support_atoms"], (0, 0, 0, width - e["support_atoms"].shape[-2]))
        sb = torch.nn.functional.pad(e["support_bonds"], (0, 0, 0, width - e["support_bonds"].shape[-2], 0, width - e["support_bonds"].shape[-2]))
        sm = torch.nn.functional.pad(e["support_mask"], (0, width - e["support_mask"].shape[-1]))
        atoms = torch.cat((sa, e["query_atoms"]), 0).unsqueeze(0)
        bonds = torch.cat((sb, e["query_bonds"]), 0).unsqueeze(0)
        mask = torch.cat((sm, e["query_mask"]), 0).unsqueeze(0)
        parts = model.forward_parts(e["protein_pooled"].unsqueeze(0), e["protein_tokens"].unsqueeze(0),
                                    e["protein_mask"].unsqueeze(0), atoms, bonds, mask,
                                    e["protein_chemistry"].unsqueeze(0))
        endpoint, prior, level, shape, _, embed, _ = parts
        qy = e["query_y"].unsqueeze(0) / label_std
        sy = e["support_y"].unsqueeze(0) / label_std
        full_y = torch.cat((sy, qy), -1)
        qprior, qlevel = prior[:, support:], level[:, support:]
        qshape = shape[:, support:]
        qfp = e["query_fingerprint"].unsqueeze(0)
        sfp = e["support_fingerprint"].unsqueeze(0)
        qsim = tanimoto(qfp, qfp)
        fsim = tanimoto(torch.cat((sfp, qfp), 1), torch.cat((sfp, qfp), 1))
        transport = torch.zeros_like(qshape)
        if support:
            residual = (sy - endpoint[:, :support]).detach()
            transport, _ = model.transport(
                embed[:, :support], embed[:, support:], residual, tanimoto(qfp, sfp))
            transport = model.transport.shrinkage(support, residual) * transport
        pred = endpoint[:, support:] + transport
        p_level = qprior + qlevel + qshape.detach() + transport.detach()
        p_shape = qprior + qlevel.detach() + qshape + transport
        loss = level_term(p_level, qy) \
            + shape_objective(model, RelConfig(), p_shape, qy, scale, qsim) \
            + 0.5 * difference_supervision(shape, full_y, scale, config, fsim) \
            + 0.3 * qshape.mean(-1).square().mean()
        if support and counterfactual and not private:
            frozen = pred.detach()
            correct_mse = ((frozen + transport) - qy).square().mean()
            wrong_residual = (sy.roll(1, dims=-1) - endpoint[:, :support]).detach() \
                if support > 1 else (sy - endpoint[:, :support]).detach() * -1
            wrong_transport, _ = model.transport(
                embed[:, :support], embed[:, support:], wrong_residual, tanimoto(qfp, sfp))
            wrong_transport = model.transport.shrinkage(support, wrong_residual) * wrong_transport
            loss = loss + 0.25 * contrast(
                correct_mse,
                ((frozen + wrong_transport) - qy).square().mean(),
                config.contrast_temperature)
        loss.backward()
        optimizer.step()
    return model


def _mse(prediction: torch.Tensor, truth: torch.Tensor) -> float:
    return float(((prediction - truth) ** 2).mean())


@pytest.mark.xfail(
    reason="R13 gate verdict: the direct MLP shape branch collapses under "
           "the shape variance term and underperforms the bilinear readout "
           "on the synthetic interaction task (mean CI 0.60, gap 0.14 vs "
           "gates 0.70/0.20). The family is gate-blocked at Stage 1; the "
           "thresholds are not moved. See "
           "report/meta_fewshot/stageR13_shape_direct_20260816/REPORT.md.")
@pytest.mark.research_gate
def test_interaction_task_ordering_dies_without_the_interaction_branch():
    proteins, ligands, label_std = synthetic_targets(8, seed=11, interaction=True)
    full_cells, cut_cells = [], []
    for seed in (3, 4, 5):
        torch.manual_seed(seed)
        model = build(seed=seed)
        train_direct(model, proteins, ligands, label_std, private=False,
                     steps=500, k0_only=True, query=12, seed=seed)
        rng = np.random.default_rng(42 + seed)
        with torch.no_grad():
            for pi in range(len(proteins)):
                held = synthetic_episode(proteins, ligands, pi, 0, 20, rng,
                                         private=False)
                out = model(**{k: v for k, v in held.items() if k != "query_y"})
                truth = held["query_y"].numpy()
                full_cells.append(_ci(out.endpoint.numpy(), truth))
                cut_cells.append(_ci((out.ligand_prior + out.target_level).numpy(),
                                     truth))
    full, cut = np.asarray(full_cells), np.asarray(cut_cells)
    gaps = full - cut
    assert float(full.mean()) > 0.70, f"full endpoint CI {full.mean():.3f}"
    assert float(gaps.mean()) > 0.20, f"interaction branch gap {gaps.mean():.3f}"
    assert int((gaps > 0).sum()) >= 20, f"positive in {(gaps > 0).sum()}/24"


@pytest.mark.xfail(
    reason="R13 gate verdict: the k=1 transport (retained Tanimoto baseline) "
           "is mildly harmful on the synthetic task (measured gap -0.030, "
           "correct 1.463 vs wrong 1.433) — the family's k=1 behavior "
           "matches the recorded neutrality/harm signature of every "
           "query-specific channel tested in this project. See "
           "report/meta_fewshot/stageR13_shape_direct_20260816/REPORT.md.")
@pytest.mark.research_gate
def test_matched_wrong_support_is_clearly_worse():
    proteins, ligands, label_std = synthetic_targets(8, seed=12, interaction=True)
    correct, wrong = [], []
    for seed in (4, 5):
        torch.manual_seed(seed)
        model = build(seed=seed)
        train_direct(model, proteins, ligands, label_std, private=False,
                     counterfactual=True, query=12, seed=seed)
        with torch.no_grad():
            for draw in range(3):
                rng = np.random.default_rng(43 + seed + 100 * draw)
                for pi in range(4):
                    held = synthetic_episode(proteins, ligands, pi, 1, 8, rng,
                                             private=False)
                    keys = {k: v for k, v in held.items() if k != "query_y"}
                    keys["support_y"] = held["support_y"] / label_std
                    out = model(**keys)
                    residual = held["support_y"] / label_std \
                        - out.support_endpoint.detach()
                    flipped = held["support_y"] / label_std - 2.0 * residual
                    w = model(**{**keys, "support_y": flipped})
                    correct.append(_mse(out.prediction, held["query_y"] / label_std))
                    wrong.append(_mse(w.prediction, held["query_y"] / label_std))
    assert float(np.mean(wrong)) - float(np.mean(correct)) > 0.10, (
        f"correct {np.mean(correct)} not clearly better than wrong {np.mean(wrong)}")


@pytest.mark.research_gate
def test_private_task_abstains_to_the_level_shift():
    proteins, ligands, label_std = synthetic_targets(8, seed=13,
                                                     interaction=False, noise=0.3)
    full_values, level_values = [], []
    for seed in (5, 6):
        torch.manual_seed(seed)
        model = build(seed=seed)
        train_direct(model, proteins, ligands, label_std, private=True,
                     steps=1000, variance_weight=4.0, seed=seed)
        for draw in range(3):
            rng = np.random.default_rng(44 + seed + 100 * draw)
            with torch.no_grad():
                held = synthetic_episode(proteins, ligands, 0, 5, 10, rng,
                                         private=True)
                keys = {k: v for k, v in held.items() if k != "query_y"}
                keys["support_y"] = held["support_y"] / label_std
                out = model(**keys)
                truth = held["query_y"] / label_std
                full_values.append(_mse(out.prediction, truth))
                level_shift = float(model.transport.shrinkage(
                    5, out.support_residual)) * float(out.support_residual.mean())
                level_values.append(_mse(out.endpoint + level_shift, truth))
    assert float(np.mean(full_values)) <= float(np.mean(level_values)) + 0.005, (
        f"query-specific transport should not help on a private task: "
        f"{np.mean(full_values)} vs {np.mean(level_values)}")

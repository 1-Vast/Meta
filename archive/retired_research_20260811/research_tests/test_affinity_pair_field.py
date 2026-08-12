import itertools

import numpy as np
import pytest
import torch

from research.meta_fewshot.affinity_pair_field import (
    AffinityContrastLoss,
    AffinityDirectedPairField,
    PairDifference,
    RectangleDifference,
    cluster_partner_necessity,
    coarse_interaction_compatibility,
    expand_slot_geometry_prior,
    exact_distance_loss,
    exact_distance_loss_per_system,
    rectangle_values,
)


def example_batch(dtype=torch.float64):
    torch.manual_seed(12)
    atom_states = torch.randn(2, 3, 4, dtype=dtype)
    residue_states = torch.randn(2, 4, 5, dtype=dtype)
    atom_mask = torch.tensor([[1, 1, 0], [1, 1, 1]], dtype=torch.bool)
    residue_mask = torch.tensor([[1, 1, 1, 0], [1, 1, 1, 1]], dtype=torch.bool)
    atom_indices = torch.tensor([[7, 11, -1], [2, 5, 9]])
    residue_indices = torch.tensor([[19, 23, 31, -1], [3, 8, 14, 21]])
    atom_classes = torch.tensor([[2, 4, -1], [0, 1, 7]])
    residue_classes = torch.tensor([[4, 2, 3, -1], [0, 1, 2, 5]])
    compatibility = coarse_interaction_compatibility(atom_classes, residue_classes).to(dtype)
    raw = torch.rand(2, 3, 4, 5, dtype=dtype)
    distance_prior = raw / raw.sum(-1, keepdim=True)
    contact_prior = distance_prior[..., :2].sum(-1)
    return dict(
        atom_states=atom_states,
        atom_mask=atom_mask,
        residue_states=residue_states,
        residue_mask=residue_mask,
        contact_prior=contact_prior,
        distance_prior=distance_prior,
        compatibility=compatibility,
        atom_indices=atom_indices,
        residue_indices=residue_indices,
    )


def test_exact_pair_field_is_traceable_masked_and_differentiable():
    field = AffinityDirectedPairField(4, 5, rank=3, dtype=torch.float64)
    batch = example_batch()
    prediction = field(**batch)
    assert prediction.typed_pair_energy.shape == (2, 3, 4, 6)
    assert prediction.typed_summary.shape == (2, 6)
    assert prediction.score.shape == (2,)
    assert prediction.distance_prob.shape == (2, 3, 4, 5)
    assert prediction.atom_indices is batch["atom_indices"]
    assert prediction.residue_indices is batch["residue_indices"]
    assert prediction.typed_pair_energy[0, 2].count_nonzero() == 0
    assert prediction.typed_pair_energy[0, :, 3].count_nonzero() == 0
    loss = prediction.score.square().mean() + exact_distance_loss(
        prediction, torch.zeros(2, 3, 4, dtype=torch.long))
    loss.backward()
    assert field.atom_projection.weight.grad is not None
    assert field.residue_projection.weight.grad is not None
    assert field.distance_residual.weight.grad is not None


def test_additive_null_has_matched_parameters_without_pair_product():
    full = AffinityDirectedPairField(4, 5, rank=3)
    additive = AffinityDirectedPairField(
        4, 5, rank=3, interaction_mode="additive")
    assert sum(value.numel() for value in full.parameters()) == sum(
        value.numel() for value in additive.parameters())
    batch = example_batch(dtype=torch.float32)
    baseline = additive(**batch).distance_prob
    changed = dict(batch)
    changed["compatibility"] = 1.0 - batch["compatibility"]
    assert torch.equal(baseline, additive(**changed).distance_prob)


def test_distance_loss_is_meaned_within_each_system():
    field = AffinityDirectedPairField(4, 5, rank=3, dtype=torch.float64)
    batch = example_batch()
    prediction = field(**batch)
    labels = torch.zeros(2, 3, 4, dtype=torch.long)
    together = exact_distance_loss_per_system(prediction, labels)
    separate = []
    for index in range(2):
        one = {name: value[index:index + 1] for name, value in batch.items()}
        separate.append(exact_distance_loss_per_system(
            field(**one), labels[index:index + 1])[0])
    assert torch.allclose(together, torch.stack(separate), atol=1e-12)


def test_pair_score_is_invariant_to_joint_exact_residue_permutation():
    field = AffinityDirectedPairField(4, 5, rank=3, dtype=torch.float64)
    batch = example_batch()
    baseline = field(**batch)
    permutation = torch.tensor([2, 0, 1, 3])
    changed = dict(batch)
    for name in ("residue_states", "residue_mask", "residue_indices"):
        changed[name] = batch[name][:, permutation]
    changed["contact_prior"] = batch["contact_prior"][:, :, permutation]
    changed["distance_prior"] = batch["distance_prior"][:, :, permutation]
    changed["compatibility"] = batch["compatibility"][:, :, permutation]
    permuted = field(**changed)
    assert torch.allclose(baseline.score, permuted.score, atol=1e-12)
    assert torch.allclose(
        baseline.typed_pair_energy[:, :, permutation],
        permuted.typed_pair_energy,
        atol=1e-12,
    )


def test_exact_identity_contract_rejects_duplicate_or_untracked_residues():
    field = AffinityDirectedPairField(4, 5, dtype=torch.float64)
    batch = example_batch()
    batch["residue_indices"] = batch["residue_indices"].clone()
    batch["residue_indices"][0, 1] = batch["residue_indices"][0, 0]
    with pytest.raises(ValueError, match="unique"):
        field(**batch)
    batch = example_batch()
    batch["atom_indices"] = batch["atom_indices"].clone()
    batch["atom_indices"][0, 2] = 99
    with pytest.raises(ValueError, match="masked"):
        field(**batch)


def test_distance_posterior_starts_from_frozen_prior_not_binary_contact():
    field = AffinityDirectedPairField(4, 5, dtype=torch.float64)
    with torch.no_grad():
        field.distance_residual.weight.zero_()
    batch = example_batch()
    prediction = field(**batch)
    active = prediction.pair_mask
    assert torch.allclose(
        prediction.distance_prob[active], batch["distance_prior"][active], atol=1e-12)
    assert torch.allclose(
        prediction.contact_prob[active],
        batch["distance_prior"][..., :2].sum(-1)[active],
        atol=1e-12,
    )


def test_slot_prior_is_only_a_prior_and_exact_residue_rows_stay_distinct():
    contact_slots = torch.tensor([[[0.2, 0.8], [0.4, 0.6]]], dtype=torch.float64)
    distance_slots = torch.tensor(
        [[[[.1, .2, .3, .2, .2], [.4, .2, .1, .1, .2]],
          [[.2, .2, .2, .2, .2], [.1, .1, .2, .3, .3]]]],
        dtype=torch.float64,
    )
    residue_slot = torch.tensor([[1, 1, 0]])
    residue_mask = torch.ones(1, 3, dtype=torch.bool)
    contact, distance = expand_slot_geometry_prior(
        contact_slots, distance_slots, residue_slot, residue_mask)
    assert torch.equal(contact[:, :, 0], contact[:, :, 1])
    assert torch.equal(distance[:, :, 0], distance[:, :, 1])

    field = AffinityDirectedPairField(2, 2, rank=2, dtype=torch.float64)
    atom_states = torch.tensor([[[1., 0.], [0., 1.]]], dtype=torch.float64)
    residue_states = torch.tensor(
        [[[1., 0.], [0., 1.], [1., 1.]]], dtype=torch.float64)
    compatibility = torch.ones(1, 2, 3, 6, dtype=torch.float64)
    prediction = field(
        atom_states, torch.ones(1, 2, dtype=torch.bool),
        residue_states, residue_mask, contact, distance, compatibility,
        torch.tensor([[4, 9]]), torch.tensor([[10, 11, 12]]),
    )
    changed_states = residue_states.clone()
    changed_states[0, 1] += 2.0
    changed = field(
        atom_states, torch.ones(1, 2, dtype=torch.bool),
        changed_states, residue_mask, contact, distance, compatibility,
        torch.tensor([[4, 9]]), torch.tensor([[10, 11, 12]]),
    )
    assert torch.equal(
        prediction.typed_pair_energy[:, :, 0], changed.typed_pair_energy[:, :, 0])
    assert not torch.equal(
        prediction.typed_pair_energy[:, :, 1], changed.typed_pair_energy[:, :, 1])


def test_measured_rectangle_removes_additive_protein_and_ligand_effects():
    protein = torch.tensor([3.0, 3.0, -2.0, -2.0])
    ligand = torch.tensor([5.0, 1.0, 5.0, 1.0])
    interaction = torch.tensor([1.5, -0.5, -2.0, 0.25])
    scores = protein + ligand + interaction
    rectangles = RectangleDifference(
        indices=torch.tensor([[0, 1, 2, 3]]),
        target=torch.tensor([0.0]),
    )
    observed = rectangle_values(scores, rectangles)
    expected = interaction[0] - interaction[1] - interaction[2] + interaction[3]
    assert observed.item() == pytest.approx(expected.item())


def test_contrast_loss_uses_only_explicit_measured_cells_and_trains():
    scores = torch.tensor([0.2, 0.8, -0.4, 0.1], requires_grad=True)
    within = PairDifference(
        left=torch.tensor([0, 2]), right=torch.tensor([1, 3]),
        target=torch.tensor([-0.7, -0.3]),
    )
    cross = PairDifference(
        left=torch.tensor([0]), right=torch.tensor([2]), target=torch.tensor([0.5]))
    rectangles = RectangleDifference(
        indices=torch.tensor([[0, 1, 2, 3]]), target=torch.tensor([-0.2]))
    output = AffinityContrastLoss()(
        scores, within_target=within, cross_protein=cross, rectangles=rectangles)
    assert set(output) == {"within_target", "cross_protein", "rectangle", "total"}
    output["total"].backward()
    assert scores.grad is not None
    assert torch.linalg.vector_norm(scores.grad) > 0

    reversed_within = PairDifference(
        left=within.right, right=within.left, target=-within.target)
    forward = AffinityContrastLoss(
        cross_protein_weight=0.0, rectangle_weight=0.0)(
            scores.detach(), within_target=within)["total"]
    reverse = AffinityContrastLoss(
        cross_protein_weight=0.0, rectangle_weight=0.0)(
            scores.detach(), within_target=reversed_within)["total"]
    assert torch.equal(forward, reverse)


def test_partner_identity_gate_is_cluster_level_and_directional():
    truth = torch.tensor([0., 1., 2., 3., 4., 5., 6., 7.])
    correct = truth + torch.tensor([.1, -.1, .2, -.2, .1, -.1, .2, -.2])
    wrong = truth + torch.tensor([1.2, -1.1, 1.0, -1.3, 1.1, -1.2, 1.3, -1.0])
    clusters = ["a", "a", "b", "b", "c", "c", "d", "d"]
    gate = cluster_partner_necessity(
        correct, wrong, truth, clusters, n_bootstrap=2000, seed=9)
    assert gate["partner_necessity"] > 0
    assert gate["lcb95_one_sided"] > 0
    reversed_gate = cluster_partner_necessity(
        wrong, correct, truth, clusters, n_bootstrap=2000, seed=9)
    assert reversed_gate["partner_necessity"] < 0
    assert reversed_gate["lcb95_one_sided"] < 0


def test_small_pair_field_learns_measured_interaction_rectangles():
    torch.manual_seed(5)
    protein_value = torch.tensor([-1.0, 0.5, 1.5])
    ligand_value = torch.tensor([-1.2, 0.3, 2.0])
    protein_states, atom_states, cell = [], [], {}
    for protein in range(3):
        for ligand in range(3):
            cell[protein, ligand] = len(protein_states)
            protein_states.append(protein_value[protein])
            atom_states.append(ligand_value[ligand])
    protein_states = torch.stack(protein_states).reshape(9, 1, 1)
    atom_states = torch.stack(atom_states).reshape(9, 1, 1)
    rectangles, targets = [], []
    for first_protein, second_protein in itertools.combinations(range(3), 2):
        for first_ligand, second_ligand in itertools.combinations(range(3), 2):
            rectangles.append([
                cell[first_protein, first_ligand],
                cell[first_protein, second_ligand],
                cell[second_protein, first_ligand],
                cell[second_protein, second_ligand],
            ])
            targets.append(
                protein_value[first_protein] * ligand_value[first_ligand]
                - protein_value[first_protein] * ligand_value[second_ligand]
                - protein_value[second_protein] * ligand_value[first_ligand]
                + protein_value[second_protein] * ligand_value[second_ligand]
            )
    measured = RectangleDifference(torch.tensor(rectangles), torch.stack(targets))
    mask = torch.ones(9, 1, dtype=torch.bool)
    distance_prior = torch.full((9, 1, 1, 5), 0.2)
    field = AffinityDirectedPairField(1, 1, rank=2)
    objective = AffinityContrastLoss(
        within_target_weight=0.0, cross_protein_weight=0.0,
        rectangle_weight=1.0,
    )
    optimizer = torch.optim.Adam(field.parameters(), lr=0.03)
    for _ in range(600):
        optimizer.zero_grad()
        scores = field(
            atom_states, mask, protein_states, mask,
            torch.full((9, 1, 1), 0.4), distance_prior,
            torch.ones(9, 1, 1, 6), torch.zeros(9, 1, dtype=torch.long),
            torch.zeros(9, 1, dtype=torch.long),
        ).score
        loss = objective(scores, rectangles=measured)["total"]
        loss.backward()
        optimizer.step()
    predicted = rectangle_values(scores.detach(), measured)
    assert torch.mean(torch.square(predicted - measured.target)) < 1e-4

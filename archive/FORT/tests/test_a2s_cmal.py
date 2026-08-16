from __future__ import annotations

import math

import pandas as pd
import torch

from research.a2s.a2s_cmal import (
    ProteinConditionedMetaAdapter,
    component_disjoint_source_partition,
    counterfactual_meta_loss,
    label_swap_support_y,
    source_mechanism_gate,
    verify_episode_measurement_identity,
)


def _example() -> tuple[ProteinConditionedMetaAdapter, dict[str, torch.Tensor]]:
    torch.manual_seed(17)
    model = ProteinConditionedMetaAdapter(
        ligand_dim=6,
        protein_dim=8,
        hidden=16,
        heads=4,
        dropout=0.0,
    ).eval()
    tensors = {
        "protein": torch.randn(2, 8),
        "support_ligand": torch.randn(2, 3, 6),
        "support_y": torch.randn(2, 3),
        "support_f0": torch.randn(2, 3),
        "query_ligand": torch.randn(2, 5, 6),
        "query_f0": torch.randn(2, 5),
    }
    return model, tensors


def test_support_free_path_is_exactly_independent_of_support() -> None:
    model, tensors = _example()
    first = model(**tensors, use_support=False)["prediction"]
    altered = dict(tensors)
    altered["support_ligand"] = tensors["support_ligand"] + 100.0
    altered["support_y"] = tensors["support_y"] - 50.0
    second = model(**altered, use_support=False)["prediction"]
    assert first.shape == (2, 5)
    torch.testing.assert_close(first, second)


def test_learned_adapter_is_support_sensitive_and_query_dependent() -> None:
    model, tensors = _example()
    with torch.no_grad():
        model.delta_head[-1].weight.fill_(0.1)
    first = model(**tensors)["delta"]
    altered = dict(tensors)
    altered["support_y"] = tensors["support_y"] + 2.0
    second = model(**altered)["delta"]
    assert not torch.allclose(first, second)
    assert first.var(dim=1).min() > 0


def test_adapter_cannot_bypass_measurements_with_a_query_only_delta() -> None:
    model, tensors = _example()
    with torch.no_grad():
        model.delta_head[-1].weight.fill_(0.1)
        support_base, _, _ = model.base_predict(
            tensors["protein"], tensors["support_ligand"], tensors["support_f0"]
        )
    neutral = dict(tensors)
    neutral["support_y"] = support_base
    prediction = model(**neutral)
    torch.testing.assert_close(prediction["delta"], torch.zeros_like(prediction["delta"]))
    torch.testing.assert_close(prediction["prediction"], prediction["base_prediction"])


def test_adapter_is_conditioned_on_recipient_protein() -> None:
    model, tensors = _example()
    with torch.no_grad():
        model.delta_head[-1].weight.fill_(0.1)
    first = model(**tensors)["prediction"]
    altered = dict(tensors)
    altered["protein"] = torch.flip(tensors["protein"], dims=(0,))
    second = model(**altered)["prediction"]
    assert not torch.allclose(first, second)


def test_label_swap_moves_residual_without_changing_correct_support_base() -> None:
    model, tensors = _example()
    donor_ligand = torch.randn(2, 3, 3, 6)
    donor_f0 = torch.randn(2, 3, 3)
    donor_y = torch.randn(2, 3, 3)
    batch = {
        **tensors,
        "negative_support_ligand": donor_ligand,
        "negative_support_f0": donor_f0,
        "negative_support_y": donor_y,
    }
    swapped = label_swap_support_y(model, batch)
    with torch.no_grad():
        correct_base, _, _ = model.base_predict(
            tensors["protein"], tensors["support_ligand"], tensors["support_f0"]
        )
        donor_base, _, _ = model.base_predict(
            tensors["protein"], donor_ligand[:, 0], donor_f0[:, 0]
        )
    torch.testing.assert_close(
        swapped - correct_base,
        donor_y[:, 0] - donor_base,
    )


def test_freeze_base_leaves_only_meta_adapter_trainable() -> None:
    model, _ = _example()
    model.freeze_base()
    assert all(not parameter.requires_grad for parameter in model.base_parameters())
    assert all(parameter.requires_grad for parameter in model.adapter_parameters())
    model.train()
    assert not model.ligand_encoder.training
    assert model.support_attention.training


def test_counterfactual_objective_scores_post_adaptation_query_ranking() -> None:
    target = torch.tensor([[3.0, 2.0, 1.0]])
    mask = torch.ones_like(target)
    base = torch.zeros_like(target)
    positive = target.clone()
    inverse = torch.flip(target, dims=(1,))
    negatives = inverse.unsqueeze(1).repeat(1, 3, 1)
    loss, detail = counterfactual_meta_loss(
        positive,
        negatives,
        base,
        target,
        mask,
        temperature=0.2,
    )
    assert loss < math.log(4.0)
    assert detail["positive_ranking_gain"].item() > 0
    assert detail["ranking_gap"].item() > 0


def test_counterfactual_objective_does_not_reward_further_wrong_support_harm() -> None:
    target = torch.tensor([[3.0, 2.0, 1.0]])
    mask = torch.ones_like(target)
    base = torch.zeros_like(target)
    positive = target.clone()
    inverse = torch.flip(target, dims=(1,))
    mildly_wrong = inverse.unsqueeze(1).repeat(1, 3, 1)
    severely_wrong = (100.0 * inverse).unsqueeze(1).repeat(1, 3, 1)
    mild_loss, mild_detail = counterfactual_meta_loss(
        positive, mildly_wrong, base, target, mask, temperature=0.2
    )
    severe_loss, severe_detail = counterfactual_meta_loss(
        positive, severely_wrong, base, target, mask, temperature=0.2
    )
    torch.testing.assert_close(mild_loss, severe_loss)
    torch.testing.assert_close(
        mild_detail["negative_ranking_gain"],
        torch.zeros_like(mild_detail["negative_ranking_gain"]),
    )
    torch.testing.assert_close(
        severe_detail["negative_ranking_gain"],
        torch.zeros_like(severe_detail["negative_ranking_gain"]),
    )


def test_episode_measurements_are_checked_against_target_and_parent() -> None:
    episodes = pd.DataFrame({
        "episode_id": ["e0"],
        "target_uid": ["t0"],
        "support_parent_uids": ['["p0"]'],
        "support_measurement_uids": ['["m0"]'],
        "query_parent_uids": ['["p1"]'],
        "query_measurement_uids": ['["m1"]'],
    })
    observations = pd.DataFrame({
        "target_uid": ["t0", "t0"],
        "compound_parent_uid": ["p0", "p1"],
        "measurement_uid": ["m0", "m1"],
    })
    assert verify_episode_measurement_identity(episodes, observations) == {
        "unique_measurements": 2,
        "episode_occurrences": 2,
    }
    observations.loc[1, "compound_parent_uid"] = "wrong"
    try:
        verify_episode_measurement_identity(episodes, observations)
    except ValueError as error:
        assert "measurement identity mismatch" in str(error)
    else:
        raise AssertionError("mismatched measurement identity was accepted")


def test_source_mechanism_gate_requires_absolute_and_wrong_support_gains() -> None:
    def snapshot(value: float) -> dict:
        metrics = {name: value for name in ("ci", "spearman", "ndcg10")}
        metrics["rmse"] = 1.0
        return {
            "absolute": {"adapted": metrics},
            "correct_support_advantage": {
                arm: {
                    name: {"mean": value}
                    for name in ("ci", "spearman", "ndcg10")
                }
                for arm in (
                    "adapted_random",
                    "adapted_protein_hard",
                    "adapted_chemical_match",
                    "adapted_label_swap",
                )
            },
        }

    assert source_mechanism_gate(snapshot(0.0), snapshot(0.1))["status"] == "PASS"
    failed = snapshot(0.1)
    failed["correct_support_advantage"]["adapted_chemical_match"]["ci"]["mean"] = -0.1
    assert source_mechanism_gate(snapshot(0.0), failed)["status"] == "FAIL"


def test_source_training_partition_keeps_homology_components_intact() -> None:
    splits = pd.DataFrame({
        "target_uid": ["t0", "t1", "t2", "t3", "t4"],
        "component_id": [0, 0, 1, 2, 3],
        "role": ["source"] * 5,
        "meta_split": ["meta_train"] * 5,
    })
    base, adapter, audit = component_disjoint_source_partition(splits, seed=17)
    assert base | adapter == set(splits.target_uid)
    assert not (base & adapter)
    assert ({"t0", "t1"} <= base) or ({"t0", "t1"} <= adapter)
    assert audit["target_overlap"] == 0
    assert audit["component_overlap"] == 0

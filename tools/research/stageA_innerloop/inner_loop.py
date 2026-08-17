"""Functional first-order inner/outer-loop adaptation on the 97-parameter readout.

The design constraint that shapes everything here: **one task must not be able
to mutate the persistent model or another task's state.** So no module is ever
written to. The adapted parameters live in a plain dict of tensors ("fast
weights"), the readout is evaluated as a pure function of that dict, and the
module's own parameters are read but never assigned.

Why the readout can be evaluated outside the model at all: the adaptable scope
is the last layer of `interaction_head`, which sits downstream of everything
expensive. `encode` already computes `embed`, `section`, `occupancy`,
`ligand_value` and `protein_value` for support *and* query in a single pass, so
the inner loop needs no additional encoder work — it re-evaluates 97 parameters
on cached features. That makes k inner steps cost k cheap readouts rather than
k full forwards.

The re-implemented readout is a correctness risk, since it duplicates three
lines of `InteractionGrammarModel.encode`. `test_readout_reproduces_the_model`
pins it: with base weights the reconstruction must equal the model's own
`endpoint` bitwise. If `encode` ever changes, that test fails rather than the
experiment silently measuring the wrong thing.

Query labels are never read here. The inner loss is computed on support
endpoints against support labels only, which is the same information the
incumbent's `transport` already consumes.
"""
from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F
from torch import Tensor

# The frozen adaptable scope: the final linear readout of `interaction_head`.
# 97 parameters — 96 in the weight, 1 in the bias. The split is load-bearing:
# the weight rotates a linear functional and therefore reorders ligands, while
# the bias is exactly a scalar level shift. Reporting them apart is how a k=1
# "adaptation" gain is separated from a recalibration.
ADAPTABLE = ("interaction_head.2.weight", "interaction_head.2.bias")


@dataclass
class AdaptationConfig:
    """Everything about adaptation that a checkpoint must carry to reload."""
    inner_steps: int = 1
    inner_lr: float = 0.1
    first_order: bool = True
    scope: tuple[str, ...] = ADAPTABLE

    def to_dict(self) -> dict:
        return {"inner_steps": int(self.inner_steps),
                "inner_lr": float(self.inner_lr),
                "first_order": bool(self.first_order),
                "scope": list(self.scope)}

    @classmethod
    def from_dict(cls, payload: dict) -> "AdaptationConfig":
        return cls(inner_steps=int(payload["inner_steps"]),
                   inner_lr=float(payload["inner_lr"]),
                   first_order=bool(payload["first_order"]),
                   scope=tuple(payload["scope"]))


@dataclass
class EncodedTask:
    """One episode's cached features. Support and query, already split."""
    support_hidden: Tensor        # [B, K, embed+task] input to the readout
    query_hidden: Tensor          # [B, Q, embed+task]
    support_additive: Tensor      # [B, K] ligand_value + protein_value
    query_additive: Tensor        # [B, Q]
    support_occupancy: Tensor     # [B, K, contact_types]
    query_occupancy: Tensor       # [B, Q, contact_types]
    support_embed: Tensor         # [B, K, embed]  (for transport)
    query_embed: Tensor           # [B, Q, embed]


def base_weights(model) -> dict[str, Tensor]:
    """The persistent parameters of the adaptable scope, by name."""
    named = dict(model.named_parameters())
    missing = [name for name in ADAPTABLE if name not in named]
    if missing:
        raise KeyError(f"adaptable scope missing from the model: {missing}")
    return {name: named[name] for name in ADAPTABLE}


def encode_task(model, protein_pooled, protein_tokens, protein_mask,
                support_atoms, support_bonds, support_mask,
                query_atoms, query_bonds, query_mask,
                protein_chemistry) -> EncodedTask:
    """One encoder pass for support+query; everything the readout will need.

    This is the expensive half and it runs exactly once per task, before any
    inner step. Gradient flows through it to the shared initialization, which
    is what makes the outer update first-order MAML rather than a detached
    two-stage fit.
    """
    support_count = support_atoms.shape[1]
    raw_atoms = torch.cat((support_atoms, query_atoms), 1)
    bonds = torch.cat((support_bonds, query_bonds), 1)
    mask = torch.cat((support_mask, query_mask), 1)

    batch, count = raw_atoms.shape[:2]
    residues, summary = model.protein_encoder(
        protein_pooled, protein_tokens, protein_mask, protein_chemistry)
    residues = model.refine_slots(residues, protein_mask)
    ligand, atom_states = model.ligand_encoder(
        raw_atoms.flatten(0, 1), bonds.flatten(0, 1), mask.flatten(0, 1))
    residue_count = residues.shape[1]
    wide_residues = residues[:, None].expand(-1, count, -1, -1).reshape(
        batch * count, residue_count, -1)
    wide_mask = protein_mask[:, None].expand(-1, count, -1).reshape(
        batch * count, residue_count)
    occupancy, mean_state, max_state = model.grammar(
        atom_states, mask.flatten(0, 1),
        model.atom_chemistry(raw_atoms).flatten(0, 1), wide_residues, wide_mask)
    wide_summary = summary[:, None].expand(-1, count, -1).reshape(
        batch * count, -1)
    embed = model.embed_norm(model.embed(torch.cat(
        (ligand, mean_state, max_state, wide_summary, occupancy), -1)))
    section = model.section_norm(model.section(embed))
    ligand_value = model.ligand_head(ligand).squeeze(-1)
    protein_value = model.protein_head(wide_summary).squeeze(-1)

    shape = (batch, count)
    hidden = torch.cat((embed, section), -1).reshape(*shape, -1)
    additive = (ligand_value + protein_value).reshape(shape)
    occupancy = occupancy.reshape(*shape, -1)
    embed = embed.reshape(*shape, -1)
    split = (support_count, count - support_count)
    support_hidden, query_hidden = torch.split(hidden, split, 1)
    support_additive, query_additive = torch.split(additive, split, 1)
    support_occupancy, query_occupancy = torch.split(occupancy, split, 1)
    support_embed, query_embed = torch.split(embed, split, 1)
    return EncodedTask(
        support_hidden=support_hidden, query_hidden=query_hidden,
        support_additive=support_additive, query_additive=query_additive,
        support_occupancy=support_occupancy, query_occupancy=query_occupancy,
        support_embed=support_embed, query_embed=query_embed)


def readout(model, hidden: Tensor, additive: Tensor, occupancy: Tensor,
            weights: dict[str, Tensor]) -> Tensor:
    """`endpoint = additive + interaction_head(hidden) + contact_weight(occ)`.

    A pure function of `weights`: nothing is assigned to any module. The first
    layer of `interaction_head` and `contact_weight` are read from the model
    because they are outside the adaptable scope and never change during an
    inner loop.
    """
    body = model.interaction_head
    first = F.linear(hidden, body[0].weight, body[0].bias)
    activated = body[1](first)
    interaction = F.linear(activated,
                           weights["interaction_head.2.weight"],
                           weights["interaction_head.2.bias"]).squeeze(-1)
    contact = model.contact_weight(occupancy).squeeze(-1)
    return additive + interaction + contact


def adapt(model, task: EncodedTask, support_y: Tensor,
          config: AdaptationConfig,
          weights: dict[str, Tensor] | None = None
          ) -> tuple[dict[str, Tensor], list[float]]:
    """Support-only inner loop. Returns fast weights and the inner loss trace.

    `support_y` is the only label read. With `inner_steps == 0` or an empty
    support the base weights are returned **unchanged and unwrapped**, so the
    k=0 path is bitwise the zero-shot path rather than approximately it.

    First order: the inner gradient is detached, so `fast = base - lr * g` has
    `d(fast)/d(base) = I` and the outer gradient still reaches the shared
    initialization through both the cached features and this identity.
    """
    weights = dict(base_weights(model) if weights is None else weights)
    if config.inner_steps <= 0 or support_y.shape[-1] == 0:
        return weights, []
    # Adaptation *is* a gradient computation, so it cannot run inside the
    # ambient `no_grad` that evaluation uses. `enable_grad` is a no-op during
    # training and is what makes support adaptation legal at inference; the
    # fast weights are detached again on the way out so an evaluation loop
    # accumulates no graph.
    inference = not torch.is_grad_enabled()
    trace: list[float] = []
    with torch.enable_grad():
        for _ in range(int(config.inner_steps)):
            prediction = readout(model, task.support_hidden,
                                 task.support_additive,
                                 task.support_occupancy, weights)
            loss = F.mse_loss(prediction, support_y)
            trace.append(float(loss.detach()))
            gradients = torch.autograd.grad(
                loss, list(weights.values()),
                create_graph=(not config.first_order) and not inference,
                allow_unused=False)
            if config.first_order:
                # Standard first-order MAML: the inner gradient is treated as a
                # constant, so `d(fast)/d(base) = I`. This severs the autograd
                # path from `support_y` through the inner gradient as well —
                # support labels still *determine* the fast weights' value and
                # still reach the prediction through `transport`, but there is
                # no differentiable path through the adaptation itself. That is
                # a property of first-order adaptation, not an implementation
                # defect, and the tests assert the functional sensitivity here
                # and the gradient path in the second-order mode.
                gradients = [gradient.detach() for gradient in gradients]
            weights = {name: value - config.inner_lr * gradient
                       for (name, value), gradient
                       in zip(weights.items(), gradients)}
    if inference:
        weights = {name: value.detach() for name, value in weights.items()}
    return weights, trace


def partial_weights(model, adapted: dict[str, Tensor], keep: str
                    ) -> dict[str, Tensor]:
    """Fast weights with only one component of the update retained.

    `keep="bias"` is the pure scalar level shift; `keep="weight"` is the pure
    shape change. Running the same episode three ways is what separates a k=1
    "adaptation" gain from a recalibration, with no extra training.
    """
    base = base_weights(model)
    if keep == "both":
        return dict(adapted)
    if keep == "bias":
        return {"interaction_head.2.weight": base["interaction_head.2.weight"],
                "interaction_head.2.bias": adapted["interaction_head.2.bias"]}
    if keep == "weight":
        return {"interaction_head.2.weight": adapted["interaction_head.2.weight"],
                "interaction_head.2.bias": base["interaction_head.2.bias"]}
    if keep == "none":
        return dict(base)
    raise ValueError(f"unknown partial-weight mode {keep!r}")


def support_query_gradient_cosine(model, task: EncodedTask, support_y: Tensor,
                                  query_y: Tensor) -> float:
    """Agreement between the support and query gradients of the readout.

    The A2 selector's second term. Both gradients are taken with respect to the
    same 97 parameters at the shared initialization, so the cosine says whether
    fitting this task's support moves the readout in the direction its query
    also wants.

    Zero or unused gradients return exactly 0.0 rather than NaN: a task whose
    gradient vanishes carries no directional evidence, and 0.0 is the honest
    encoding of "no information", not a silent drop.

    Query labels enter **only here**, on `meta_train` tasks, to score a task
    for training-time sampling. They never reach a model input, an inference
    path, or a deployment decision.
    """
    weights = base_weights(model)
    values = list(weights.values())

    def gradient_of(hidden, additive, occupancy, truth):
        prediction = readout(model, hidden, additive, occupancy, weights)
        loss = F.mse_loss(prediction, truth)
        parts = torch.autograd.grad(loss, values, retain_graph=True,
                                    allow_unused=True)
        return torch.cat([
            (torch.zeros_like(value) if part is None else part).reshape(-1)
            for value, part in zip(values, parts)])

    if support_y.shape[-1] == 0 or query_y.shape[-1] == 0:
        return 0.0
    support_gradient = gradient_of(
        task.support_hidden, task.support_additive, task.support_occupancy,
        support_y)
    query_gradient = gradient_of(
        task.query_hidden, task.query_additive, task.query_occupancy, query_y)
    left = float(support_gradient.norm())
    right = float(query_gradient.norm())
    if left < 1e-12 or right < 1e-12:
        return 0.0
    value = float((support_gradient * query_gradient).sum() / (left * right))
    if value != value:                       # NaN guard, explicit
        return 0.0
    return max(-1.0, min(1.0, value))

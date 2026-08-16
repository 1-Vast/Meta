"""Deterministic, label-safe SIMA-DTA episode construction."""

from __future__ import annotations

from collections import defaultdict
from typing import Mapping, Sequence

import torch

from .contract import AffinityRow, Episode


def _device() -> torch.device:
    if not torch.cuda.is_available():
        raise RuntimeError("SIMA-DTA numerical protocol requires CUDA")
    return torch.device("cuda")


def _vector(value: Sequence[float]) -> torch.Tensor:
    result = torch.as_tensor(value, dtype=torch.float64, device=_device())
    if result.ndim != 1 or result.numel() == 0 or not torch.isfinite(result).all():
        raise ValueError("support-design vectors must be finite one-dimensional arrays")
    return result


def selectsupport(
    rows: Sequence[AffinityRow],
    candidate_indices: Sequence[int],
    support_size: int,
    design_vectors: Mapping[int, Sequence[float]],
    *,
    chemical_components: Mapping[int, str] | None = None,
    chemical_component_cap: int = 1,
    ridge: float = 1e-3,
) -> tuple[int, ...]:
    """Greedily select support from covariates only; affinity values are never read."""

    if support_size <= 0:
        raise ValueError("support_size must be positive")
    if chemical_component_cap <= 0:
        raise ValueError("chemical_component_cap must be positive")
    candidates = tuple(
        sorted(
            candidate_indices,
            key=lambda index: (
                rows[index].scaffold_key,
                rows[index].ligand_parent_key,
                index,
            ),
        )
    )
    if len(candidates) < support_size:
        raise ValueError("insufficient candidates for requested support size")
    if any(index not in design_vectors for index in candidates):
        raise ValueError("every candidate needs a label-free design vector")

    matrix = torch.stack([_vector(design_vectors[index]) for index in candidates])
    dimension = matrix.shape[1]
    inverse = torch.eye(dimension, dtype=torch.float64, device=_device()) / ridge
    positions = {index: position for position, index in enumerate(candidates)}
    selected: list[int] = []
    selected_scaffolds: set[str] = set()
    selected_components: dict[str, int] = defaultdict(int)
    for _ in range(support_size):
        options: list[int] = []
        for index in candidates:
            if index in selected or rows[index].scaffold_key in selected_scaffolds:
                continue
            component = chemical_components[index] if chemical_components else ""
            if chemical_components and selected_components[component] >= chemical_component_cap:
                continue
            options.append(index)
        if not options:
            raise ValueError("scaffold or chemical-neighbour constraints exhaust support candidates")
        remaining = [index for index in candidates if index not in selected]
        if len(remaining) <= 1:
            raise ValueError("support selection leaves no query candidates")
        remaining_matrix = matrix[[positions[index] for index in remaining]]
        option_matrix = matrix[[positions[index] for index in options]]
        transformed = option_matrix @ inverse
        quadratic = (transformed * option_matrix).sum(dim=1)
        query_gram = remaining_matrix.T @ remaining_matrix
        baseline = torch.trace(inverse @ query_gram)
        reduction = ((transformed @ query_gram) * transformed).sum(dim=1) / (1.0 + quadratic)
        own = quadratic / (1.0 + quadratic)
        scores = (baseline - reduction - own) / (len(remaining) - 1)
        chosen = min(
            zip(scores.detach().cpu().tolist(), options),
            key=lambda item: (
                item[0],
                rows[item[1]].scaffold_key,
                rows[item[1]].ligand_parent_key,
                item[1],
            ),
        )[1]
        vector = matrix[positions[chosen]]
        transformed_vector = inverse @ vector
        inverse = inverse - torch.outer(transformed_vector, transformed_vector) / (
            1.0 + vector @ transformed_vector
        )
        selected.append(chosen)
        selected_scaffolds.add(rows[chosen].scaffold_key)
        if chemical_components:
            selected_components[chemical_components[chosen]] += 1
    return tuple(sorted(selected))


def buildregistry(
    rows: Sequence[AffinityRow],
    *,
    homology_components: Mapping[str, str],
    design_vectors: Mapping[int, Sequence[float]],
    support_size: int,
    chemical_components: Mapping[int, str] | None = None,
    chemical_component_cap: int = 1,
) -> tuple[Episode, ...]:
    """Build one deterministic strict scaffold-cold episode per target stratum."""

    by_target: dict[str, list[int]] = defaultdict(list)
    for index, row in enumerate(rows):
        if row.split_role != "meta_test":
            continue
        by_target[row.target_key].append(index)

    episodes: list[Episode] = []
    for target_key in sorted(by_target):
        candidate_indices = by_target[target_key]
        endpoints = {rows[index].endpoint for index in candidate_indices}
        if len(endpoints) != 1:
            raise ValueError("filter endpoint strata before building episodes")
        support = selectsupport(
            rows,
            candidate_indices,
            support_size,
            design_vectors,
            chemical_components=chemical_components,
            chemical_component_cap=chemical_component_cap,
        )
        support_scaffolds = {rows[index].scaffold_key for index in support}
        support_components = (
            {chemical_components[index] for index in support} if chemical_components else set()
        )
        query = tuple(
            index
            for index in sorted(candidate_indices)
            if index not in support
            and rows[index].scaffold_key not in support_scaffolds
            and (not chemical_components or chemical_components[index] not in support_components)
        )
        if not query:
            raise ValueError("strict support/query closure leaves no query rows")
        episode = Episode.create(
            target_key=target_key,
            support_indices=support,
            query_indices=query,
            homology_component=homology_components[target_key],
            support_scaffolds=[rows[index].scaffold_key for index in support],
            query_scaffolds=[rows[index].scaffold_key for index in query],
            provenance_components=[
                rows[index].document_or_provenance_key for index in support + query
            ],
        )
        episode.assertrows(rows)
        episodes.append(episode)
    return tuple(episodes)


def buildwrong(
    episode: Episode,
    rows: Sequence[AffinityRow],
) -> tuple[int, ...]:
    """Match support chemistry to a different target without touching affinity values."""

    wrong: list[int] = []
    for support_index in episode.support_indices:
        ligand_key = rows[support_index].ligand_parent_key
        candidates = [
            index
            for index, row in enumerate(rows)
            if row.target_key != episode.target_key
            and row.ligand_parent_key == ligand_key
            and row.endpoint == rows[support_index].endpoint
        ]
        if not candidates:
            raise ValueError("cannot construct chemistry-matched wrong-target support")
        wrong.append(min(candidates))
    if len({rows[index].target_key for index in wrong}) == 0:
        raise ValueError("wrong-target support is empty")
    return tuple(wrong)

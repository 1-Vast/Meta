"""Deployment-coverage audit: could MMP ever serve as an inference mechanism?

T1 asks whether the transformation graph is identifiable for *estimation*. This
asks a different and equally decisive question: at deployment, when a governed
episode hands the model k support ligands and a query, how often does any
support-query pair even form a valid MMP?

    C_k = P(at least one support-query pair forms a valid MMP)

Estimated on the frozen, label-blind nested k = {1,2,3,5} episode banks that the
governed protocol itself uses. **No query label is read** -- the MMP relation is
a function of structure only.

A low C_k does not invalidate the transformation space as a *training signal*;
it does mean MMP cannot be a universal reference-based inference mechanism, and
the report must say so.

Run:
    python -m tools.research.stageT_mmp.t1_coverage
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

import numpy as np

if __package__ in {None, ""}:  # pragma: no cover
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from scripts.internal_validation import partition_components
from tools.research.stageT_mmp.mmp import fragment, transformation
from tools.research.stageT_mmp.observations import load_governed

HERE = Path(__file__).resolve().parent

SUPPORT_SIZES = (1, 2, 3, 5)
QUERY_SIZE = 16
DRAWS = 1
BANK_SEED = 20260820
# Frozen novelty cut, matching the Stage S Phase 0 quantiles.
NOVELTY_TERCILES = (0.30136987566947937, 0.5606504082679749)


def _pairs_form_mmp(data, support: tuple[int, ...], query_index: int,
                    coarse: bool) -> bool:
    """Does any support ligand form a single-cut MMP with the query ligand?"""
    query_smiles = data._ligand_smiles.get(data.cells[query_index]["ligand_id"])
    if not query_smiles:
        return False
    query_pieces = fragment(query_smiles)
    if not query_pieces:
        return False
    query_by_core = defaultdict(list)
    for piece in query_pieces:
        query_by_core[piece.core].append(piece)
    for cell in support:
        smiles = data._ligand_smiles.get(data.cells[cell]["ligand_id"])
        if not smiles:
            continue
        for piece in fragment(smiles):
            for other in query_by_core.get(piece.core, ()):
                built = transformation(piece, other)
                if built is not None:
                    return True
        if coarse:
            # The coarse relation only relaxes the attachment context, so it is
            # a superset; the exact test above already covers the strict case.
            for piece in fragment(smiles):
                for other in query_pieces:
                    if (piece.core == other.core
                            and piece.r_group != other.r_group
                            and piece.coarse_context == other.coarse_context):
                        return True
    return False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "T1_COVERAGE.json")
    args = parser.parse_args()

    data, seal = load_governed()
    if not seal["isolation"]["physically_isolated"]:
        raise SystemExit("refusing to audit without the physical split view")
    fit, internal = partition_components(data)
    component_of = {cell["target_id"]: cell["protein_group_40"]
                    for cell in data.cells}
    internal_set = set(internal)

    banks = data.fixed_nested_episode_banks(
        "meta_train", SUPPORT_SIZES, QUERY_SIZE, DRAWS, BANK_SEED)

    # Ligand novelty against the fit components, label-blind.
    fit_ligands = sorted({cell["ligand_id"] for cell in data.cells
                          if cell["protein_group_40"] in set(fit)})
    table = data.fingerprints
    reference = np.stack([table[key].numpy() for key in fit_ligands])
    reference_counts = reference.sum(axis=1)
    novelty_cache: dict[str, float] = {}

    def novelty(ligand_id: str) -> float:
        if ligand_id not in novelty_cache:
            row = table[ligand_id].numpy()
            intersection = reference @ row
            union = reference_counts + row.sum() - intersection
            with np.errstate(divide="ignore", invalid="ignore"):
                value = np.where(union > 0,
                                 intersection / np.maximum(union, 1e-12), 0.0)
            novelty_cache[ligand_id] = float(value.max()) if value.size else 0.0
        return novelty_cache[ligand_id]

    report: dict = {
        "schema": "MetaSieve.StageT.T1Coverage.v1",
        "definition": ("C_k = P(at least one support-query pair forms a valid "
                       "single-cut MMP), on the frozen nested episode banks"),
        "bank": {"support_sizes": list(SUPPORT_SIZES), "query_size": QUERY_SIZE,
                 "draws": DRAWS, "seed": BANK_SEED,
                 "constructor": "QPSMPData.fixed_nested_episode_banks"},
        "labels_read": False,
        "meta_test": seal,
        "coverage": {},
    }

    for size in SUPPORT_SIZES:
        exact_hits: list[float] = []
        coarse_hits: list[float] = []
        by_component: dict[str, list[float]] = defaultdict(list)
        by_novelty: dict[str, list[float]] = defaultdict(list)
        by_population: dict[str, list[float]] = defaultdict(list)
        for spec in banks[size]:
            for query_index in spec.query:
                exact = _pairs_form_mmp(data, spec.support, query_index, False)
                coarse = exact or _pairs_form_mmp(
                    data, spec.support, query_index, True)
                exact_hits.append(float(exact))
                coarse_hits.append(float(coarse))
                by_component[spec.component].append(float(exact))
                value = novelty(data.cells[query_index]["ligand_id"])
                bucket = ("novelty_low" if value < NOVELTY_TERCILES[0]
                          else "novelty_mid" if value < NOVELTY_TERCILES[1]
                          else "novelty_high")
                by_novelty[bucket].append(float(exact))
                population = ("internal" if component_of[spec.target] in internal_set
                              else "fit")
                by_population[population].append(float(exact))
        component_means = [float(np.mean(v)) for v in by_component.values() if v]
        report["coverage"][str(size)] = {
            "queries_scored": len(exact_hits),
            "C_k_exact": float(np.mean(exact_hits)) if exact_hits else 0.0,
            "C_k_coarse": float(np.mean(coarse_hits)) if coarse_hits else 0.0,
            "C_k_exact_component_equal_weight": (
                float(np.mean(component_means)) if component_means else 0.0),
            "components": len(by_component),
            "by_novelty": {name: float(np.mean(v))
                           for name, v in sorted(by_novelty.items())},
            "by_population": {name: float(np.mean(v))
                              for name, v in sorted(by_population.items())},
        }

    values = [report["coverage"][str(s)]["C_k_exact"] for s in SUPPORT_SIZES]
    report["interpretation"] = {
        "max_C_k": max(values),
        "verdict": (
            "MMP can serve as a reference-based inference mechanism for at most "
            f"{max(values):.1%} of governed queries at k<=5. Below that share it "
            "is a TRAINING signal, not a universal deployment mechanism, and no "
            "artifact may present it as one."),
    }

    args.output.write_text(json.dumps(report, indent=1, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({"output": str(args.output),
                      "C_k_exact": {s: round(report["coverage"][str(s)]["C_k_exact"], 4)
                                    for s in SUPPORT_SIZES},
                      "C_k_coarse": {s: round(report["coverage"][str(s)]["C_k_coarse"], 4)
                                     for s in SUPPORT_SIZES}}, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

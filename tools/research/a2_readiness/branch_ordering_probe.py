"""A2 readiness probe: does frozen A0 carry protein-conditioned SAR ordering?

No training, no gradient, forward passes only. `meta_test` is physically
sealed by an explicit `include_meta_test=False`.

## Why this probe decides the A2 family

`NEXT_RESEARCH_PLAN_A2_MOMENT_20260816.md` rests on one premise: that the
frozen A0 complex representation `e0(P,L)` contains a *protein-conditioned*
SAR coordinate system that at most five support labels can identify. A2-min
builds `z = A_phi(e0)` and forms `c_S = mean_i r_i z(P,L_i)`, then predicts
`delta = eta(k) <c_S, z(P,L_q)>`.

If the frozen endpoint's within-target ordering is entirely explained by its
**protein-blind** ligand branch, then the protein conditioning in `e0` does
not vary usefully across ligands of one target, and A2-min would be
identifying a coordinate in a subspace that carries no protein-specific
ordering signal. The premise would be false before any training.

A0's endpoint decomposes exactly (`model/interaction_grammar.py::encode`):

    endpoint = ligand_value(L) + protein_value(P) + interaction(P, L)

* `protein_value` depends only on the protein, so it is **constant within a
  target** and drops out of any centered within-target statistic. It cannot
  affect ordering at all.
* `ligand_value` is protein-blind: the same ligand receives the same value
  for every target. It can order ligands, but identically for all proteins.
* `interaction` is the only branch that can order ligands *differently* for
  different proteins.

So the probe measures, per target, the centered correlation with truth of:

  `full`         the endpoint actually deployed;
  `ligand_only`  `ligand_value` alone (the protein-blind ordering);
  `interaction`  `endpoint - ligand_value - protein_value`.

and the **incremental** ordering the interaction branch supplies over the
protein-blind branch. Three quantities decide the premise:

1. `r(full)` vs `r(ligand_only)` — if equal, the interaction branch adds no
   ordering;
2. `r(interaction)` alone — whether the branch orders at all;
3. the **protein counterfactual**: recompute the interaction branch with a
   similarity-matched wrong protein from the same split. If the interaction
   ordering is unchanged, whatever it encodes is not protein-specific.

Point 3 is the control that separates "the interaction branch orders" from
"the interaction branch orders *because of this protein*", and it is the
exact failure mode that made seven query-specific channels inert.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
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
    component_bootstrap, component_target_mean,
)
from scripts.stageR6_compare_arms import SUPPORT_SIZES            # noqa: E402
from tools.research.a2_readiness._arms import build_arm           # noqa: E402
from scripts.train_level_shape import matched_donors, normalized  # noqa: E402
from scripts.train_qpsmp import (                                 # noqa: E402
    COMPACT_LIGAND_BANK, CORPUS, LIGAND_BANK, PROTEIN_BANK,
    compact_episode, training_label_scale,
)


def centered_correlation(prediction: np.ndarray, truth: np.ndarray) -> float:
    """Within-target Pearson correlation; the quantity R14 localized."""
    p = prediction - prediction.mean()
    y = truth - truth.mean()
    denominator = float(np.sqrt((p ** 2).mean()) * np.sqrt((y ** 2).mean()))
    return float((p * y).mean() / denominator) if denominator > 1e-12 else 0.0


def concordance(prediction: np.ndarray, truth: np.ndarray) -> float:
    rows, cols = np.triu_indices(len(truth), 1)
    delta_y = truth[rows] - truth[cols]
    comparable = delta_y != 0
    if not comparable.any():
        return float("nan")
    signed = np.sign(delta_y[comparable]) * (
        prediction[rows] - prediction[cols])[comparable]
    return float((signed > 0).mean() + 0.5 * (signed == 0).mean())


def protein_inputs(data: QPSMPData, target: str, device: str, dtype):
    pooled, tokens, mask = data.protein_for_target(target)
    chemistry = data.protein_chemistry_for_target(target)
    return (pooled.to(device, dtype).unsqueeze(0),
            tokens.to(device, dtype).unsqueeze(0),
            mask.to(device, dtype).unsqueeze(0),
            chemistry.to(device, dtype).unsqueeze(0))


def branch_values(model, data, spec, episode, protein_target, device, dtype):
    """Zero-shot branch decomposition for one episode's queries.

    `adapt=False` makes this exactly the k=0 endpoint, so no support label is
    read and the transport contributes nothing.
    """
    atoms = torch.cat((episode.support_atoms, episode.query_atoms), 0)
    bonds = torch.cat((episode.support_bonds, episode.query_bonds), 0)
    mask = torch.cat((episode.support_mask, episode.query_mask), 0)
    pooled, tokens, protein_mask, chemistry = protein_inputs(
        data, protein_target, device, dtype)
    output = model(
        pooled, tokens, protein_mask,
        atoms[:0].to(device, dtype).unsqueeze(0),
        bonds[:0].to(device, dtype).unsqueeze(0),
        mask[:0].to(device, dtype).unsqueeze(0),
        episode.support_y[:0].to(device, dtype).unsqueeze(0),
        episode.query_atoms.to(device, dtype).unsqueeze(0),
        episode.query_bonds.to(device, dtype).unsqueeze(0),
        episode.query_mask.to(device, dtype).unsqueeze(0),
        adapt=False, protein_chemistry=chemistry,
        support_fingerprint=episode.support_fingerprint[:0].to(device, dtype).unsqueeze(0),
        query_fingerprint=episode.query_fingerprint.to(device, dtype).unsqueeze(0))
    endpoint = output.zero_shot.squeeze(0).float().cpu().numpy()
    ligand = output.ligand_only.squeeze(0).float().cpu().numpy()
    additive = output.additive.squeeze(0).float().cpu().numpy()
    # additive = ligand_value + protein_value, so interaction = endpoint - additive
    # and protein_value = additive - ligand_value (constant within a target).
    return {"full": endpoint, "ligand_only": ligand,
            "protein_only": additive - ligand,
            "interaction": endpoint - additive}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm", action="append", required=True,
                        help="name=path/to/checkpoint.pt (repeat for seeds)")
    parser.add_argument("--split-directory", type=Path, required=True)
    parser.add_argument("--evaluation-seed", type=int, default=73101)
    parser.add_argument("--query-size", type=int, default=16)
    parser.add_argument("--bootstrap-draws", type=int, default=9999)
    parser.add_argument("--device", default="cuda")
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()

    # The seal is opt-out in QPSMPData; state it explicitly here.
    data = QPSMPData(CORPUS, PROTEIN_BANK, LIGAND_BANK, COMPACT_LIGAND_BANK,
                     split_directory=arguments.split_directory,
                     include_meta_test=False)
    scale = training_label_scale(data)
    donors = matched_donors(data, "meta_val", donor_pool="meta_val",
                            whitening_pool="meta_train")
    banks = data.fixed_nested_episode_banks(
        "meta_val", SUPPORT_SIZES, arguments.query_size, 1,
        arguments.evaluation_seed, None)
    specs = banks[0]                      # k=0 only: this is a zero-shot probe

    rows: list[dict] = []
    for item in arguments.arm:
        name, _, path = item.partition("=")
        model, kind, seed = build_arm(path, name, data, arguments.device)
        dtype = next(model.parameters()).dtype
        with torch.no_grad():
            for spec in specs:
                episode = compact_episode(
                    normalized(data.materialize(spec), scale))
                truth = (episode.query_y.numpy() * scale.scale + scale.mean)
                correct = branch_values(model, data, spec, episode,
                                        spec.target, arguments.device, dtype)
                wrong = branch_values(model, data, spec, episode,
                                      donors[spec.target], arguments.device, dtype)
                row = {"arm": name, "seed": seed, "target": spec.target,
                       "component": spec.component}
                for branch, values in correct.items():
                    row[f"r_{branch}"] = centered_correlation(
                        values * scale.scale, truth)
                    row[f"ci_{branch}"] = concordance(values, truth)
                    row[f"sd_{branch}"] = float(
                        (values * scale.scale).std())
                row["r_interaction_wrong_protein"] = centered_correlation(
                    wrong["interaction"] * scale.scale, truth)
                row["r_full_wrong_protein"] = centered_correlation(
                    wrong["full"] * scale.scale, truth)
                # How much the interaction branch moves when only the protein
                # changes; zero means the branch is protein-inert for ordering.
                interaction_shift = (correct["interaction"] - wrong["interaction"])
                centered_shift = interaction_shift - interaction_shift.mean()
                row["interaction_protein_shift_sd_pk"] = float(
                    (centered_shift * scale.scale).std())
                # POSITIVE CONTROL. A null on ordering is only interpretable if
                # the protein swap demonstrably changed something. These are
                # uncentered magnitudes: the level must move a lot, otherwise
                # the donor was not actually substituted and the null is a bug.
                row["level_protein_shift_pk"] = float(abs(
                    (correct["protein_only"] - wrong["protein_only"]).mean()
                    * scale.scale))
                row["endpoint_protein_shift_pk"] = float(abs(
                    (correct["full"] - wrong["full"]).mean() * scale.scale))
                row["interaction_uncentered_shift_pk"] = float(abs(
                    interaction_shift.mean() * scale.scale))
                row["sd_truth"] = float(truth.std())
                rows.append(row)
        del model
        if arguments.device.startswith("cuda"):
            torch.cuda.empty_cache()

    def weighted(subset: list[dict], field: str) -> float:
        return component_target_mean(
            (r["component"], r["target"], r.get(field)) for r in subset)

    report: dict = {
        "schema": "MetaSieve.A2ReadinessBranchOrdering.v1",
        "split": "meta_val", "k": 0,
        "split_directory": str(arguments.split_directory),
        "split_assignment_sha256": data.split_manifest["assignment_sha256"],
        "meta_test": {"evaluated": False, "included": False,
                      "seal": "physical: QPSMPData(include_meta_test=False)"},
        "population": {"targets": len(specs), "query_size": arguments.query_size},
        "arms": {},
    }
    fields = ["r_full", "r_ligand_only", "r_interaction",
              "ci_full", "ci_ligand_only", "ci_interaction",
              "sd_full", "sd_ligand_only", "sd_interaction",
              "r_interaction_wrong_protein", "r_full_wrong_protein",
              "interaction_protein_shift_sd_pk", "level_protein_shift_pk",
              "endpoint_protein_shift_pk", "interaction_uncentered_shift_pk",
              "sd_truth"]
    for name in sorted({r["arm"] for r in rows}):
        subset = [r for r in rows if r["arm"] == name]
        cell = {field: weighted(subset, field) for field in fields}
        # Paired contrasts on the two questions that decide the A2 premise.
        cell["interaction_ordering_increment"] = component_bootstrap(
            [(r["component"], r["target"], r["r_full"] - r["r_ligand_only"])
             for r in subset], arguments.bootstrap_draws, 20260816)
        cell["interaction_protein_specificity"] = component_bootstrap(
            [(r["component"], r["target"],
              r["r_interaction"] - r["r_interaction_wrong_protein"])
             for r in subset], arguments.bootstrap_draws, 20260816)
        report["arms"][name] = cell

    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    arguments.output.write_text(json.dumps(report, indent=1), encoding="utf-8")
    with arguments.output.with_suffix(".rows.jsonl").open(
            "w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row) + "\n")

    print(f"{'arm':<10}{'r_full':>9}{'r_lig':>9}{'r_inter':>9}"
          f"{'sd_lig':>9}{'sd_inter':>10}{'r_int_wrongP':>14}")
    for name, cell in report["arms"].items():
        print(f"{name:<10}{cell['r_full']:>9.3f}{cell['r_ligand_only']:>9.3f}"
              f"{cell['r_interaction']:>9.3f}{cell['sd_ligand_only']:>9.3f}"
              f"{cell['sd_interaction']:>10.3f}"
              f"{cell['r_interaction_wrong_protein']:>14.3f}")
    print("\npaired component bootstrap:")
    for name, cell in report["arms"].items():
        inc = cell["interaction_ordering_increment"]
        spec = cell["interaction_protein_specificity"]
        print(f"  {name}")
        print(f"    r(full) - r(ligand_only)        {inc['mean']:+.4f} "
              f"[{inc['lo']:+.4f}, {inc['hi']:+.4f}]  "
              f"{'RESOLVED' if inc['lo'] > 0 else 'unresolved'}")
        print(f"    r(inter) - r(inter|wrong P)     {spec['mean']:+.4f} "
              f"[{spec['lo']:+.4f}, {spec['hi']:+.4f}]  "
              f"{'RESOLVED' if spec['lo'] > 0 else 'unresolved'}")
    print("\npositive control — did the protein swap change anything at all?")
    for name, cell in report["arms"].items():
        print(f"  {name}: level shift {cell['level_protein_shift_pk']:.4f} pK | "
              f"endpoint shift {cell['endpoint_protein_shift_pk']:.4f} pK | "
              f"interaction shift (uncentered) "
              f"{cell['interaction_uncentered_shift_pk']:.4f} pK | "
              f"interaction shift (centered) "
              f"{cell['interaction_protein_shift_sd_pk']:.4f} pK")
    print(f"\nwrote {arguments.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

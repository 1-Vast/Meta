# A2S-CFES C0B structural semantic gate decision

Date: 2026-08-02
Branch at execution: `research/a2s-conformational-free-energy-state-20260802`
Runner: `research/a2s_cfes_semantic_gate.py`
Artifacts: `reports/active/a2s_cfes_semantic_gate_2026-08-02.json`,
`..._records_2026-08-02.parquet`, `..._weights_2026-08-02.pt`

**FACT — decision: `CFES_C0B_SEMANTICS_NOT_ADMITTED_STOP_CFES`.**

This closes the conformational free-energy state branch at its first real-data
gate, exactly as the preregistration required. No affinity label was opened;
the gate ran entirely on outcome-free PLINDER structural annotations.

## What was asked

C0B asked the narrowest necessary question behind CFES: on protein-, pocket-,
ligand- and provenance-disjoint structural splits, does an explicit
ligand-by-pocket term predict the observed eight-type contact profile better
than matched additive and protein-free predictors, and is that increment
destroyed when the physical pairing is broken?

## Binding results (217 audit clusters, 5 seeds)

| Contrast | Mean | 95% interval | Required | Pass |
|---|---:|---:|---|---|
| `cross` − `additive` | +0.00172 | [+0.00040, +0.00323] | lower > 0 | yes |
| `cross` − `no_cross_capacity` | −0.00177 | [−0.00711, +0.00397] | lower > 0 | **no** |
| `cross` − `frozen_random_cross` | −0.00017 | [−0.00103, +0.00074] | lower > 0 | **no** |
| `cross` − `ligand_only` | **−0.09166** | [−0.14337, −0.04675] | lower > 0 | **no** |

The multiplicative term does beat the additive model, but it does not beat a
parameter-matched residual without a multiplicative term, and it does not beat
the same bilinear interface with **frozen random projections**. The increment is
therefore capacity, not physics.

The decisive number is the last row. A ligand-only predictor beats the full
ligand-plus-pocket cross model by 0.092 on the standardized contact loss. Adding
the pocket representation actively *hurts* prediction of the pocket's own contact
profile.

## Destruction controls

| Destruction | Destroyed effect | Removal fraction | Required |
|---|---:|---:|---|
| `pocket_shuffle` | −0.00025 | 1.14 | ≥ 0.70 |
| `ligand_shuffle` | −0.00008 | 1.05 | ≥ 0.70 |
| `structure_transplant` | +0.00012 | 0.93 | ≥ 0.70 |
| `residue_randomization` | +0.00137 | 0.20 | ≥ 0.70 |

Removal fractions near or above 1.0 do not confirm the mechanism here. They say
the tiny surviving increment is indistinguishable from noise in both the intact
and the destroyed arm — and `residue_randomization` removes only 20 %, i.e.
scrambling the amino-acid identities of the pocket barely changes the result.

Consistency also failed: the cross-minus-additive effect was negative in audit
fold 3 (−0.00076) and in seed 1733 (−0.00084), and only 3 of 8 contact
coordinates improved, with **hydrogen bonds (−0.00005) and hydrophobic contacts
(+0.00013) both essentially null** — the two coordinates the preregistration
named as mandatory. The improving coordinates were metal complexes, salt bridges
and water bridges, which are rare and dominated by a few systems.

## Harness validity

**FACT.** The synthetic positive control passed: with a planted rank-4
ligand-by-pocket term at the measured dimensions, sample count and noise scale,
the cross-minus-additive lower bound was positive, all five audit folds were
positive, and pocket/ligand shuffling each removed ≥ 70 % of the effect. The
state-duplication and state-order-permutation invariances were exact no-ops. The
firewall check passed: no affinity column, no PLINDER `test`/`removed` row, and
all 467 A2S accessions excluded.

So the negative is a measurement, not a broken harness.

## Interpretation

**INFERENCE.** Pocket composition — normalized neighbouring-residue identity,
chemistry groups, size and coarse sequence position — carries no transferable
ligand-conditional interaction semantics beyond what the ligand alone implies.
The contact profile of a protein–ligand complex is, in this representation,
predictable from the ligand.

This is the same shape of failure as Gate G3 (pooled ESM-2 predicts the
per-target head at −0.019) and PIRS (learned protein-conditioned coordinates not
admitted), now measured on *experimental structures* rather than sequence. The
CFES premise required state-specific scoring to be a real, transferable
capability before any population logit could be inferred from k ≤ 5 labels. It
is not available at the coarse pocket-composition level, and the preregistration
forbids rescuing it with more capacity, rank, or epochs.

## Consequences

1. **FACT.** CFES stops. Gates C1, C2 and C3 are not run. No conformational
   population operator is trained.
2. **FACT.** Per the branch charter's stop rule, the successor must start from a
   different question, not a larger structural encoder.
3. **INFERENCE.** Three independent representations of the protein — pooled
   sequence (G3), learned sequence segments (PIRS), and experimental pocket
   composition (C0B) — have now each failed to carry target-conditional
   interaction information on this substrate. Protein-side conditioning is
   retired as a *source* of the adaptation object until some representation
   demonstrates incremental value on a label-free semantic task first.
4. **FACT.** Nothing is promoted to `model/` or `script/`. No breakthrough.

## Successor

The successor branch is `research/a2s-transfer-object-20260802`. It does not
propose another representation. It first measures the two premises that all nine
prior mechanisms assumed and none tested: whether the per-target headroom is a
chemical object at all, and whether the adaptation object transfers between
targets. See `A2S_TRANSFER_OBJECT_GATE_T0_DECISION_2026-08-02.md`.

Artifact hashes:

- JSON content: see `content_sha256` in the artifact
- verdict block: `{"pass": false, "verdict": "CFES_C0B_SEMANTICS_NOT_ADMITTED_STOP_CFES"}`

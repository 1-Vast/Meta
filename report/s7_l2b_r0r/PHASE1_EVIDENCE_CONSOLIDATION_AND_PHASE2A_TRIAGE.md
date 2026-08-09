# S7/L2B Phase 1 evidence consolidation and Phase 2A triage

Date: 2026-08-10.

> **SUPERSEDED IN PART BY PHASE 2A (2026-08-10).** Every number below stands and
> nothing is rewritten. Two *interpretations* below are corrected by
> `PHASE2A_SYNTHESIS.md` / `PHASE2A_VERDICT.json`:
>
> 1. `GENERIC_POCKET_SIGNAL_DOMINANT` and `LIGAND_CONDITIONING_WEAK` describe
>    **B5**, not the corpus. Against a real alternative ligand of the same exact
>    construct and a replicate noise floor measured from the data, the **labels**
>    are strongly ligand-conditioned: ΔJ `+0.258` [LCB `+0.234`], Spearman
>    ρ `+0.322`. `TEACHER_GENERIC_POCKET_ONLY` is refuted.
> 2. The decision rule "teacher ligand-conditioned increment weak → stop" did not
>    fire. The executed verdict is
>    `LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING`, whose mandated
>    action is one ligand-conditioned **residue residual** head — not the
>    "residue-first typed mechanism field" sketched in this document, and not a
>    pair-coupling head.
>
> The 92.5% wrong-ligand retention quoted here used an *arbitrary foreign*
> ligand, which Phase 2A classifies as a corruption control rather than a
> biological negative.

This document is the compact current source of truth for the completed S7/L2B
Phase 0 integrity work and Phase 1 frozen-ESM2 B5 discriminator. It does not
replace the immutable machine-readable artifacts listed below.

## Current verdict

```text
PHASE0_INTEGRITY_CONTRACT_PASS
B5_FROZEN_ESM2_DEVELOPMENT_GATE_PASS_6_OF_6
EXACT_RESIDUE_LOCALISATION_IDENTIFIED_IN_DEVELOPMENT
GENERIC_POCKET_SIGNAL_DOMINANT
LIGAND_CONDITIONING_WEAK
EXACT_RESIDUE_ATOM_COUPLING_NOT_IDENTIFIED
AFFINITY_DIRECTION_NOT_TESTED
K_SHOT_SECTION_NOT_IDENTIFIED
BIOLOGICAL_STATISTIC_NOT_ADMITTED_TO_Z
PHASE2A_AUDIT_NOT_REGISTERED
```

The B5 result is a structural development PASS. It is not an affinity result,
an independent time-forward confirmation, or an admission into the production
statistic `z`.

## Chronology and integrity

| Item | Frozen evidence |
|---|---|
| Unified preregistration | commit `ce186f4` |
| Phase 0 integrity artifacts | commit `139effd`, before B5 scoring |
| Phase 1 B5 result | commit `623602e` |
| Atom contract | 375,311 positions checked; 14,585 admitted; 4 quarantined |
| Sealed predictions | 52,062,975 held-out cells per arm, hashed float16 |
| Determinism | same-data state dictionaries bit-identical |
| ESM2 | `esm2_t33_650M_UR50D`, revision `08e4846e`, weight SHA-256 `c874668852c7275a159e2c7ceb6069671d7b1ba2c7b52f59600b34ce0f721008` |
| Label firewall | no ChEMBL, BindingDB, DAVIS or recipient affinity value read |

## Phase 1 result

Held-out A contains 2,409 complexes in 196 protein closure components. B5 and
B4 differ only in residue features; the atom branch, pair head, rank, sampler,
optimizer, budget, split, evaluation mask and tie policy are matched.

| Arm | Tie-aware component-macro AP |
|---|---:|
| B0 prevalence | 0.003185 |
| BM5 motif shuffle | 0.004508 |
| BP5 wrong protein | 0.004640 |
| BL ligand-only | 0.005719 |
| BX5 wrong ligand | 0.019679 |
| B4 non-PLM residue | 0.023251 |
| **B5 frozen ESM2** | **0.069601** |

All six registered contrasts pass. The smallest one-sided lower bound is
`+0.040392` for B5 minus B4, above the frozen `0.02` practical margin.

## What B5 identified

The sealed marginal decomposition localizes the improvement to residues:

```text
residue marginal: B5 - B4 = +0.177221, LCB95 = +0.160064
atom marginal:    B5 - B4 = -0.009949, interval spans zero
```

Changing only the residue representation therefore resolves the earlier B4
representation deficit for residue localization. It does not establish an
affinity mechanism.

Wrong-ligand residue AP is `0.245328` against B5 `0.265114`; approximately
92.5% of residue localization survives the ligand substitution. B5 mainly
identifies a generic protein pocket. Pair-level AP can still be produced by an
additive residue-pocket term plus ligand-atom propensity, so the six-Gate PASS
does not by itself identify exact residue-atom coupling.

## Data and confirmation boundary

MONN supports a nearly document-closed cohort after removing four shared
documents. It cannot provide a meaningful time-forward confirmation: only two
additional-PDB entries satisfy the frozen 2019-01-01 cutoff and none satisfy a
2024-01-01 cutoff. A time-forward test requires a separately governed source,
for example post-2024 RCSB coordinates with a locally pinned interaction
teacher and an explicit teacher-concordance audit.

## Required Phase 2A before any new model

Phase 2A is an audit-only stage. It must be preregistered and committed before
execution and must not train a new neural architecture.

1. Decompose sealed B5 pair logits under the actual mask and weights:
   `G = Projection_W(additive marginals) + C`.
2. Report pair AP for full `G`, additive projection and coupling residual `C`.
3. Compare `C` with degree-preserving rewiring, wrong-protein and wrong-ligand
   controls. Rewiring remains an evaluation null, never a training negative.
4. Census exact proteins with multiple ligands and scaffolds, and quantify how
   much true residue labels vary across ligands for the same protein.
5. Compare a protein-only residue oracle with a protein-plus-ligand residue
   oracle using component-level inference.

Decision rule:

```text
teacher ligand-conditioned increment weak
  -> current corpus cannot identify ligand-specific residue selection; stop.

teacher increment strong but B5 coupling weak
  -> authorize one preregistered residue-first typed mechanism field.

B5 coupling already survives all attribution controls
  -> skip architectural repair and move to sealed structural confirmation.
```

No Phase 2A outcome automatically authorizes affinity labels, DAVIS, few-shot
adaptation, production `z`, or changes to CSMO/Band/the frozen law operator.

## Immutable evidence index

| Artifact | SHA-256 |
|---|---|
| `P1_B5_GATE.json` | `99083ba91657cbfd6ad3c9124290e41907f706269259a5b4dc3c2075ffbdb22f` |
| `P1_MARGINAL_DECOMPOSITION.json` | `95527d46e7ca754e303acd218d7b803dd7dac33e28747d1ff0484aa33de15222` |
| `P1_B5_REPORT.md` | `19c9c2052c56dcd18df61c538bf274f47e40813483b356d433c52bb38984a8aa` |
| `P0_SEALED_PREDICTION_MANIFEST.json` | `db02961bfd2603713f1d1b68849ac754880ec6816e4352dfb01ded34ba6aaeed` |
| `PUBLICATION_TIME_CLOSURE_AUDIT.json` | `ce12e081e08b0f8acdae605b774764c1e840e6a10fb8b34d1df40495cb0d964c` |

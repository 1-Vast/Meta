# Phase 2A literature-grounded research plan

Date: 2026-08-10.

This is a research plan, not a preregistration or execution authorization.

> **EXECUTED 2026-08-10.** The plan below was superseded by the formal
> registration `research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2A.md`
> (SHA-256 `4e01401d…`) and carried out. Terminal verdict:
> **`LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING`**. Of the six
> candidate verdicts tabulated below, `TEACHER_GENERIC_POCKET_ONLY` was
> **refuted** — the labels are ligand-conditioned at the residue level. See
> `PHASE2A_SYNTHESIS.md`.

## Problem

Frozen ESM2 B5 passed all six development Gates for exact-residue localisation,
but a wrong ligand retained about 92.5% of residue AP. The current pair score may
therefore be explained by a generic protein pocket marginal plus ligand atom
propensity. Exact ligand-conditioned residue-atom coupling is not identified.

This ambiguity has two distinct causes that must be separated before training:

1. MONN labels may contain little within-protein, across-ligand variation;
2. the labels may contain coupling that B5 fails to recover.

## Literature boundary

Atom-residue pair maps and joint interaction/affinity training already exist in
[MONN](https://www.sciencedirect.com/science/article/pii/S2405471220300818).
Self-supervised protein representations in that family already exist in
[SPE-MONN](https://pmc.ncbi.nlm.nih.gov/articles/PMC9756525/) and
[CPAC](https://pmc.ncbi.nlm.nih.gov/articles/PMC9486597/). ESM2-based
ligand-conditioned site models are represented by
[LaMPSite](https://openreview.net/forum?id=MASqXhMTy7) and
[ProMoSite](https://pmc.ncbi.nlm.nih.gov/articles/PMC13228254/). Thus another
PLM, attention block, or generic pair head is not a defensible primary repair.

Strict two-entity separation remains necessary because ligand memorization is a
documented failure mode in affinity models
([Volkov et al.](https://www.nature.com/articles/s42256-023-00756-9)).
[DataSAIL](https://www.nature.com/articles/s41467-025-58606-8) and
[CleanSplit](https://www.nature.com/articles/s42256-025-01124-5) support keeping
protein and ligand similarity controls explicit rather than merging them into a
single unusable giant component.

## Phase 2A audit sequence

1. **Data census.** Count exact protein/construct groups with at least two
   scaffold-distinct ligands and report independent component depth.
2. **Teacher conditionality.** Measure within-protein residue-mask variation,
   wrong-ligand retention, Jaccard/PR overlap, and between-ligand variance.
3. **Weighted marginal projection.** On the actual observation mask, fit the
   additive residue/atom null and compute its orthogonal coupling residual.
   Ordinary double-centering is invalid for sparse or unequally weighted pairs.
4. **Coupling attribution.** Compare empirical coupling with wrong partners and
   multiple degree-preserving bipartite rewires. Rewiring is an evaluation null,
   not a biological non-binder or training negative.
5. **Label semantics.** Audit binary-label disagreement and potential
   positive-unlabelled contamination before considering a soft or PU teacher.

## Mutually exclusive verdicts

| Verdict | Only permitted next action |
|---|---|
| `PHASE2A_DATA_NOT_IDENTIFIABLE` | acquire a new governed structural supervision corpus; do not train |
| `TEACHER_GENERIC_POCKET_ONLY` | retain B5 only as generic-pocket evidence and close the exact-coupling claim on MONN |
| `LABEL_SEMANTICS_AMBIGUOUS` | rebuild a continuous/soft audited teacher, then repeat Phase 2A |
| `LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING` | register one ligand-conditioned residue residual head |
| `EDGE_COUPLING_PRESENT_B5_ABSENT` | register one marginal-orthogonal pair coupling head |
| `EDGE_COUPLING_ALREADY_IDENTIFIED` | do not repair; seek one sealed independent structural confirmation |

Only one trainable repair may be authorized by the audit. No verdict authorizes
simultaneous PLM, geometry, typed-interaction, attention, or affinity changes.

## Interface to MetaSieve mathematics

B5 and any residual head are measurement frontends, not `z`. A structural
candidate must later pass a source-affinity Gate against ligand-only and wrong
protein controls. Only then may the support design identify a low-dimensional
section: usable adaptation directions are limited to the support row space,
with rank, conditioning, query coverage, and abstention reported.

The frozen law operator remains unchanged:

```text
admitted biological statistic
  -> support-identifiable bounded summary z
  -> F(z)
  -> B(z)F(z)
  -> K(beta)
```

The bioinformatics contribution is therefore the marginal-audited measurement
and admission interface. The mathematical contribution remains the
identifiability-constrained statistic-to-law construction.

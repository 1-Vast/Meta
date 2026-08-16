# Protein-ligand interaction identifiability decision

Date: 2026-07-31

> **Superseded by** [`innovation_gate_decision_2026-07-31.md`](innovation_gate_decision_2026-07-31.md) and the strict v2 audit JSONs. The original rectangle totals below included same-homology pairs and raw pKd rows removed by the active registry policy; use the v2 artifacts for decisions.

## Result

The TRAIN-only audit used raw ChEMBL-37 activity rows restricted to the 559
registered TRAIN targets. pKi and pKd were kept separate. A rectangle is an
exact canonical ligand pair shared by two targets in one source document; a
rectangle was never treated as an independent biological observation.

| Quantity | All | pKi | pKd |
| --- | ---: | ---: | ---: |
| document-local exact-ligand rectangles | 17,120,772 | 15,346,780 | 1,773,992 |
| target-pair/document units | 39,166 | 12,785 | 26,381 |
| homology-pair bootstrap units | 7,996 | 3,588 | 7,115 |
| unit median absolute delta | 0.367 | 0.200 | 0.514 |
| unit order-reversal fraction | 0.385 | 0.358 | 0.398 |
| strict same-assay rectangles | 0 | 0 | 0 |

The apparent rectangle count is therefore not a strict same-assay crossed
experiment. The pKi unit difference-in-differences standard deviation is 1.34,
while the propagated replicate-noise scale is 1.47 at the 90th percentile and
2.26 at the 95th percentile. The median-noise ratio is not an adequate decision
statistic because replicate noise is heavy-tailed.

## Decision

`AUDIT_RECTANGLES_EXIST_STRICT_ASSAY_COMPARABILITY_LIMITED`

This is a **no-go for protein-conditioned architecture expansion on the current
affinity registry**. It is not a claim that protein-ligand interaction is absent
in biology. It says the current observations do not provide a clean, assay-
comparable interaction estimand that can justify adding TADAM/AdaMBind task
conditioning, a deeper cross-attention block, a graph ligand encoder, or an
auxiliary protein loss.

The decision is consistent with the existing AnchorDelta result: joint training
improved aggregate ranking, but correct and wrong proteins were indistinguishable
under component bootstrap. AdaMBind's easy-to-hard scheduler cannot repair this
because it reallocates meta-training tasks; it does not add crossed assays or
target-specific labels.

The existing AnchorDelta path was also smoke-tested in the `drug` CUDA
environment after the audit (`8` gate targets, `6` independent components).
The smoke produced pairwise accuracy `0.4031` for AnchorDelta versus `0.4023`
under wrong-protein replacement, consistent with the established failure mode.

## Reopening condition

Reopen only after a new data protocol supplies source/assay-comparable crossed
measurements (or an independently validated interaction label). Then run the
low-capacity residual probe from the supplied audit: protein-free antisymmetric
ligand-pair baseline, correct protein, homology-matched wrong protein, random
wrong protein, and an oracle target-ID ladder under nested target/homology-cold
splits. A model change is admissible only if the correct-protein arm beats both
protein-free and matched-wrong-protein controls with component-level confidence
intervals above zero.

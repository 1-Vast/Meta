# Final completion audit — 2026-07-26

## Terminal claim

Category **② `SIGNAL_PRESENT_EVIDENCE_INSUFFICIENT`** is proven. Category ① is not authorized.
Category ③ is too pessimistic because the strict A1 probe contains a statistically positive
protein-dependent effect; it simply fails the substantive-effect and coarse-taxonomy controls.

## Requirement-to-evidence audit

| Requirement | Authoritative evidence | Audit status |
| --- | --- | --- |
| Audit prior failures; do not rename/capacity-scale them | `history.md`; frozen negative routes and reopening conditions | Proven |
| At most three autonomous candidates; attachment ideas excluded | `kirhub_autonomous_preregistration.md`; A1/A2/A3; U1–U3 separately ledgered | Proven |
| At most two candidates enter experiment | Only A1 entered; A2 gate-blocked, A3 data-blocked | Proven |
| Pre-result hypothesis, data flow, controls and gates | Autonomous preregistration plus pre-strict-run firewall correction | Proven |
| Literature support for standard machinery | DOI registry in firewall correction: ESM-2, ECFP, Murcko, Tanimoto, k-mer comparison, kernel regression | Proven |
| Public/traceable data and contamination statement | KirHub workbook/hash/license record, KLIFS snapshot, Wikidata/ChEMBL structure sources, frozen ESM provenance | Proven |
| Target identity and homology isolation | 324 full-sequence 4-mer-containment connected components; component-only folds | Proven |
| Ligand identity, scaffold and high-similarity isolation | 79 connected components joining parent identity, Murcko equality or Morgan Tanimoto >=0.50; maximum between components 0.478873 | Proven |
| Assay/document/publication isolation | Impossible inside the single KirHub release; explicitly excluded from the scope of positive inference | Limitation, not claimed |
| Low-cost mechanism test before learned model | H0, then nonparametric A1; no neural/learned innovation run | Proven |
| Destructive/matched controls | ligand-only, protein shuffle, matched random protein and KLIFS-group centroid | Proven |
| Independent-unit statistics and power | 308 homology-component units; grouped bootstrap; MDE80 +0.0160 at SD 0.10 | Proven |
| Stop when mechanism gate fails | +0.0290 < +0.030 and true ESM does not beat group centroid; A2 not run | Proven |
| No confirmation/sealed-label reuse | All authoritative JSON flags false; no multi-seed/confirmation execution | Proven |
| Reproducibility | 38 relevant tests pass; strict component registry rebuilds exactly; strict JSON, compile and CUDA checks pass | Proven |
| Required final reporting | `task.md`, `history.md` and strict decision report contain candidate source, experiments, statistics, dependencies, excluded explanations, root cause and reopening conditions | Proven |

## Final scientific interpretation

Under full-sequence homology and identity/scaffold/Tanimoto chemical isolation, frozen ESM improves
component-macro ranking over ligand-only by +0.0290 [+0.0083,+0.0497] and beats shuffled/random
protein controls. The effect is real but one thousandth below the preregistered +0.030 minimum and
does not beat the KLIFS-group centroid (-0.0110 [-0.0318,+0.0099]). It is therefore compatible with
coarse kinase taxonomy rather than proven sequence-resolved transferable SAR.

The most credible bottleneck is not posterior mathematics or insufficient model capacity. It is
the lack of independent factorial measurements spanning target homology, chemical components and
assay/document/publication sources, together with replicate-level mutation noise. Reopening
requires such data; relaxing thresholds, adding a Transformer, increasing seeds or training longer
is inadmissible.


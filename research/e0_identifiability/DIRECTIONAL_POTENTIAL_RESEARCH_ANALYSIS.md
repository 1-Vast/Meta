# Directional Potential Proposal: Research Analysis

Updated: 2026-08-07

## Executive Verdict

The proposal points at the correct unresolved boundary but combines stages that
are not jointly authorized:

```text
synthetic objective closure
-> structural directional/type identifiability
-> reference-state structural log-odds
-> source OOF delta-affinity
-> few-shot mechanism adaptation
-> biological statistic admission to z
```

Only the first stage was executable under the current label-free research
contract. E0R2 passed it. This proves that the frozen 240-dimensional statistic,
corrected residual/difference objective, train design and deterministic solver
are mutually consistent. It does not prove that real affinity direction or a
directional interaction potential is identified.

## Claim Audit

| Proposal claim | Verdict | Evidence or missing test |
|---|---|---|
| A signed geometry-to-affinity layer is the next useful hypothesis | Partially supported | Historical global readouts failed, but this is not proved to be the unique missing layer. |
| P1B has already lost affinity direction | Unsupported | Real affinity direction has not been tested; the synthetic teacher is recoverable from frozen P1B geometry and chemistry. |
| Type and orientation are required | Biochemically plausible | PLIP and directional molecular models show why angles/types can matter, but a held-out structural Gate is required here. |
| `-log(p_bound/p_ref)` is an affinity energy | Incorrect | It is a reference-dependent structural log-density ratio. Its offset, scale and sign depend on the reference. |
| Distance-distribution expectation preserves uncertainty | Overstated | `E[u(d)]` propagates one expectation, not the full distribution or an uncertainty certificate. |
| OOF ligand residual is required | Supported as protocol | It prevents in-sample ligand-prior leakage, but the residual still contains assay noise and model misspecification. |
| Cross-fitting makes the design R-learner orthogonal | Not established | The analogy is useful; no Neyman-orthogonal score or causal estimand has been defined. |
| Point plus residual-difference loss is semantically aligned | Supported | E0R1 confirmed the old residual/total mismatch. Pair differences do not create new labels and require frozen task weighting. |
| `d <= k` makes the few-shot adapter identifiable | False as stated | Rank, conditioning and query row-space coverage are also necessary. Ridge fills missing directions with a prior. |
| Adapter coordinates are named mechanisms | Not without anchoring | `U a` has a rotation gauge. Named H-bond or ionic coordinates require a fixed biological basis. |
| Adapter inverse Hessian is uncertainty | Conditional | It needs an explicit noise model, prior and parameter set; otherwise it is not the theory's outer certificate. |
| The mechanism vector can be inserted into `z` unchanged | Engineering-compatible only | It must first be observable, bounded, fixed-dimensional, query-label-free and gauge-invariant, then pass affinity and transfer Gates. |

## E0R2 Result

| Metric | Result | Gate |
|---|---:|---:|
| train RMSE | `3.188e-8` | `<=1e-6` |
| train maximum absolute error | `2.007e-7` | descriptive |
| corrected objective | `1.567e-15` | `<=1e-12` |
| full-gradient L2 | `6.150e-17` | `<=1e-8` |
| correct CI | `0.99737` | `>=0.80` |
| correct-minus-ligand | `+0.51283` | `>=0.10` |
| correct-minus-deranged | `+0.38289` | `>=0.10` |

Verdict: `SYNTHETIC_OBJECTIVE_DESIGN_SOLVER_IDENTIFIED`.

The eight holdout tasks are a reused development diagnostic. The derangement
control also inherits concentrated wrong-protein reuse. These CI values localize
the synthetic numerical failure; they do not estimate biological transfer.

## Next Registered Gates

### T-DIR: Structural Directionality

Compare contact/distance, chemistry, frozen local states and an oracle
coordinate ceiling on homology/scaffold-held-out structures. Preserve explicit
residue identity before comparing it with the 128-slot bank. Each interaction
channel must exceed prevalence and deranged controls with a component-bootstrap
positive lower bound. A PASS establishes only structural type identifiability.

### T-REF: Reference-State Log-Odds

Freeze one chemistry-marginal-matched native/reference construction, binning,
pseudocount and sign convention before scoring. Wrong proteins cannot be treated
as non-binders. A PASS may be called `STRUCTURAL_LOG_ODDS_IDENTIFIED`, never
affinity energy identified.

### E-AFF: Source OOF Delta-Affinity

Only after structural channels pass, use the governed ChEMBL closure folds to
construct a cross-fitted ligand prior and train the frozen mechanism feature on
point residual plus within-task residual differences. The source Gate remains:

```text
correct-minus-ligand >= 0.03
correct-minus-deranged >= 0.03
component-bootstrap 95% LCB > 0 for both
```

### M-FS: Identified Few-Shot Section

Only after the population mechanism passes E-AFF. Report support design rank,
condition, query row-space residual and an explicitly justified outer radius.
Abstain when the query mechanism direction is not identified by support.

## Theory Boundary

The frozen theory accepts an observable compact statistic `z(S,Q,gamma)` and
then constrains the law-valued output through

```text
z -> F(z) in simplex -> B(z)F(z) -> K(beta).
```

It does not prove that a proposed interaction potential contains affinity
information, and it does not supply a pairwise-ranking theorem. Deep integration
requires biology to define admitted, bounded and intervention-sensitive
coordinates of `z`; it cannot be claimed merely because a scalar potential can
be concatenated to the operator input.

## Primary Sources Checked

- ITScoreAff: https://pubmed.ncbi.nlm.nih.gov/37822259/
- PLIP: https://pmc.ncbi.nlm.nih.gov/articles/PMC4489249/
- DimeNet: https://arxiv.org/abs/2003.03123
- Noise-contrastive estimation: https://proceedings.mlr.press/v9/gutmann10a
- NequIP: https://www.nature.com/articles/s41467-022-29939-5
- MACE: https://papers.neurips.cc/paper_files/paper/2022/hash/4a36c3c51af11ed9f34615b81edb5bbc-Abstract-Conference.html
- PLANET v2.0 preprint: https://arxiv.org/abs/2601.07415


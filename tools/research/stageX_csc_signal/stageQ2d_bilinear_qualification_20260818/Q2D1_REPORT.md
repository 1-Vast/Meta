# Stage Q2d-1 report — isomorphic bilinear positive control (final)

Frozen prereg SHA: 4f7e80027b9b82564bf1ea262d360813f23ffe77f8378ca54b490cc6752fa3d1.
Artifacts: Q2D1_BILINEAR.json (ladder), Q2D1_CLOSED_FORM_DIAGNOSTIC.json
(information bound), Q2D_LITERATURE_CORRECTIONS.md.

## Ladder results (medians over truth seeds 0-2)

| phase | exact_bilinear dz / sp | oracle_latent dz / sp | exact gate |
|---|---|---|---|
| A z-scale identity | 0.493 / 0.072 | 0.662 / 0.417 | FAIL |
| B + sigmoid %-scale loss | 0.508 / 0.058 | 0.674 / 0.425 | FAIL |
| C + 70% missingness | 0.509 / 0.129 | 0.687 / 0.431 | FAIL |
| D + interval censoring | 0.487 / 0.090 | 0.672 / 0.438 | FAIL |
| E + main-effect competition | 0.467 / -0.006 | 0.699 / 0.461 | FAIL |

Negative arms failed the gate in every phase (additive_only dz 0.000,
no_interaction_head dz 0.000, ligand_only/shuffled/random 0.44-0.53).

## Closed-form information bound (diagnostic only; deployment stays
## gradient-trained)

- in-sample (train holdout): dz 0.951 / sp 0.927 — the planted interaction is
  fully recoverable where rows AND columns are observed;
- half-cold (train parents x val ligands): dz 0.500 / sp 0.026 — recovery
  collapses to chance the moment the ligand is unseen.

## Attribution

1. The learner and optimization are NOT the failure: the closed-form bound
   reaches 0.95 in-sample, and the endpoint ladder (A->E) barely moves the
   exact-bilinear result (0.493 -> 0.467). Sigmoid, missingness, censoring
   and main-effect competition are each ruled out as the killer step.
2. The truth generator is structurally inconsistent with the mechanism
   under test. I(p,l) = (P U)(L V)^T with U, V iid random matrices makes the
   unseen-ligand factors (L V) unidentifiable from ligand features: learning
   V (2048 x 4) plus A (1700 x 4) needs ~15k parameters against 6,429 train
   observations, so cold-ligand extrapolation is information-theoretically
   impossible for any learner. The in-sample 0.95 and half-cold 0.50 are the
   two faces of the same fact.
3. Consequence: the planted harness as generated tests an ID-space
   interaction, not the transferable 'pocket-residue physicochemical x ligand
   substructure interaction field' that the mechanism claim requires. Q2d-2
   is NOT STARTED (frozen decision rule: Q2d-1 must pass first).

## Next preregistration (Q2d-1b)

Regenerate the planted truth with FEATURE-CONDITIONED low-dimensional
factors: U = f(pocket z-scales / per-position physicochemical descriptors,
~10-85 dims), V = g(ligand substructure descriptors, ~100-200 dims), rank 4,
double-centred, standardized to tau*=1.0, then rerun the same ladder
(A-E), the same arms and the same 0.30/0.70/0.05 gate. This makes the
factor maps identifiable at the given sample size and directly tests the
bounded mechanism claim. The existing frozen gates are not modified.

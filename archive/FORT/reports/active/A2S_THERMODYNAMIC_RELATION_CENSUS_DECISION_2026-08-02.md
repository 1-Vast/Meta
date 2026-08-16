# A2S Thermodynamic-Relation Census Decision

Date: 2026-08-02  
Branch executed: `research/a2s-thermodynamic-cycle-20260802`  
Artifact: `reports/active/a2s_thermodynamic_relation_census_2026-08-02.json`  
Runner: `research/a2s_thermodynamic_relation_census.py`  
Tests: `tests/test_a2s_thermodynamic_relation_census.py`

**FACT - decision:** `EXPLICIT_TRANSFORMATION_GRAMMAR_NOT_ADMITTED`.

## Firewall And Mechanical Result

**FACT.** The run used only source `fit` and source `probe` labels. It did not
request `locked` or recipient labels. The unit suite returned `6 passed`.

**FACT.** The census covered 52,677 measurements and 48,594 unique molecules.
RDKit fragmentation produced 342,879 fragment records and 278,096 matched-pair
edges: 190,195 in `fit` and 87,901 in `probe`.

## Transformation Transfer

**FACT.** The `fit` role contained 98,026 distinct transformations. Only 4,427
were repeated across at least three independent fit targets/components and
therefore qualified as robust source priors.

**FACT.** Those priors matched 17,070 probe edges across 89 probe targets and 82
homology components. Every matched edge was chemically local. The held-probe
low-similarity stratum, defined by Tanimoto below 0.35, contained zero robust
fit-known transformations, zero targets, and zero components.

**FACT.** On all 17,070 matched probe edges, the fixed source transformation
prior did not improve held-component direction:

| Contrast | Mean gain | Component-bootstrap lower 95% | Components |
|---|---:|---:|---:|
| pairwise proper loss | +0.0097 | -0.0059 | 82 |
| transformation direction | -0.0181 | -0.0644 | 77 |

**INFERENCE.** A repeated MMP rule is not a stable cross-target response rule in
this corpus. The positive point estimate for proper loss is unresolved, while
direction is negative and both intervals include zero.

## Passive Few-Shot Coverage

**FACT.** In all-query cells, fit-known transformation coverage rose from 9.1%
at k=1 to 20.9% at k=3 and 28.8% at k=5. This is a local-chemistry opportunity
and remains a baseline stratum.

**FACT.** In the required low-similarity stratum, fit-known coverage was exactly
zero at k=1, k=3, and k=5. MMP connectivity of any kind was 0.047%, 0.096%, and
0.161%, respectively. No low-similarity episode had at least eight fit-known
queries, so the registered floor of 47 powered components was missed by 47.

**INFERENCE.** An explicit finite SAR grammar cannot be the main passive A2S-DTA
mechanism. It has no deployment path to the chemically distant queries where
TRACE/KRR fail and where a successor is required to add value.

## Thermodynamic-Cycle Audit

**FACT.** For any scalar score function `f_t(x)`, defining
`Delta_t(a,b)=f_t(b)-f_t(a)` gives antisymmetry and exact cycle closure
algebraically. These constraints do not by themselves reduce the dimension of
the target-specific score function or add support information.

**FACT.** The provisional TCRS predictor
`mu_t(x) + z_t^T phi_t(x)` is a low-dimensional target code over learned basis
functions. After integration of the pair potentials it belongs to the same
functional family as the already-tested BIR/IDA global-code arm.

**INFERENCE.** Calling the basis functions perturbation potentials does not make
TCRS scientifically distinct. A distinct successor would need either observed
cross-target relation coverage, which this gate rejects, or a representation in
which the basis is jointly protein-ligand conditioned and demonstrably changes
held-target k-shot recoverability.

## Consequences

1. **FACT.** Stop the explicit MMP grammar and the provisional TCRS formulation
   as primary mechanisms. Do not increase vocabulary size, network capacity, or
   training epochs to rescue this gate.
2. **FACT.** Retain MMP, Matsy-style series transfer, ActFound-style pairwise
   prediction, scaled KRR, and TRACE as chemically local baselines.
3. **FACT.** Retain assay identity as a stratification/control variable, not an
   adaptation state.
4. **HYPOTHESIS.** The remaining representation-shaped question must be tested
   with target-conditioned protein-ligand interaction coordinates, not another
   target-independent ligand basis.
5. **FACT.** No code is promoted to `model/` or `script/`; no major breakthrough
   has been achieved.


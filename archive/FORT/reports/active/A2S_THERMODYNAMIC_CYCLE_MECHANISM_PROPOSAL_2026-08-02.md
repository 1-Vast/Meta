# A2S Thermodynamic-Cycle Meta-Adaptation Proposal

Date: 2026-08-02
Status: pre-implementation hypothesis and source-only gate specification
Active branch: `research/a2s-thermodynamic-cycle-20260802`

## Exploration Objective

Learn a transferable, support-conditioned target adaptation state from abundant
source targets that uses only `k in {1,3,5}` recipient measurements to make a
small number of justified, query-dependent ranking interventions for a strictly
unseen target, while preserving the support-free DTA prediction when evidence
is absent or invalid.

The exploration is not allowed to turn into calibration, retrieval, active
learning, uncertainty estimation, kernel weighting, or unrestricted function
generation. Those remain baselines or auxiliary controls.

## Sequential Failure Diagnosis

**FACT.** TRACE found that support-residual transport is effectively null below
nearest-support Tanimoto about 0.35 and useful mainly for local chemistry.

**FACT.** MODE found target-specific ranking headroom in the distant stratum,
but naive discrete response modes were not identifiable at k=3/5. Its synthetic
mode control also failed, so the experiment did not distinguish representation
failure from absence of target states.

**FACT.** RIP found a large hindsight ceiling for sparse ranking interventions,
but observed reliability/margin variables could not identify useful actions.

**FACT.** The hotspot-sparsity falsification showed that apparent coordinate
sparsity was not stable under rotations or held-out components. The named basis
retained less top-8 mass than arbitrary rotations, and effective state size was
about eight rather than two or three.

**FACT.** Exact-assay grouping reduced residual SD from 1.2647 to 0.8275 pKi,
but the assay-coherent EB gate failed at k=3. Only one underpowered k=5 distant
cell passed, while the synthetic control passed.

**INFERENCE.** The evidence rejects five adaptation objects as primary claims:
local support transport, naive target modes, reliability-only selection,
coordinate-sparse absolute heads, and assay-conditioned global heads.

**INFERENCE.** TRACE is insufficient as the final contribution because its
learned object is pair/query reliability for transporting measured residuals.
It does not learn a compact target state that changes the response to chemical
perturbations outside the support neighborhood.

## Successor Biological Principle

**FACT.** Binding affinity is a free-energy state function. Relative effects
within a coherent target/assay satisfy antisymmetry and thermodynamic cycle
closure up to measurement error:

\[
\Delta_t(a,b) = y_t(b)-y_t(a),\qquad
\Delta_t(a,b)=-\Delta_t(b,a),
\]

\[
\Delta_t(a,b)+\Delta_t(b,c)+\Delta_t(c,a)=0.
\]

**FACT.** Matched molecular pair/series studies show that some substituent
orders transfer across medicinal chemistry series. Activity-cliff and
nonadditivity studies show that these effects are conditional: strong
nonadditivity often accompanies binding-mode changes, and leakage-free cliff
prediction frequently depends on one or a few assay-neighbor relationships.

**HYPOTHESIS.** The transferable object is not an absolute target SAR head. It
is a small response state describing how an unseen target modifies a shared
repertoire of chemical perturbation potentials. Source tasks teach the
perturbation potentials; k support measurements identify at most one state
coordinate at k=1 and at most two centered coordinates at k=3/5.

This principle is generalizable in form because relative free-energy changes,
antisymmetry, and cycle closure do not depend on a particular target family or
fingerprint. Whether the current public corpus contains enough coherent
cross-target perturbation evidence is an empirical question and is gated below.

## Literature And Novelty Boundary

Primary precedents that constrain the claim:

- MMP-cliffs replace opaque fingerprint similarity with defined substructure
  changes: https://doi.org/10.1021/ci3001138
- Matsy transfers preferred orders in matched molecular series:
  https://doi.org/10.1021/jm500022q
- DeltaDelta learns congeneric potency ranking:
  https://doi.org/10.1039/c9sc04606b
- DeepDelta directly learns molecular property differences and can scaffold
  hop: https://doi.org/10.1186/s13321-023-00769-x
- PBCNet learns relative affinities for congeneric ligands:
  https://doi.org/10.1038/s43588-023-00529-9
- ActFound combines pairwise bioactivity learning and meta-learning:
  https://doi.org/10.1038/s42256-024-00876-w
- FS-CAP predicts a new assay from a support-set embedding:
  https://doi.org/10.1021/acs.jcim.4c00485
- Activity-cliff transfer can improve DTI representations:
  https://doi.org/10.1021/acs.jcim.5c00484
- Strong nonadditivity is associated with binding-mode changes after accounting
  for experimental error: https://doi.org/10.1021/acs.jcim.5b00018
- Heterogeneous public Ki measurements have an estimated SD of 0.54 pKi:
  https://doi.org/10.1021/jm300131x
- Assay context improves proteochemometric models but does not remove all
  heterogeneity: https://doi.org/10.1021/acs.jcim.5c00603
- Leakage-controlled activity-cliff prediction from bioactivity profiles is
  often driven by nearest assay relations:
  https://doi.org/10.1186/s13321-026-01210-9

**FACT FROM RELEASED CODE.** ActFound adapts a high-dimensional final layer by
MAML-style inner updates, but its released `system_base.py` skips that update
when support size is `<=5`. At the A2S budget it still predicts with support
anchors, learned pairwise differences, and similarity weighting. FS-CAP uses an
averaged support-set context and a generic decoder.

**NOVELTY BOUNDARY.** Pairwise learning, few-shot assay prediction, MMP rules,
cycle losses, and meta-learning are not individually novel. A defensible
increment would be the combination of:

1. an explicitly 1-2 degree-of-freedom target response state at k<=5;
2. a learned support-to-state operator that is active at this budget;
3. source-learned perturbation potentials constrained by antisymmetry/cycles;
4. a hard STOP state and bounded intervention budget;
5. target-disjoint, homology-component evaluation showing that this state is
   load-bearing beyond ActFound/FS-CAP/KRR and fixed analytic solvers.

## Candidate Mechanisms

### C1. Thermodynamic-Cycle Response State (TCRS)

**Learned object.** A bounded 1-2 dimensional target state controlling a small
set of source-learned perturbation potentials.

**Transfer rationale.** Source targets repeatedly expose how defined chemical
changes alter binding under different interaction regimes. Antisymmetry and
cycle closure share structure across all targets.

**k-shot argument.** k labels contain only k-1 independent centered contrasts.
Use one state coordinate at k=1 (sign/STOP only) and no more than two at k=3/5.

**Difference from prior art.** Unlike MAML/ANIL/ActFound, no high-dimensional
weight vector is adapted. Unlike FS-CAP/MetaDTA, no unrestricted set context is
decoded. Unlike KRR, labels infer a target response state rather than weights
over support labels.

**Verdict.** Highest-potential hypothesis, conditional on the source-only gate.

### C2. Finite Assay-Disentangled SAR Grammar

**Learned object.** A finite set of auditable MMP clauses with activation,
composition, reversal, and unsupported states.

**Transfer rationale.** Matched molecular series exhibit preferred substituent
orders across programs.

**k-shot argument.** One or two clauses could be selected from k=3/5 contrasts.

**Risk.** Passive random support may have negligible MMP coverage, and global
rules fail under binding-mode nonadditivity. This route is admitted only if the
MMP connectivity census is strong.

### C3. Discrete STOP/SWAP Rank-Edit Policy

**Learned object.** A policy selecting STOP or at most B adjacent ranking swaps.

**Transfer rationale.** Source tasks teach which residual patterns justify a
small ranking correction.

**k-shot argument.** The action space is finite and B can be <=k-1.

**Risk.** It is candidate-list dependent and is closer to support-conditioned
learning-to-rank than inductive affinity adaptation. Distractor, subset, and
library-size invariance are hard requirements.

### C4. Pocket-Conditioned Perturbation Prior

**Learned object.** A protein/pocket prior over perturbation-response states,
with k-shot evidence updating one or two logits.

**Transfer rationale.** Conserved interaction microenvironments can preserve
responses to functional-group changes across homologous targets.

**k-shot argument.** Support only chooses among a small prior repertoire.

**Risk.** Previous pooled ESM and residue-contact student routes failed held
domain semantic controls. This cannot be primary without a new pocket-level
positive control and protein shuffle/zero sensitivity.

### C5. Assay-State Transition Adapter

**Learned object.** A small construct/readout state that changes query response
potentials.

**Transfer rationale.** Measurement context changes the observed thermodynamic
state and reduced residual noise in the exact-assay audit.

**k-shot argument.** Discrete context states are low-cardinality.

**Risk.** The exact-assay gate failed at k=3 and assay metadata may be absent at
deployment. Keep as a conditioning variable, not the primary state.

### C6. Cycle-Constrained Learned Optimizer

**Learned object.** A learned update rule for a pairwise scoring head, constrained
to preserve antisymmetry and cycle closure.

**Transfer rationale.** Source episodes teach how gradients from pairwise
support evidence should change relative-affinity predictions.

**k-shot argument.** Only credible if the update is projected to <=2 effective
degrees of freedom.

**Risk.** Without that projection it is MAML/ANIL/ActFound with another loss and
violates the information budget. It is a baseline/ablation, not the selection.

## Selected Hypothesis: TCRS

Selection is provisional. It has the strongest biological basis and the
clearest distinction from similarity weighting, but no core model will be
implemented unless the relational source-only gate passes.

Let `mu_t(x)` be the frozen support-free prediction and
`r_i = y_i - mu_t(x_i)`. A learned relation encoder emits antisymmetric support
edge evidence:

\[
e_{ij}=E_\theta(p_t,x_i,x_j,r_j-r_i,m_i,m_j),\qquad e_{ji}=-e_{ij}.
\]

The adaptation operator emits a bounded state and a hard action bit:

\[
(a_t,z_t)=A_\psi(\{e_{ij}\}_{i<j}),\qquad
a_t\in\{0,1\},\quad z_t\in[-1,1]^d,
\]

with `d=1` at k=1 and `d<=2` at k=3/5. Null, removed, or unsupported evidence
forces `a_t=0` and `z_t=0` exactly.

Each learned response channel is a ligand potential `phi_j`; differences are
therefore antisymmetric and cycle-consistent by construction:

\[
d_j(a,b)=\phi_j(p_t,x_b)-\phi_j(p_t,x_a).
\]

The adapted score is

\[
\hat y_q=\mu_t(x_q)+\hat b_t+
a_t c\tanh\left(\frac{\sum_{j=1}^d z_{t,j}\phi_j(p_t,x_q)}{c}\right).
\]

`b_t` is a separately reported rank-null calibration channel. The bounded rank
term is the claimed intervention. At k=1 it may use only a one-bit polarity or
STOP state; any stronger one-shot state is prohibited.

## Trainable Path And Objective

Trainable modules:

1. relation/perturbation encoder `E_theta`;
2. response potentials `phi_theta`;
3. support-to-state operator `A_psi`;
4. hard STOP and bounded-action head.

Source residuals must be produced out of fold by the complete frozen base.
Training uses only fit components; model selection uses held-out source
components. A minimal objective is

\[
L = L_{pairwise\ proper}+\lambda_q L_{query\ rank}
  +\lambda_{cycle} L_{cycle}+\lambda_{state} L_{budget}
  +\lambda_{cf} L_{counterfactual\ no-op}.
\]

`L_cycle` tests triangle closure in coherent target/assay groups. `L_budget`
limits effective state dimension and intervention magnitude. Counterfactual
arms use residual permutation at k>=3, norm-matched sign/transplant controls at
k=1, wrong-target support, support removal, and random support.

## Source-Only Admission Gate Before Core Implementation

### Gate R0: relation and MMP coverage

- Build transformation vocabulary from fit compounds only.
- Measure exact-MMP and learned-relation coverage in fit/probe separately.
- Report passive k=1/3/5 support-to-query coverage by similarity stratum.
- Stop the explicit grammar route if fewer than 47 held-out components support
  the preregistered low-similarity contrast or if coverage is concentrated in a
  single target family.

### Gate R1: transferable relation signal

- Fit relation representation on fit components only.
- Test whether it improves pairwise proper loss/CI on probe components,
  especially below Tanimoto 0.35.
- Compare with Tanimoto KRR, scaled KRR/TRACE, ActFound-style anchor-difference,
  FS-CAP-style context, descriptor difference, random relation features, and
  assay-only effects.

### Gate R2: state complexity and k-shot identifiability

- Abundant-target oracle states must have effective dimension <=2 on probe.
- The support-to-state operator at k=3 and k=5 must recover a stable state under
  support resampling and beat fixed ridge/KRR/sparse solvers on the same
  representation.
- k=1 may activate only one bit/coordinate and must pass sign-flip,
  norm-matched transplant, and wrong-target controls.
- A synthetic cycle-state positive control must pass before any real-data null
  is interpreted.

### Gate R3: load-bearing adaptation

At k=3 and k=5 in low-similarity queries, require all of:

1. positive component-bootstrap lower 95% bound for CI/NDCG gain over base;
2. gain over scaled KRR/TRACE and fixed analytic state estimation;
3. correct support over wrong target, residual permutation, random support,
   support removal, coordinate permutation, and frozen/random operator;
4. positive target-level confidence interval and reduced negative-transfer rate;
5. no support-free/source-head degradation.

For k=1, require non-degradation plus a positive effect only if the one-bit
state passes all available destructive controls. A k=1 assignment claim cannot
be made because within-support permutation is undefined.

## Mandatory Ablations

1. Remove `A_psi` and set `z=0`.
2. Replace `A_psi` with closed-form ridge/KRR and sparse fixed solvers.
3. Freeze or randomize `A_psi`.
4. Remove cycle/antisymmetry constraints.
5. Randomize or rotate response coordinates.
6. Wrong-target support.
7. Residual-label permutation.
8. Random support assignment.
9. Support removal/residual-null.
10. Protein shuffle and protein zero.
11. Query permutation, distractor insertion, subset, and library-size checks.
12. Low/high support-query Tanimoto and exact-assay/cross-assay strata.

## Stop Rule And Promotion Rule

Stop TCRS without increasing capacity if relation coverage is inadequate, the
fit-learned representation does not transfer to probe components, effective
state dimension exceeds two, k=3/5 states are unstable despite a passing
synthetic control, gains remain local/KRR-equivalent, or an analytic solver
matches the learned operator.

Only after all gates pass may core code be promoted to `model/` and utilities
to `script/` (or the repository's reconciled script path), tested, hashed, and
frozen. Promoted files must not be modified afterward.


# A2S Protein-Conditioned Interaction Response State Proposal

Date: 2026-08-02  
Status: source-only preregistration; no mechanism admitted  
Active branch: `research/a2s-interaction-state-20260802`

## Exploration Objective

**HYPOTHESIS.** Learn a transferable, support-conditioned adaptation state from
abundant source targets that uses only `k in {1,3,5}` measurements to apply a
bounded, query-dependent ranking intervention to a strictly unseen target.
Reliability, retrieval, calibration, and local interpolation remain auxiliary
modules or baselines, not the contribution.

**FACT.** The active objective remains:

> Learn a transferable few-shot meta-adaptation mechanism for A2S-DTA.

## Sequential Diagnosis

**FACT.** TRACE learns pair/query reliability for transporting observed support
residuals. Its effect is essentially null below support-query Tanimoto 0.35, so
it does not learn a target adaptation state outside local chemistry.

**FACT.** MODE found real within-target, scaffold-disjoint ranking headroom, but
source target heads had a nearly flat spectrum: the top three directions held
34.7% of variance and rank-2 source projection retained -6% of full-head gain.
The empirical-Bayes learning-curve knee was near k=10, not k<=5.

**FACT.** RIP found a useful hindsight intervention ceiling, but observable
margins and posterior uncertainty could not identify the useful actions.

**FACT.** The hotspot branch's apparent coordinate sparsity was stronger after
arbitrary rotations, was unstable across splits, and required about eight
coordinates. The biological-hotspot interpretation was withdrawn.

**FACT.** Exact-assay grouping reduced residual SD from 1.2647 to 0.8275 pKi but
did not admit a k=3/5 distant-query state. Explicit MMP priors had zero passive
low-similarity coverage and failed held-component direction transfer.

**FACT.** BIR meta-learned `U` by differentiating through an empirical-Bayes
solve, but `phi(x)=g_frozen(x)U` was a protein-free linear projection of frozen
ligand features and the support-free base was fixed. Its state scale was about
0.18-0.19 pKi against noise about 0.997 pKi; the certificate fired on 0% of
episodes.

**INFERENCE.** TCRS is not a new escape from BIR. Once its pair potential is
integrated, its score is `mu(x)+z^T phi(x)`, and antisymmetry/cycle closure are
automatic for every scalar potential. The MMP census supplies no distant-query
relation coverage that would distinguish it.

**INFERENCE.** One scientifically open distinction remains. BIR and MODE used
target-independent ligand directions. They did not test whether a source-only
learned coordinate `phi(p,x)`, conditioned jointly on local protein sequence
environments and ligand motifs, can point in different ligand directions for
different proteins while requiring only one or two recipient coefficients.

## Biological Principle And Generalizability

**FACT.** Binding specificity is compositional: ligand functional groups are
interpreted by local chemical environments in a protein, and the same ligand
motif can contribute differently in different environments. ESM-2 segment
features and 2D ligand features are locally available for every fit/probe target
and compound in the source substrate.

**HYPOTHESIS.** A support-free DTA model captures the average interaction law,
while its unseen-target ranking error can be corrected by recalibrating at most
two reusable protein-ligand interaction channels. The channels transfer; the
small recipient state contains only their deviations.

**INFERENCE.** This principle can explain why a global ligand subspace is flat:
`phi(p_1,x)` and `phi(p_2,x)` need not define the same direction over compounds.
It does not assume a shared target head, a binding hotspot, an observed pose, or
an MMP connection.

**FACT.** The principle is not yet established. Earlier pooled-ESM head
prediction and fixed KLIFS pocket-coordinate tests did not beat their strong
protein-free controls. Segment conditioning is therefore load-bearing and must
beat protein-zero, protein-shuffle, pooled-protein, and random-coordinate arms.

## Candidate Mechanisms

| Candidate | Learned object | k<=5 identifiability | Transfer argument | Decision |
|---|---|---|---|---|
| **C1. Protein-Conditioned Interaction Response State (PIRS)** | two target-conditioned interaction channels plus a learned support-to-state operator | no rank state at k=1; one coordinate at k=3; two at k=5 | local protein environments and ligand motifs are reusable while channel activation is target/query specific | **selected for source-only gate** |
| C2. Low-rank learned gradient program | one or two learned update directions | dimension is explicit | source tasks teach useful gradient directions | defer: biologically unstructured and closest to Meta-SGD/LEO |
| C3. Discrete response mode | one mode/STOP logit | finite state is nominally compatible | source targets share response classes | reject: MODE G2/A4 found no small shared mode structure |
| C4. Sparse determinant selector | two or three active coordinates | sparse recovery could fit k=5 | targets reuse different determinants | reject: rotation and split controls falsified privileged sparsity |
| C5. Thermodynamic perturbation state | one or two MMP response gains | centered support contrasts | relative free-energy responses can transfer | reject: no distant passive coverage and no held-component direction gain |
| C6. Bounded rank-edit policy | STOP or a few swaps | finite actions | source tasks teach limited interventions | defer: RIP observables failed and the action is candidate-list dependent |
| C7. Assay-state adapter | a small context state | assay categories are finite | coherent assays reduce noise | auxiliary only: noise fell, assignment-dependent adaptation did not appear |

**INFERENCE.** C1 is the only candidate that is both untested in the current
record and aimed directly at the failed representation assumption. Its value is
high only if segment-conditioned channels pass the protein destructions and
move the k=3/5 learning curve left on held source components.

## Mathematical Formulation

Let `mu(p,x)` be the frozen, support-free, out-of-fold DTA prediction and
`r_i=y_i-mu(p,x_i)` the support residual. A source-learned interaction encoder
emits two scalar channels:

\[
\phi_\theta(p,x)=(\phi_1(p,x),\phi_2(p,x))\in\mathbb R^2.
\]

`phi` combines ligand motif features with a ligand-conditioned pooling of 32
ESM-2 sequence segments. It is an interaction representation, not a claim of a
physical contact map.

Remove the rank-null level from support evidence:

\[
\tilde r_i=r_i-\bar r_S,\qquad
\tilde\phi_i=\phi(p,x_i)-\bar\phi_S.
\]

The learned adaptation operator emits a bounded state:

\[
z_t=A_\psi(\{(\tilde\phi_i,\tilde r_i)\}_{i=1}^k,p,k),
\qquad z_t\in[-1,1]^2.
\]

The budget mask is fixed before training:

\[
M_1=(0,0),\qquad M_3=(1,0),\qquad M_5=(1,1).
\]

The adapted prediction is

\[
\hat y_q=\mu(p,x_q)+\hat b_t+
c\tanh\left(\frac{(M_k\odot z_t)^\top
[\phi(p,x_q)-\bar\phi_S]}{c}\right).
\]

`b_hat` is a separately scored shrunk level channel and cannot support a ranking
claim. The rank correction is exactly zero at k=1, after support removal, and
when centered residual evidence is exactly zero.

**INFERENCE.** k=1 contains no within-target contrast after the unknown level is
removed, so any claimed one-shot ranking state would be confounded. At k=3,
there are two independent centered observations and only one adapted
coordinate. At k=5, there are four independent centered observations and at
most two coordinates. No high-dimensional target latent is hidden behind a
sparsity penalty.

## Trainable Path

**Stage R0 - representation admission.** Train `phi_theta` on fit components by
backpropagating held-query pairwise proper loss through a fixed empirical-Bayes
state solve. The solve is an instrument and baseline, not the final mechanism.
The frozen base is unchanged.

**Stage R1 - learned adaptation.** Only after R0 passes, freeze the admitted
representation and train `A_psi` on fit episodes. Parameterize the learned state
as a zero-initialized correction to the empirical-Bayes state so setting the
correction to zero recovers that baseline exactly. PIRS is admitted only if the
learned operator improves held-probe ranking beyond the fixed solve.

**Stage R2 - optional joint shaping.** Only after R1 passes, test joint base and
channel training under an explicit support-free non-degradation constraint.
This stage cannot rescue an R0 or R1 failure.

The source objective is

\[
L=L_{pairwise\ proper}+\lambda_{rank}L_{rank}
+\lambda_{state}L_{state\ stability}
+\lambda_{cf}L_{correct>deranged}
+\lambda_{null}L_{exact\ no-op}.
\]

All weights and thresholds are selected on fit components only. Probe is a
held-component development gate. `locked` and recipient labels remain sealed.

## Prior-Art Boundary

**FACT.** ActFound learns pairwise differences and meta-learns across assays,
but its released implementation skips inner-loop parameter updates when support
size is at most five. At this budget its prediction remains support-anchor and
similarity based.

**FACT.** FS-CAP averages support encodings into an unrestricted assay context;
MHNfs enriches predictions using context retrieval; MetaDTA/CNP/ANP and graph
neural processes decode broad task contexts. None gives a measured two-degree
recipient state with structural k=1 silence.

**FACT.** AdaMBind uses MAML-style adaptation and task scheduling. ADKF-IFT and
deep-kernel transfer meta-learn representations around analytic per-task GP
solves. They are mandatory conceptual baselines, but do not establish this
protein-conditioned, low-cardinality intervention on target-disjoint k<=5 DTA.

**FACT.** HyperPCM already uses a protein-conditioned hypernetwork for zero-shot
unseen-target DTI. DrugBAN and PSICHIC already learn joint or interpretable
protein-ligand interaction representations, including sequence-derived
interaction fingerprints. Protein conditioning, bilinear interaction encoding,
and interaction channels are therefore prior art rather than PIRS contributions.
HyperPCM does not infer a measured recipient state from k<=5 labels, while
DrugBAN and PSICHIC are support-free predictors; these distinctions are relevant
only if the PIRS support state and learned operator pass their causal gates.

**INFERENCE.** PIRS differs from MAML/ANIL/AdaMBind because no network weights
are adapted; from KRR/GP/ridge because the final support-to-state map is learned
and must beat the analytic solve; from retrieval/TRACE because predictions are
not weighted sums of support labels; from FS-CAP/MetaDTA because the state has a
declared dimension and acts only through two auditable channels; and from BIR
because the channels depend jointly and nonlinearly on protein segments and the
query ligand rather than on a target-independent linear ligand projection.
R0 alone cannot support an innovation claim: it can admit only a representation
for R1. The prospective contribution is the structurally budget-matched,
support-identified intervention state and an operator that improves over
analytic inference on that same representation.

Relevant literature:

- ActFound: https://doi.org/10.1038/s42256-024-00876-w
- FS-CAP: https://doi.org/10.1021/acs.jcim.4c00485
- AdaMBind: https://doi.org/10.1038/s41467-026-70554-5
- MHNfs: https://doi.org/10.1021/acs.jcim.4c02373
- Graph neural processes for molecules: https://doi.org/10.1186/s13321-024-00904-2
- ADKF-IFT: https://arxiv.org/abs/2205.02708
- BioBridge: https://doi.org/10.1002/advs.202506404
- Negative-transfer meta-learning: https://doi.org/10.1038/s41598-025-22058-3
- HyperPCM: https://doi.org/10.1021/acs.jcim.3c01417
- DrugBAN: https://doi.org/10.1038/s42256-022-00605-1
- PSICHIC: https://doi.org/10.1038/s42256-024-00847-1

## Preregistered Source-Only Gates

### R0A - harness positive control

Inject a two-channel target state into fit/probe features at measured noise.
Require monotone recovery from k=1 to k=3 to k=5, exact k=1 rank silence, and a
k=5 correct-vs-wrong CI lower bound above 0.005. A failed synthetic control
invalidates the harness; no real-data null is interpreted.

### R0B - representation generalization

Train representation statistics and weights on fit only. On scaffold-disjoint
probe targets, require all of:

1. full-data oracle use of the two learned channels has CI-gain lower 95% above
   0.005 in pooled and support-query Tanimoto-below-0.35 cells;
2. k=3 and k=5 empirical-Bayes state gains are monotone, with k=5 lower 95%
   above 0.005;
3. the learned channels beat matched target-independent, pooled-protein,
   random, protein-zero, protein-shuffle, and segment-transplant channels;
   an orthogonal rotation of both dense coordinates must be prediction-invariant;
4. correct support beats wrong-target support and centered-residual
   derangement; and
5. at least 47 independent components contribute to the pooled primary gate.

The oracle requirement tests whether the representation contains a
low-dimensional adaptation object. The k-shot requirement separately tests
whether that object is recoverable at the deployment budget.

### R1 - learned operator

If and only if R0 passes, require `A_psi` to beat the identical admitted
representation with empirical Bayes, ridge, KRR, and a frozen/random operator
at k=3 and k=5. Require positive component-bootstrap lower bounds for absolute
gain over the frozen base and for learned-minus-analytic gain. The correct arm
must beat wrong target, random support, residual derangement, coordinate
permutation, support removal, and label-noise dose controls. Gains confined to
Tanimoto at least 0.35 fail the mechanism.

### R2 - complete breakthrough

A major breakthrough requires R0 and R1 to pass on fit-to-probe transfer, no
support-free degradation, positive CI/NDCG@10/pairwise-proper results, reduced
negative-transfer rate, and exact no-op under support removal. Only then may
core code be copied once to `model/`, utilities once to `script/`, tested,
hashed, frozen, and left unchanged.

## Ablation Ladder

| Rung | Configuration | Exact fallback / claim |
|---|---|---|
| A0 | frozen base | support-free reference |
| A1 | + shrunk level only | rank-null calibration |
| A2 | + target-independent linear channels + EB | BIR/IDA-style baseline |
| A3 | + segment-conditioned channels + EB | representation-shaping claim |
| A4 | + learned `A_psi` correction | learned adaptation beyond analytic EB |
| A5 | + counterfactual objective | correct evidence specificity |
| A6 | + bounded hard no-op | intervention safety, not the core gain |

If A3 fails, A4-A6 are not run. If A4 is matched by EB, PIRS is a representation
baseline rather than a meta-adaptation mechanism and is rejected.

## Maximum Scientific Risk

**HYPOTHESIS.** Coarse ESM segments may not localize binding environments, and
scalar affinities may contain no transferable signal linking a random k<=5
support set to chemically distant query ranking. A fit-trained network can also
manufacture apparent correctability by degrading the support-free base.

**FACT.** Protein destructions, fixed-base R0, low-similarity reporting, and the
support-free non-degradation rule are designed to expose those failures. If R0
fails with a valid synthetic control, the current passive data do not justify
another compact-state architecture; the honest deliverable is the measured
upper bound and the need for diagnostic/prospective measurements or richer
binding-state observations.

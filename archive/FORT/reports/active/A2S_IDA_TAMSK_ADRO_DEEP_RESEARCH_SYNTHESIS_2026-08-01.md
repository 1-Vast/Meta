# A2S-DTA Deep Research Synthesis: IDA, TAMSK, and ADRO

Date: 2026-08-01

Scope: source-only research audit. No locked-source or recipient affinity label was opened.

Decision: `NO_FINAL_MECHANISM_ADMITTED`
Epistemic tags: **FACT**, **INFERENCE**, and **HYPOTHESIS** have their usual strict meanings below.

## Executive decision

**FACT:** The balanced ChEMBL v2 source diagnostic remains
`NO_GO_INFORMATION_NOT_ADMITTED`. Its real label and assignment intervals cross
zero, while its synthetic label channel is strongly detectable and a high-data
source oracle shows positive ranking headroom.

**FACT (reported, not independently verified):** A separate BindingDB branch
reports that fixed Tanimoto KRR produces a large k=5 ranking gain and that TAMSK
adds a much smaller gain over a static kernel mixture. The TAMSK source code,
data and split manifests, per-query predictions, component metrics, bootstrap
draws, gate definitions, tests, and frozen hashes are not present in the
accessible workspace.

**INFERENCE:** These results are not contradictory. They concern different
estimands. The ChEMBL v2 task is passive, chemically distant, frequently
cross-assay, and provenance-closed. The BindingDB task randomly samples support
within a target, lacks assay/document closure, and is likely much more favorable
to local ligand SAR smoothing. The most coherent current explanation is:

> Support-label assignment information is stratum-dependent. It can be strong
> in chemically local or provenance-linked episodes and absent or too weak to
> detect in the current global scaffold-cold passive construction.

**INFERENCE:** TAMSK is therefore valuable as an information probe and mandatory
local-SAR baseline, but it cannot be the requested final mechanism: its query
correction is a convex mixture of closed-form KRR experts, its learned object is
essentially a task-conditioned bandwidth/MKL selector, and protein conditioning
is not load-bearing in the described router.

**HYPOTHESIS:** IDA asks the most important question that remains inside the
fixed passive task: can source meta-training shape the base residual geometry so
that k=3/5 labels recover one or two stable rank-relevant directions? The current
IDA equations do not yet establish this. They contain scale, rotation, k=1,
base-degradation, OOF, and closed-form-posterior loopholes that must be closed
before implementation can count as a mechanism test.

**HYPOTHESIS:** ADRO remains the highest-potential long-term scientific idea,
because it learns what to measure as well as how to interpret the result. It is
not the immediate answer to the present prompt: the exploration contract fixes
passive support, whereas ADRO changes the observation process. It must remain a
separate future task and dataset track.

No model should be promoted to `model/`, no recipient labels should be opened,
and no positive-mechanism claim should be pushed from the present evidence.

## 1. Materials and evidence status

The four supplied documents were read in full or by complete section coverage:

| Document | SHA-256 | Evidence role |
|---|---|---|
| `A2S_DTA_INDEPENDENT_COMMITTEE_REVIEW_2026-08-01.md` | `85B4844F54E9C8CE1E71E682035F2D48D0299DFEEEE5451471EDB0D3C8650483` | Detailed CMAL code/run audit and pre-v2 recommendations |
| `A2S_EXPLORATION_PROMPT.md` | `C40E729B05795055CD3D1377DBD61048CB5AA0B45F0A1035F5FFCBE0BF749EFF` | Fixed passive-task and admissibility contract |
| `A2S_IDA_MECHANISM_DESIGN.md` | `81C52D2B997BF0B8F2584A21805A69746B5A9C9BB3A49FF1FF7F03BB3C064368` | New identifiability-shaping hypothesis |
| `TAMSK_SOURCE_ONLY_FEASIBILITY_HANDOFF_2026-08-01.md` | `415FCFB466012DD61E4FF1EE4DA182A20B8CDA9F164A037C256AB371AFB78361` | Document-only report of a separate BindingDB branch |

Repository evidence cross-checked:

- `reports/active/A2S_POST_REVIEW_V2_GATE_DECISION_2026-08-01.md`
- `reports/active/a2s_source_information_gate_v2_2026-08-01.json`
- `reports/active/a2s_source_information_gate_lock_v2_2026-08-01.json`
- `reports/active/A2S_FINAL_PI_META_MECHANISM_REDESIGN_2026-08-01.md`
- `research/a2s_cmal.py`
- `research/a2s_information_gate.py`

**FACT:** The independent committee packet did not include the later balanced v2
diagnostic. Its CMAL code defects remain useful facts, but its statement that
support information was not measurable is superseded by the balanced v2
instrument. V2 recovered synthetic signal under balanced OOF folds and still
returned a real-data no-go.

**FACT:** Searches of `D:\FORT`, all accessible Git refs, the supplied Desktop
directory, `.codex/attachments`, and the local Temp directory found only the
TAMSK handoff text. No executable TAMSK artifact was found. The TAMSK metrics
must therefore be labeled `REPORTED / UNVERIFIED`, not independently reproduced
results.

## 2. Authoritative ChEMBL v2 result

**FACT:** V2 retained 68,782 of 185,591 source pKi rows. Fit/probe/locked roles
contain 222/110/107 homology components and 36,609/16,068/16,105 rows. Target,
homology, document, and assay overlap across roles is zero. Five OOF folds each
hold approximately 20% of fit rows.

| k | Probe components | Delta-label rank-loss 95% CI | Delta-assignment rank-loss 95% CI | Synthetic Delta-label | High-data oracle LCB |
|---:|---:|---:|---:|---:|---:|
| 1 | 69 | `[-0.00426, +0.01271]` | undefined | `+0.31635` | `+0.02647` |
| 3 | 45 | `[-0.00865, +0.01185]` | `[-0.00087, +0.01366]` | `+0.22976` | `+0.06979` |
| 5 | 36 | `[-0.00867, +0.01665]` | `[-0.01862, +0.00359]` | `+0.38702` | `+0.07097` |

**INFERENCE:** Query ranking is improvable, and the probe can recover a strong
label channel. What is not admitted is stable incremental information in the
present passive support assignment under this model class and available power.
This is a bounded empirical null, not an information-theoretic impossibility.

**FACT:** Locked-source and recipient labels remain sealed. The v2 report is a
source development diagnostic, not confirmatory recipient evidence.

## 3. Why the reported TAMSK positive does not reverse v2

### 3.1 The tasks differ

| Property | ChEMBL balanced v2 | BindingDB TAMSK branch |
|---|---|---|
| Records | 68,782 retained source pKi rows | 24,815 deduplicated pairs reported |
| Targets | Role-closed source registry | 264 reported |
| Support policy | Passive ordered construction | Random k within target |
| Support-query geometry | Nearest Tanimoto about 0.223; about 81% scaffold-cold queries | Not quantified in the handoff; target-random sampling favors local analogs |
| Provenance control | Target/homology/document/assay role overlap zero | Assay, construct, qualifier, replicate, and document/source closure unavailable |
| Information result | No positive real lower bound | Fixed KRR reported strongly positive at k=5 |
| Confirmation status | Locally executable development diagnostic | Document-only; underlying artifacts absent |

**INFERENCE:** A support compound may be globally absent from source training and
still be a close analog of a query compound on the same recipient target.
"Globally new scaffold query" does not establish support-query scaffold coldness.
The BindingDB handoff does not report recipient support-query Murcko overlap,
nearest Tanimoto, MMP connectivity, or shared assay/document provenance.

### 3.2 Reported TAMSK effects

**FACT (reported, unverified):** The handoff gives the following k=5 locked-test
macro values over 34 evaluable protein components:

| Method | RMSE | CI | Spearman | NDCG@5 |
|---|---:|---:|---:|---:|
| Frozen base | 1.7363 | 0.5141 | 0.0372 | 0.3777 |
| Fixed KRR | 1.4344 | 0.5579 | 0.1589 | 0.4661 |
| Static mixture | 1.4357 | 0.5599 | 0.1643 | 0.4705 |
| TAMSK | 1.4243 | 0.5612 | 0.1678 | 0.4742 |

**FACT (reported, unverified):** TAMSK minus static mixture paired component
effects are CI `+0.00147 [0.00034, 0.00268]`, Spearman
`+0.00360 [0.00057, 0.00688]`, NDCG@5
`+0.00408 [0.00096, 0.00765]`; the RMSE interval crosses zero.

**INFERENCE:** Fixed KRR supplies almost all of the useful adaptation. The learned
router contributes about 3% of the fixed-KRR CI gain relative to the frozen base.
That may be statistically nonzero in the reported bootstrap, but no practical
effect floor was preregistered and the raw evidence is unavailable.

### 3.3 TAMSK mechanism classification

The reported TAMSK prediction is

\[
\hat y_q=f_0(p,x_q)+\sum_{m=1}^{4}w_m(S)\,
k^{(m)}_{qS}(K^{(m)}_{SS}+\lambda_m I)^{-1}r_S.
\]

**FACT:** The query corrections are closed-form residual KRR experts. The router
reads centered residual-kernel alignment, leave-one-support-out consistency, and
k, and outputs one task-level mixture vector. It has no explicit protein input
in the handoff definition.

**INFERENCE:** TAMSK qualifies operationally as amortized task adaptation, but
its scientific object is learned kernel-scale/model selection. It is closest to
centered-alignment MKL, MetaVRF/MetaKernel, and ADKF-IFT, not a new adaptation
family. Query dependence is supplied by KRR similarity; the learned router is
not itself query-specific.

**INFERENCE:** Matching four experts to `k-1=4` does not prove identifiability.
The four convex weights have three degrees of freedom, their expert outputs are
correlated, and at k=3 only two centered residual contrasts exist. Expert count
and observation rank are not interchangeable quantities.

**INFERENCE:** A derangement no-op loss is useful, but it is not a structural
no-op guarantee. With positive convex mixture weights, the model may have no
exact zero action unless a separate hard gate or null expert exists. Source code
is required to determine the actual behavior.

**Decision:** `MANDATORY_BASELINE_AND_INFORMATION_PROBE`, not final mechanism.

## 4. IDA: the right question, an incomplete mechanism

### 4.1 What IDA genuinely changes

For centered rank features, the relevant observation is

\[
H_k r_S=H_k\Phi_Sz_t+H_k\epsilon,
\qquad H_k=I-k^{-1}\mathbf 1\mathbf 1^\top.
\]

A scale- and coordinate-aware information matrix is

\[
J_S=\sigma^{-2}\Sigma_z^{1/2}\Phi_S^\top H_k\Phi_S
\Sigma_z^{1/2},
\]

with effective information dimension

\[
d_{\mathrm{info}}(S)=\operatorname{tr}[J_S(I+J_S)^{-1}]\le k-1.
\]

**INFERENCE:** Learning `Phi` can change the spectrum of `J_S`; IDA therefore can
change the geometry of a restricted latent model. It does more than fit a new
estimator on frozen features. It cannot increase the `k-1` contrast ceiling or
create missing assay, chemistry, or provenance relations.

The honest core claim is therefore:

> **HYPOTHESIS:** Source meta-training can make one or two rank-relevant residual
> directions empirically recoverable on unseen source components without
> degrading the strongest support-free predictor.

It is not "identifiable by construction."

### 4.2 Mathematical defects that must be repaired

1. **Scale loophole - INFERENCE.** The proposed inverse-Gram penalty and code L2
   penalty can be reduced by scaling `Phi` up and `z` down. Hard basis
   normalization, rank-space centering, code whitening, and orthogonality or an
   equivalent invariant parameterization are required.
2. **Rotation dependence - INFERENCE.** Per-coordinate `rho_j` changes under a
   basis rotation even when predictions do not. Admission must use the invariant
   information spectrum and held-out code recovery, not training coordinates.
3. **Intercept contamination - INFERENCE.** Conditioning must be computed on
   `H_k Phi_S`, not on an uncentered matrix containing the constant column.
4. **k=1 is not silent - FACT/INFERENCE.** With one observation, rank coefficients
   and intercept are confounded. Ridge shrinkage does not force rank coefficients
   to zero, and `Delta_psi` can emit a rank code. A hard graph-level mask
   `z_rank=0` and `g_rank=0` is mandatory at k=1.
5. **k=3 is saturated - FACT/INFERENCE.** Intercept plus two rank directions use
   all three observations and leave no residual degree of freedom to test
   evidence consistency. A safer fixed budget is
   `m_active(1)=0`, `m_active(3)=1`, `m_active(5)<=2`.
6. **Training-code optimism - INFERENCE.** Free abundant-source codes are learned
   objects, not ground truth. Their dispersion can be manufactured. Recovery
   must be measured from k-shot support against an independent full-data oracle
   code on unseen source components.
7. **Support conditioning is not query utility - INFERENCE.** A well-conditioned
   support basis may have no useful action on query compounds. The diagnostic
   must include the query predictive covariance/action, not only
   `(Phi_S^T Phi_S)^{-1}`.
8. **Collapse loophole - INFERENCE.** Penalizing low-reliability code directions
   can collapse all codes into the base rather than concentrate useful signal.
   A positive held-out rank result is the only evidence that this did not occur.

### 4.3 Base degradation can manufacture adaptation

**INFERENCE:** Joint training can deliberately weaken the support-free base and
place otherwise predictable signal in `Phi z`. Adaptation then looks strong only
relative to its own damaged base. L2 on source codes does not prevent this.

IDA must satisfy all of the following:

- compare with an independently optimized, equal-budget support-free reference;
- preregister a one-sided no-support non-inferiority margin;
- beat the strongest reference, not only IDA's internal base;
- reject the run if base degradation is comparable to the adaptation gain;
- use component-held-out support-free risk, not a universal training-RMSE floor,
  for early stopping.

### 4.4 The OOF proposal is not strict OOF

**INFERENCE:** A shared affinity encoder trained on all component labels has read
the held-out fold, even if its readout head is excluded. Head-only cross-fitting
does not make representation-derived residuals strictly OOF.

Required correction:

- every label-trained layer of a fold model excludes that fold;
- the common basis uses fixed label-free inputs or a representation whose
  provenance is explicitly separated from OOF scalar residual evidence;
- probe targets use only fit-role-trained frozen models;
- OOF/probe residual distribution parity is reported.

### 4.5 Closed-form core and novelty boundary

IDA R3 computes

\[
\hat z_{EB}=(\Phi_S^\top\Phi_S+\Lambda)^{-1}\Phi_S^\top r_S.
\]

**FACT:** This is low-rank ridge / a linear-Gaussian empirical-Bayes posterior
mean. `L_ident` changes the learned representation but not the closed-form inner
adaptation class.

**INFERENCE:** If the later learned correction `Delta_psi` has no independent
held-out gain, IDA is an identifiability-regularized deep-feature EB baseline.
It does not meet the project's final non-closed-form mechanism requirement.

If `Delta_psi` does add value, the operator is a finite set-to-code amortized
inference network, adjacent to CNP/ANP and MetaFun. Its defensible increment is
the load-bearing budget-shaped coordinate system, not the existence of a set
encoder.

**Decision:** `CONDITIONAL_DIAGNOSTIC_BRANCH`. Do not call it the final mechanism
unless both the shaping increment and the learned-over-EB increment pass.

## 5. ADRO under the corrected fixed-task contract

**HYPOTHESIS:** ADRO learns a measurement policy and a finite response-mode
executor. It attacks the information-acquisition bottleneck rather than only
reparameterizing the residual field.

**INFERENCE:** This remains the highest-potential long-term biological direction,
especially if its executor uses a small, protein-conditioned, assay-anchored
response basis and a structural null mode.

**FACT:** The current exploration contract fixes passive support episodes. ADRO
requires a legal candidate pool and sequentially selected measurements. It is a
different task contract and cannot be used to rescue the current passive result.

ADRO may proceed only on a separate track with same-assay candidate pools,
strict label accounting, random/diversity/D-optimal/active-KRR baselines, and
preferably prospective recipient measurements. Retrospective ChEMBL simulation
would support only a finite-pool acquisition claim, not real campaign causality.

## 6. Revised research question

The previous binary question, "does support information exist?", is now too
coarse. The most valuable unknown is:

> Under which predeclared chemical, assay, and provenance relations does correct
> support assignment contain transferable ranking information, and can source
> meta-training either shape a passive residual coordinate system or design an
> active diagnostic measurement so that this information exceeds strong local
> analytic baselines?

This separates three distinct objects:

1. **Local information admission:** fixed KRR and derangement determine whether
   target-local SAR information exists in a stratum.
2. **Passive geometry shaping:** corrected IDA tests whether source learning can
   make a globally useful low-dimensional correction recoverable without
   weakening the base.
3. **Active information creation:** ADRO tests whether measurement selection can
   make a finite response state distinguishable under the same label budget.

TAMSK addresses (1) and a narrow form of task-adaptive analytic routing. IDA
addresses (2). ADRO addresses (3). Combining all three now would erase the
load-bearing claim and is not justified.

## 7. Minimum next experiment: two gates, no model sprawl

### Gate 0 - recover or retire the TAMSK evidence

Before citing the BindingDB result as experimental evidence, require:

- exact source and `kappa_gamma` implementation;
- data, preprocessing, split, and episode manifests;
- code/data/config/environment hashes;
- per-query predictions and support IDs;
- component metrics and bootstrap indices;
- complete seven-gate JSON and all nineteen tests;
- a clean-room rerun reproducing the summaries within frozen tolerances.

Failure means the TAMSK numbers remain a hypothesis-generating handoff only.
Reimplementing a plausible TAMSK from prose does not reproduce the old result.

### Gate 1 - information-locality crossover

Use one frozen base, one fixed KRR implementation, one true residual derangement,
and one component-bootstrap pipeline in a 2x2 comparison:

| Corpus | Episode A | Episode B |
|---|---|---|
| ChEMBL | Original v2 passive/ordered | BindingDB-style random within target |
| BindingDB with restored metadata | Current random | Provenance-closed and assay/locality matched |

Predeclare and report exact-assay, same-document, same-scaffold, MMP-connected,
scaffold-cold, and support-query Tanimoto strata.

**HYPOTHESIS:** If the effect follows locality/shared provenance rather than the
dataset name, TAMSK is a local SAR/provenance smoother. If it survives closure
and scaffold-cold or exact-transformation strata, a narrower transferable claim
becomes credible.

This gate is cheaper and more informative than training either IDA or ADRO. It
determines which task the next model would actually solve.

### Gate 2 - corrected IDA R2/R3 falsification

Run only after Gate 1 defines the admissible passive stratum. R2 and R3 must use
the same base, architecture, parameter count, optimizer, episodes, EB solve, and
active dimensions. The only change is an invariant identifiability-shaping loss.

Required fixed design:

- balanced v2 fit/probe roles; locked source remains sealed;
- strict component-cross-fitted residuals;
- centered, normalized, whitened/orthogonal basis;
- `m_active={0,1,2}` for k={1,3,5};
- strongest independent support-free reference;
- synthetic positive controls and preregistered MDE;
- at least 2,000 paired component bootstrap draws.

Measure on unseen probe components:

\[
\Delta_{shape}=L_{rank}(R2)-L_{rank}(R3),
\]

\[
\Delta_{abs}=L_{rank}(f_{ref})-L_{rank}(R3),
\]

\[
\Delta_{assign}=L_{rank}(R3,\widetilde S)-L_{rank}(R3,S).
\]

Also measure invariant `d_info`, k-shot versus full-data-oracle code recovery,
support-half reproducibility, query correction variance, changed-pair fraction,
and no-support base non-inferiority.

All three deltas must have paired component 95% lower bounds above a
preregistered practical/MDE threshold at k=3 and k=5. k=1 rank action must be
exactly zero. Residual-null must be bit-level no-op. Derangement, protein
shuffle/zero, target shuffle, and label noise must destroy the claimed gain;
query permutation, distractor, subset, and library-size controls must pass.

If R3 passes, only then test

\[
\Delta_{learned}=L_{rank}(EB\text{-only})-
L_{rank}(EB+\Delta_\psi).
\]

If `Delta_learned` does not clear the same material threshold, IDA remains a
closed-form baseline and cannot be the final paper mechanism.

### Stop rules

Stop IDA before locked-source evaluation if any is true:

1. training `tau/sigma` rises but unseen-component code recovery does not;
2. `Delta_shape`, `Delta_abs`, or `Delta_assign` fails;
3. support-free base degradation explains the gain;
4. k=1 produces a rank action;
5. G0/derangement/protein-shuffle preserves the gain;
6. the correction is an episode constant, base-score calibration, or ligand-only
   similarity action;
7. available component count cannot detect the preregistered minimum effect.

Locked source is opened once only after the protocol, MDE, seeds, baselines,
controls, and stop rules are frozen. Recipient labels remain sealed until that
locked-source evaluation passes.

## 8. TAMSK external confirmation kill criteria

If TAMSK is reconstructed as a baseline or narrow application candidate, it must
pass all of the following on a genuinely untouched corpus:

1. **Artifact gate:** complete immutable artifacts and clean-room reproduction.
2. **Information gate:** correct fixed KRR beats frozen base and derangement by
   a preregistered material lower bound at k=3/5; noise destroys the effect
   monotonically.
3. **Nearest-method gate:** TAMSK beats static MKL, CKA-linear/NNLS, nested
   support-only gamma/lambda selection, LOO expert selection, hierarchical
   empirical-Bayes MKL, and feasible ADKF-IFT.
4. **Router necessity:** weights vary reproducibly across components and predict
   held-out oracle expert utility; router permutation or a deterministic selector
   must not tie it.
5. **Primary outcome:** co-primary CI and NDCG@10 lower bounds exceed both a
   practical threshold and MDE80, with RMSE non-inferiority. Merely exceeding
   zero is insufficient.
6. **Shortcut/invariance:** residual-null exact no-op; derangement no benefit;
   support/query permutation, distractor, subset, and library-size stable.
7. **Protein claim:** matched protein-shuffle/zero must cause a material routing
   decrement. Otherwise the result is ligand/activity adaptation, not
   protein-conditioned DTA.
8. **External replication:** pass on a second independent source or endpoint
   before any sealed A2S recipient is opened.

Failure of the nearest-method or router-necessity gate classifies TAMSK as a
useful meta-tuned KRR baseline, not a failed experiment.

## 9. Novelty boundary and relevant prior art

The following are mandatory comparison neighborhoods, not optional citations:

- Cortes, Mohri, and Rostamizadeh, centered-alignment MKL, JMLR 2012:
  <https://jmlr.org/papers/v13/cortes12a.html>
- MetaVRF, ICML 2020:
  <https://proceedings.mlr.press/v119/zhen20a.html>
- MetaKernel, task-adaptive meta-kernel learning, DOI
  `10.1109/TPAMI.2022.3154930`.
- ADKF-IFT, meta-learned molecular deep kernels:
  <https://arxiv.org/abs/2205.02708>
- FS-CAP, few-shot continuous compound activity prediction, DOI
  `10.1021/acs.jcim.4c00485`.
- MetaDTA:
  <https://openreview.net/forum?id=yzlif16IASM>
- AdaMBind, task-adaptive DTA:
  <https://www.nature.com/articles/s41467-026-70554-5>
- CNP and MetaFun for finite support-to-code/function inference:
  <https://proceedings.mlr.press/v80/garnelo18a.html> and
  <https://proceedings.mlr.press/v119/xu20i.html>.

**INFERENCE:** IDA's potential novelty is not a low-rank code, EB solve, or set
encoder. It is a load-bearing, invariant objective that demonstrably reshapes
unseen-target residual geometry at a declared label budget while preserving the
support-free predictor.

**INFERENCE:** TAMSK's potential narrow contribution is not a Tanimoto kernel,
MKL, CKA, or meta-kernel. At most it is the jointly validated protocol of
residual-assignment evidence, a finite action family, and derangement-aware
routing. That contribution disappears if CKA/LOO/nested-CV selection ties it.

## 10. Final committee ranking

| Direction | Immediate role | Scientific potential | Current admission |
|---|---|---:|---|
| CMAL | Failed mechanically active baseline | Low | Frozen |
| TAMSK | Local-SAR information probe and mandatory KRR/MKL baseline | Narrow/moderate | Reported, not reproducible locally |
| Corrected IDA | Passive geometry-shaping falsification branch | Moderate/high for ML methodology | Not admitted |
| ADRO | Separate active-measurement research track | Highest long-term | Not testable under current passive contract |
| TAMSK+IDA hybrid | None | Low incremental clarity | Rejected for now |

The key positive development is conceptual, not yet algorithmic: the evidence
now points to an **information-admission boundary** controlled by chemical
locality and provenance, while IDA identifies a legitimate way to test whether
that boundary can be moved by meta-training. Neither constitutes the required
learned transferable adaptation mechanism until the registered gates pass.

Final status: `RESEARCH_DIRECTION_CLARIFIED; KEY_POSITIVE_MECHANISM_NOT_YET_ESTABLISHED`.

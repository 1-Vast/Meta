# Post-R3 multi-agent deep-model route selection

**Date:** 2026-07-28  
**Status:** design selection only; no new model run, affinity-label read, provider call, acquisition, or
assay.  
**Program state preserved:** `R3_ESTIMATOR_INSENSITIVE_NO_DECISION`;
`SIGNAL_PRESENT_EVIDENCE_INSUFFICIENT`; `sealed_test_consumed=false`.

## 1. Question and binding assumptions

This review asks for a learnable deep model that could address strict simultaneous target-cold and
ligand-cold affinity prediction after the K-LBP R3 no-decision. It does not treat a data design,
closed-form estimator, LLM compiler, or larger encoder as a model candidate by itself.

The following facts remain binding:

1. R3 failed before an effect statistic was interpretable. It did not show that the R1 coordinate has
   zero affinity signal.
2. The current identified label support is kinase-only. Prediction still requires approximately 423
   independent multi-family components, `PA2 >= 0.5 pK`, and at least 40 scaffold-diverse query ligands
   per target after every firewall.
3. Pooled ESM-2 reparameterization, pose/pocket rescue, support-posterior/meta-learning, unrestricted
   RAG, direct LLM prediction, and Transformer/Mamba growth remain closed.
4. Task A and Task B remain separate. The selected primary route is Task A: sequence plus ligand
   structure at inference. LEXOR-MC can later define a separate Task B arm only after MC0--MC2.

There is a formal support limitation behind item 2. Two data-generating processes can agree on every
kinase training observation and assign opposite ligand orderings to an unseen protein family. Their
training likelihoods are identical. No deep architecture can identify which extrapolation is correct
without multi-family supervision or an independently validated target-side variable.

## 2. R3 failure diagnosis

K-LBP parameterizes the correction as `Theta = gamma * a * c^T`. At `gamma = 0`,

```text
d mu / d a = 0
d mu / d c = 0
d mu / d gamma depends on arbitrary, unidentified a and c
```

Thus `a` and `c` exist only under the alternative and the Fisher information is singular at the null.
Requiring every alternating-GLS direction fit to converge in the S1-null regime is a non-regular
execution condition. The observed 3/5 converged held folds are consistent with this structural problem;
more iterations or another seed would not repair it. The frozen R3 result stands and must not be amended.
A successor must use a parameterization with a unique null and an inference method valid at that null.

## 3. Candidate screening

| candidate | new identifying information | learnable deep core | strict dual-cold compatible | cheapest falsification | decision |
| --- | --- | --- | --- | --- | --- |
| Randomized Multi-fidelity Assay-Monotone Operator Network (`R-MAON`) with a regular-null operator | randomized, inactive-retaining single-dose/censored observations plus independent dose-response calibration | sequence encoder, ligand encoder, direct interaction operator | yes, conditional on a compliant prospective panel | topology/power audit, then null-score calibration | **selected** |
| Cross-Provenance Replicate-View DTA (`CPRV-DTA`) | same cell measured by genuinely independent laboratories | shared interaction bottleneck plus monotone site heads | yes, but approximately doubles measurement cost | blinded cross-site residual-rank pilot | conditional runner-up |
| LEXOR-MC fixed-channel scalar network | source-bound qualitative mechanism evidence | fixed mechanism channels plus learned scalar gates | Task B only | MC0--MC2 and slot-coverage audit | retain as separate Task B candidate, not stacked initially |
| state-marginal energy/MoE model | no validated new labels | learnable state router and energy heads | nominally | state destruction | reject: free-routing equivalence and weak state identifiability |
| causal/proximal matrix completion on the current graph | no valid instrument/proxy or randomized anchor | neural bridge | no | proxy-validity audit | reject: assumptions unavailable; cannot create components |
| active acquisition policy alone | chooses measurements but is not the predictor | trainable policy | only after assays | random-policy comparison | retain only as a cost tool; not an innovation |
| larger attention/diffusion/Mamba model | none | yes | unproven | ligand-only and target destruction | reject: reopens a closed architecture-only route |

The candidates are not stacked. `CPRV-DTA` becomes reasonable only if two independent laboratories are
already available. LEXOR-MC remains an optional Task B information source and does not supply labels.

## 4. Selected model: R-MAON

### 4.1 Biological and statistical hypothesis

A systematic low-fidelity panel contains target-specific ligand reordering even when its absolute
single-dose values are assay-dependent. Within a fixed target and assay, an unknown monotone measurement
map preserves order outside frozen saturation/tie bands. A randomized, inactive-retaining panel makes
that order estimable without the historical hit-selection mechanism. Independent dose-response
measurements calibrate the latent ordering to one endpoint-specific affinity scale.

This continues the internally admissible AMOB mechanism on a new substrate; it does not rename the
`OPEN_DATA_INSUFFICIENT_FOR_AMOB` result. AMOB's external-data signal was real but could not pass the
assay/document firewall. The prospective design changes that missing information.

### 4.2 Standard encoders and exact interaction decomposition

Use compact, standard encoders rather than another architecture campaign:

```text
p_t = P_alpha(sequence_t) - E_component-balanced_train P_alpha(sequence)
l_d = L_beta(ligand_d)   - E_scaffold-balanced_train L_beta(ligand)
z_td = p_t^T Theta l_d

y_hat_td = b_eta(d) + w_bar^T l_d + z_td
```

`P_alpha` is a small sequence CNN and `L_beta` is either the existing descriptor MLP or a
parameter-matched standard ligand encoder. They are not claimed as innovations. Double centering gives
an exact functional-ANOVA interaction: `z` cannot represent a target intercept or a target-common
ligand-potency term. `b_eta` and the shared-global term remain explicit, strong nulls.

The low-fidelity loss uses only reliable within-target, within-assay order. A cross-fitted
ligand-only low-fidelity nuisance is a frozen offset in the pair score; it is not subtracted from the
observed percent-inhibition scale:

```text
L_ord = sum_(t,a,i,j) softplus(
    -r_taij * (
        b_lo(d_i) - b_lo(d_j) + z_ti - z_tj
    ) / tau_a
)
```

Here `r_taij` is the raw observed order outside the frozen tie/saturation band, `b_lo` is fitted without
the held component, and `tau_a` is assay-owned rather than an entity feature. This preserves the
unknown-monotone order contract while blocking `z` from relearning generic potency. pKi and pKd are
never pooled. High-fidelity dose-response records use an endpoint-specific continuous likelihood.
Selection probabilities are recorded; the randomized tranche defines the estimand, while any adaptive
tranche is training-only and propensity-adjusted.

### 4.3 Innovation accounting

**I1: assay-monotone multi-fidelity supervision.** The model learns a shared latent interaction from
randomized low-fidelity order and independent high-fidelity calibration without treating percent
inhibition as affinity or generating pseudo-affinity labels. Removing I1 means training on high-fidelity
records alone; replacing it with raw low-fidelity regression is the mandatory standard control.

**I2: regular-at-null direct target-chemical operator.** `Theta` is parameterized directly, with
`lambda * ||Theta||_* + epsilon * ||Theta||_F^2`, rather than as `gamma*a*c^T`. `Theta=0` is unique and
bit-identically recovers the shared-global null. Before interaction training, the maximum rank-1 score
at the null is tested by the operator norm of component-aggregated score matrices. Homology-component
and ligand-scaffold multiplier weights handle dependence; empirical `V_t` is not inverted.

No third model innovation is introduced. Binding-profile sentinels, two-site replication, randomized
sampling, EVMG, and M2LC are data or control infrastructure, not extra model modules.

### 4.4 Why this is a deep model

`P_alpha`, `L_beta`, and `Theta` are learned. The null-score gate freezes cross-fitted encoder features
only long enough to make the pre-training test regular. If that gate passes, the one-seed pilot trains
the compact encoders and direct operator under `L_ord + L_high_fidelity`; encoder-gradient detachment is
a required ablation. A convex direct-operator layer does not make the full predictor a closed-form
method.

## 5. Required prospective substrate

The model is not authorized on the current observational graph. A compliant substrate must combine:

1. **Mechanism-orthogonal multi-family roster:** mechanism or sequence variation must recur both across
   families and within families; taxonomy, study depth, and target popularity cannot determine the
   roster.
2. **Randomized multi-fidelity measurement:** retain inactive, censored, and failed curves; record the
   inclusion probability for every Stage-A and Stage-B cell. At least 40 high-fidelity query ligands per
   target remain in a preassigned random tranche.
3. **Binding-profile sentinel firewall:** a chemically disjoint sentinel library is held by a
   custodian and used only to create immutable target blocks. Profiles and correlations never enter the
   model.
4. **Independent provenance:** a blinded bridge pilot and a genuinely independent confirmation
   lineage. Plates, batches, and nominal documents are not provenance families.
5. **Final size:** approximately 423 or more independent multi-family components after all eight
   firewalls, with `PA2 >= 0.5 pK` and sufficient query depth. The 70--155 component range is only a
   train-only mechanism pilot.

## 6. Minimum kill gate

The first gate is deliberately cheaper than building the model or commissioning the full panel.

### G0-A: label-blind topology and information audit

Freeze the target roster, ligand libraries, randomization probabilities, site allocation, and
binding-profile blocks. Using no affinity outcome:

* MC/profile/homology blocking must leave a planned 70--155-component mechanism pilot and a credible
  expansion to approximately 423 predictive components;
* every target must retain at least 40 preassigned scaffold-diverse high-fidelity query ligands;
* the planned crossed design must have 80% power at the frozen 0.03 paired ranking floor under the
  empirical noise sensitivity analysis.

Failure returns:

```text
RMAON_G0_TOPOLOGY_OR_POWER_STOP
```

No assay or model implementation follows a failure.

### G0-B: regular-null synthetic score calibration

Before reading any new real label, re-use the frozen R3 S1/S2/S3/S5 constructions, 200 replicates, seed
1729, and empirical covariance-scale distribution. Cross-fitted deep features are frozen. Let `S_c` be
the centered score matrix for homology component `c` and
`T = ||sum_c S_c||_op`. Calibrate `T` by component/scaffold Rademacher multipliers.

Pass requires:

* null and S5 false rejection rate `<= 0.10`;
* power `>= 0.80` for the frozen effect magnitude in S1 and S2;
* S3 power no more than 0.15 below S1;
* median recovered operator magnitude divided by truth in `[0.8, 1.25]`;
* no optimizer-convergence criterion at the null.

Failure returns:

```text
RMAON_G0_NULL_SCORE_OR_RECOVERY_STOP
```

G0-A and G0-B are the complete minimum kill gate. A pass authorizes only a small blinded
cross-provenance reliability pilot, not affinity-model training or a prediction claim.

## 7. First paid and first model gates

**A0 reliability pilot:** at least 12 targets spanning at least six families and 16 ligands, randomized
across two independent provenance lineages. Primary contrast:

```text
rho(correct-target, target-centered residual order)
  - rho(wrong-target, target-centered residual order)
```

Its component/site bootstrap LCB95 must exceed zero, cross-site rank reliability must meet the frozen
reliability threshold, and measured dispersion must support the later 0.03 contrast. Otherwise:
`RMAON_A0_ORDINAL_SIGNAL_OR_POWER_STOP`.

**M1 train-only mechanism pilot:** only after A0, on 70--155 independent components. The direct-operator
score for the true target representation must beat wrong-target, matched-random, taxonomy-centroid,
ligand-only ordinal, raw-low-fidelity-regression, and high-fidelity-only controls. Failure stops before a
full deep training run.

**Predictive gate:** only on the approximately 423-component prospective substrate. Primary endpoint is
paired component-macro Spearman versus both B0 and shared-global, with LCB95 above zero and point gain at
least 0.03. It must collapse under target, ligand-order, assay, and source destruction and reproduce in
the independent provenance lineage. Sealed evaluation remains last.

## 8. Novelty boundary and recent literature

The claim is not "first multi-fidelity DTA", "first bilinear DTA", or "first cold-start neural model".
Theisen et al. already combined single-dose and dose-response data and prospectively measured new pairs,
but evaluated compound/compound-cluster cold splits rather than simultaneous target/ligand cold
isolation (Nature Communications 2024, DOI `10.1038/s41467-024-52055-5`). CS-DTA's reported strict
settings are separate warm, protein-cold, and drug-cold settings, not this program's simultaneous
eight-axis dual cold (DOI `10.3389/fchem.2026.1834317`). Recent noisy/inexact matrix-completion theory
improves estimation given side information; it cannot create missing target support (arXiv
`2605.17189`; DOI `10.1080/01621459.2024.2335591`).

The defensible prospective novelty, if every gate passes, is the exact combination:

> a randomized inactive-retaining unknown-monotone assay bridge, an exactly main-effect-orthogonal deep
> interaction, a regular-at-null direct operator test, and simultaneous target/ligand/provenance/profile
> cold evaluation.

That statement still requires a formal systematic novelty review before publication.

## 9. Decision

```text
SELECT_RMAON_FOR_G0_PREREGISTRATION__NO_MODEL_OR_PREDICTIVE_STAGE_AUTHORIZED
```

Immediate next work is a separate G0 preregistration and score-calibration runner. It must not alter the
frozen K-LBP R3 artifact or run R4. Prospective acquisition remains separately authorized and cannot be
inferred from this design selection.

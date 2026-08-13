# QPSMP Theory for Cold-Target Few-Shot Drug-Target Affinity Prediction

## 1. Objective

The objective is to learn one shared deep model from historical protein tasks and adapt it to an
unseen protein using only a small support set. The active core module is
**QPSMP-HyperSAR**, an amortized conditional meta-potential.

QPSMP is not an auxiliary diagnostic head. It is the retained scalar affinity path used for
zero-shot prediction, interaction quotients, and few-shot adaptation.

### 1.1 Active HyperSAR amendment

This amendment supersedes the affine-section/ridge adaptation claims in Sections 4.2, 6, 7, 8,
G3, and the old minimal implementation sketch wherever they conflict with the active model. Those
parts remain as historical comparator theory only.

For support residuals $r_i=Y_i-f_{\theta,0}(p,l_i,a)$, define a support-anchored level and
relative-shape channel

$$
f_{\mathrm{level}}(l\mid S)=\bar Y_S+g_\eta(\bar r,\operatorname{MAD}(r),\log(1+k))
\left[f_{\theta,0}(p,l)-\overline{f_{\theta,0}}{}_S\right],
\qquad 0\le g_\eta\le1.
$$

Thus $g_\eta=0$ is exactly the support-mean baseline, while $g_\eta=1$ transfers the zero-shot
relative ligand shape after cancelling its target-level bias.

The centered support evidence is encoded by a permutation-invariant neural map

$$
t_S=\mathcal H_\psi\left(
\{h_\theta(p,l_i,a),e_\xi(l_i,a),r_i-\bar r,|r_i-\bar r|\}_{i=1}^k
\right).
$$

The implementation includes an explicit residual-bound support moment

$$
m_S=\frac{\sum_i(r_i-\bar r)\,B_\psi[h_\theta(p,l_i),e_\xi(l_i)]}
{\sum_i|r_i-\bar r|+\varepsilon},
\qquad t_S=\tanh(\mathcal H_\psi(S)+C_\psi m_S).
$$

Consequently, a joint permutation of support triples preserves the code while a label-only
permutation generally changes it. This makes label--ligand binding part of the architecture rather
than an optional pattern that an unconstrained set encoder may ignore.

$t_S$ is structurally multiplied by the centered-evidence norm, so constant residual evidence and
$k=1$ give $t_S=0$ exactly. Unlike the retired affine section, $t_S$ is not constrained to the row
span of current support embeddings; it is amortized prior-dependent inference learned across source
targets.

Following the relative-recognition principle of PBCNet2.0, the support code is a reference rather
than the final adapter input. Each query obtains

$$
t_{S,l}=t_S\odot\left[1+\tanh R_\psi\left(
q_\psi(h_\theta(p,l),e_\xi(l))-a_\psi(t_S),
q_\psi(h_\theta(p,l),e_\xi(l))\odot a_\psi(t_S)
\right)\right].
$$

This lightweight Siamese comparison connects the support reference to both the drug and
protein--ligand interaction branches without assuming that every deployment pair has a resolved
complex structure. Because $t_{S,l}$ is multiplicatively anchored by $t_S$, $k=1$ and constant
centered evidence still give exact zero structural modulation.

To avoid forcing every ligand-specific effect through one target bottleneck, the same module also
uses a learned reference--query matching correction

$$
c_S(l)=\sum_i \operatorname{softmax}_i\left(
\langle q_\psi(h(p,l),e(l)),s_\psi(h(p,l_i),e(l_i))\rangle/\sqrt d
\right)(r_i-\bar r).
$$

This is an amortized neural attention operator, not a ridge or closed-form parameter solve. It is
invariant to joint support permutation, changes under label-only permutation, and is exactly zero
for $k=1$. The final SAR term is the sum of endpoint change induced by $t_{S,l}$ and $c_S(l)$.

Both terms are multiplied by a support-only reliability gate. The gate compares leave-one-out
matching error on centered support residuals with the zero-residual null and maps the normalized
improvement through learned shared scalars. It never reads query labels. Thus the relative SAR
hypothesis is attenuated when the observed support set does not internally validate it, while all
gate parameters remain trained jointly by the ordinary episodic objective.
The shared matching metric additionally receives a leave-one-support-out centered-residual
reconstruction loss. This is an auxiliary loss inside the same training stage, not a separate
pretraining or adaptation phase.

Only the last declared interaction blocks are modulated. For block $m$,

$$
\Delta W_m(S)=\epsilon U_m\operatorname{diag}(\tanh(A_mt_S))V_m,
\qquad W_m^{(S)}=W_m+\Delta W_m(S),
$$

with shared low-rank bases $U_m,V_m$, shared gate $A_m$, and fixed scale $\epsilon$. Deployment
performs no target-specific gradient update, optimization, ridge regression, or linear solve. The
transient code is discarded after the episode. The active scalar predictor is

$$
f_{\theta,\psi}(p,l,a\mid S)
=f_{\theta,\Delta\Theta(t_S)}(p,l,a)+\alpha_\eta(S).
$$

Quotient preservation remains structural because every delta and rectangle is still computed from
this single scalar endpoint. Task adaptation is now amortized conditional inference, not affine
section identification. The bottom-layer strong-convexity regularizer in the separate law-valued
operator theory is not a few-shot adaptation algorithm.

This document specifies a candidate model class and the conditions under which mathematical
generalization statements apply. It does not assert that protein-specific interaction transfer has
already been established empirically.

## 2. Task and input contract

A task is indexed by a protein $p$ and a legal measurement context $a$. Its latent or observed
affinity surface is

$$
f_{p,a}:l\mapsto Y_{p,a}(l).
$$

An episode contains

$$
I=(x_p,x_l,a,S),
$$

where:

- $x_p$ is the amino-acid sequence and any structure information legally available at deployment;
- $x_l$ is the query molecular graph;
- $a$ contains the measurement type, units, and declared assay context;
- $S=\{(x_{l_i},O_i)\}_{i=1}^k$ is the support set for the same unseen protein task;
- $O_i$ is a point affinity or a closed censoring interval.

The query label is never an input. Target identifiers may be used for grouping, deduplication, and
splitting, but not as predictive features.

## 3. Encoded model inputs

The external model accepts sequence, molecular graph, context, and support. Internally it computes

$$
H_\eta(x_p)=(h_{\eta j}(x_p))_{j=1}^{m_p},
\qquad
e_\xi(x_l,a)\in\mathbb R^{d_l}.
$$

$H_\eta$ is a position-preserving protein token bank and $e_\xi$ is a ligand-context embedding.
Every support ligand is passed through the same ligand encoder as the query ligand.

Two training regimes are permitted:

1. **Frozen encoders:** $\eta$ and $\xi$ are fixed and only the QPSMP meta-learner is trained.
2. **End-to-end episodic training:** $\eta$, $\xi$, and the QPSMP parameters receive gradients from
   the source meta-objective. Their complexity must then be included in the hypothesis class.

In either regime, all representations must be computable for a new protein without target lookup,
query labels, or persistent target memory.

The current `QPSMPBioModel` code is a fixed standardized-Ki-context specialization of this general
interface. Until an explicit context encoder is implemented and audited, its claims are restricted
to that endpoint/context convention and do not include cross-assay transport.

## 4. The QPSMP scalar potential

### 4.1 Ligand-conditioned protein localization

One admissible localization map is

$$
\alpha_{\theta j}(p,l,a)
=\operatorname{softmax}_{j\in\operatorname{mask}(p)}
\left(
\frac{\langle K_\theta h_{\eta j}(p),Q_\theta e_\xi(l,a)\rangle}{\sqrt d}
\right),
$$

$$
r_\theta(p,l,a)
=\sum_j\alpha_{\theta j}(p,l,a)V_\theta h_{\eta j}(p),
\qquad
z_\theta(p,l,a)=\psi_\theta(r_\theta(p,l,a),e_\xi(l,a),a).
$$

$\psi_\theta$ is an origin-preserving crossed map with
$\psi_\theta(0,e,a)=\psi_\theta(r,0,a)=0$. A projected Hadamard product or a
bounded tensor-product map is admissible; uncrossed protein or ligand main
features may not enter the interaction heads.

The attention formula is not mandatory. A deployment-computable pocket mask, structure network, or
another permutation-compatible localizer may replace it. The retained scalar and section heads
defined below may not be replaced by an unrelated quotient head.

### 4.2 Shared baseline and interaction channels

The support-free nuisance baseline is

$$
b_{\theta,0}^{\mathrm{add}}(p,l,a)
=\mu_\theta(a)+u_\theta(p,a)+v_\theta(l,a).
$$

Level calibration and query-shape reliability are separate shared scalars. With
$f_{\theta,0}=b_{\theta,0}^{\mathrm{add}}+s_\theta$, define

$$
\alpha_\theta(S)=\omega_{\theta,k}\,\overline{Y-f_{\theta,0}}{}_S,
\qquad
\omega_{\theta,k}=\frac{k\tau_\theta^2}{\sigma_\theta^2+k\tau_\theta^2},
\qquad \sigma_\theta^2>0,\quad\tau_\theta^2>0.
$$

The current implementation learns the positive variance parameters in the meta-train standardized
label space. This is the zero-prior-mean plug-in empirical-Bayes special case of the
random-intercept model in the pure theory. The residual mean is invariant to a permutation of the
support ligand-label binding, and $\alpha_\theta(S)$ does not depend on the query ligand. Raw labels
enter the learned support operator only through the centered residual quotient defined below; query
labels never enter the support state or deployment input.

Define the retained zero-shot scalar potential and section basis by

$$
s_\theta(p,l,a)=w_\theta(a)^\top z_\theta(p,l,a),
$$

$$
\phi_\theta(p,l,a)=W_\theta(a)z_\theta(p,l,a)
\in\mathbb R^{r_{\mathrm{task}}}.
$$

For a transient support state $c$, the complete scalar predictor is

$$
f_{\theta,c}(p,l,a\mid S)
=f_{\theta,0}(p,l,a)
+\omega_{\theta,k}\overline{Y-f_{\theta,0}}{}_S
+\left[\phi_\theta(p,l,a)-\overline\phi_{\theta,S}\right]^\top c.
$$

$c$ is the only support-derived SAR/interaction state. All neural parameters are shared across
tasks.

For attribution, the implementation must expose mutually exclusive channels:
`f_zero = additive + cross_zero_shot` and
`f_few = zero_shot + calibrated_level + centered_SAR_adaptation`. The G2 interaction estimand is evaluated on
the crossed channel, while complete-scalar utility is reported separately. An
additive protein-level arm is a nuisance control, not evidence of
target-conditioned chemistry. The primary learned module must retain the
endpoint scalar path while reporting these channels separately.

## 5. Quotient-preserving outputs

The combined model-and-section interface should return

$$
\operatorname{QPSMP}_\theta
\bigl(p,l,a,T_{\mathrm{level}}(S),\widehat c_\theta(S)\bigr)
=\left(
\widehat f_{p,S}^{\mathrm{meta}}(l),
m_{p,S}^{\mathrm{sec}}(l),
\phi_\theta(p,l,a),
r_{p,S}^{\mathrm{sec}}(l)
\right).
$$

All differences are derived from retained endpoint scalar values:

$$
\Delta_\theta(p;l,l'\mid S)
=\widehat f_{p,S}^{\mathrm{meta}}(l')
-\widehat f_{p,S}^{\mathrm{meta}}(l),
$$

$$
\mathcal R_\theta((p,S),(p',S');l,l')
=\Delta_\theta(p;l,l'\mid S)
-\Delta_\theta(p';l,l'\mid S').
$$

A protein-side swap exchanges the complete task-support pairs $(p,S)$ and $(p',S')$. Both ligand
endpoints and both tasks in a rectangle must share the declared measurement context $a$ and the
same ligand-baseline parameterization; otherwise additive baseline cancellation is not asserted.

The architecture enforces:

- ligand-reversal antisymmetry of $\Delta_\theta$;
- task-support and ligand antisymmetry of $\mathcal R_\theta$;
- zero circulation of every delta field derived from the scalar potential;
- cancellation of the shared additive baseline in the rectangle;
- zero quotient prediction when either contrast is zero.

Independent delta or rectangle heads are prohibited. Quotient preprocessing may use fixed
origin-preserving linear scaling, but not a quotient intercept or sample-mean centering.

These invariants do not establish correct-protein specificity. That property requires an empirical
admission test.

## 6. Few-shot adaptation

For point-valued support, define the support-free prediction

$$
\widetilde b_i=f_{\theta,0}(p,l_i,a),
$$

the residual vector $r_S=(Y_i-\widetilde b_i)_{i=1}^k$, and

$$
\Phi_S=(\phi_\theta(p,l_i,a)^\top)_{i=1}^k.
$$

For $k\ge1$, let

$$
H_k=I_k-\frac1k\mathbf1\mathbf1^\top.
$$

### 6.1 Primary learned meta-adapter

The active primary path uses a learned row-span support-set operator. Let

$$
\widetilde r=H_k r_S,
\qquad A=H_k\Phi_S,
$$

and let a shared permutation-equivariant neural map produce evidence-bound weights
$w_\theta(A,\widetilde r)$. Define

$$
\eta=H_k\left[\widetilde r\odot\tanh w_\theta(A A^\top,\widetilde r)\right],
\qquad
\widetilde c=\frac1kA^\top\eta,
\qquad
c_\theta^{\mathrm{neural}}(S)
=R\frac{\widetilde c}{R+\|\widetilde c\|_2}.
$$

$w_\theta$ is implemented by a shared support self-attention network whose token statistics include
the centered row $A_i$, centered residual, diagonal of $AA^\top$, and row mean of $AA^\top$.
Consequently
$c_\theta^{\mathrm{neural}}(S)=0$ whenever the centered support evidence is
zero, $c_\theta^{\mathrm{neural}}(S)\in\operatorname{row}(A)$, and the adaptive
information dimension is at most $\operatorname{rank}(A)\le k-1$. The map is
permutation invariant, invariant to a constant shift of the support residuals,
bounded by $R$, and independent of query labels. Query loss backpropagates
through every unfrozen declared trainable module. For $k=1$, centered SAR
evidence is zero and only the separate level channel may adapt.

The deployed few-shot scalar is

$$
f_{\mathrm{few}}(q\mid S)
=f_{\theta,0}(q)+\alpha_\theta(S)
+\left[\phi_\theta(q)-\overline\phi_{\theta,S}\right]^\top
c_\theta^{\mathrm{neural}}(S).
$$

The zero-shot query shape is retained without an evidence gate. The structural identity
$c_\theta^{\mathrm{neural}}(S)=0$ is already the exact SAR no-op, so multiplying the SAR term by a
second evidence gate is prohibited in the primary training path. A state-norm score may be reported
as a diagnostic or calibrated later for selective prediction, but it does not alter this endpoint.
For $k=1$, the endpoint retains query-dependent zero-shot shape and only the residual level may adapt.

### 6.2 Analytic comparator and section diagnostic

For an identical learned scalar family, a future analytic transient comparator is

$$
\widehat c_\theta(S)
=\arg\min_c
\left\{
\frac1k\|H_k(r_S-\Phi_S^{\mathrm{eff}}c)\|_2^2
+\lambda\|c\|_2^2
\right\},
\qquad\lambda>0.
$$

Writing $A=H_k\Phi_S$ and $d=H_kr_S$ gives the differentiable solve

$$
\widehat c_\theta(S)
=\left(\frac1kA^\top A+\lambda I\right)^{-1}
\frac1kA^\top d.
$$

The matrix has minimum eigenvalue at least $\lambda$. A Cholesky or linear solve can therefore be
used inside the computation graph. If $A$ and $d$ are $C^1$ in the shared parameters, then the state
is $C^1$ and its implicit derivative is

$$
\frac{d\widehat c}{d\theta}
=-[\nabla_{cc}^2\mathcal L_{\mathrm{in}}]^{-1}
\nabla_{\theta c}^2\mathcal L_{\mathrm{in}},
$$

with inverse-Hessian norm at most $1/(2\lambda)$.

The identifiable SAR-state dimension satisfies

$$
\operatorname{rank}(H_k\Phi_S^{\mathrm{eff}})
\le\min\{k-1,r_{\mathrm{task}}\}.
$$

Thus one absolute support observation can calibrate level but cannot identify a centered SAR
direction. The zero-support branch is defined separately by $S=\varnothing$ and $c=0$.

Closed interval observations may be handled by a continuous convex interval loss plus positive
ridge. Uniqueness remains, but the closed-form solve generally does not. This
analytic state is a baseline, identifiability diagnostic, and potential certificate
helper. It is not the primary learned meta-adapter or the core innovation. The current
flat-feature `QPSMPCore` has independent heads and is only a legacy diagnostic; it does not
implement this identical-family bridge.

## 7. Episodic bilevel training

Training samples an independent source component, then a task, then a disjoint support-query split:

$$
C\sim\Pi_{\mathrm{src}},
\qquad
\tau\mid C,
\qquad
(S,Q)\mid\tau,C.
$$

The inner solve reads only $S$. Query regression and delta losses read only $Q$. A legal rectangle
auxiliary loss reads source rectangles whose measurement cells are excluded from the corresponding
support sets.

The population objective is

$$
\mathcal J_{\mathrm{pop}}(\theta)
=\mathbb E_C\mathbb E_{\tau\mid C}\mathbb E_{S,Q\mid\tau,C}
\mathcal L_Q(\theta,c_\theta^{\mathrm{neural}}(S))
+\lambda_{\mathrm{shape}}\mathcal L_{\mathrm{shape,pop}}(\theta)
+\lambda_\Delta\mathcal L_{\Delta,\mathrm{pop}}(\theta)
+\lambda_R\mathcal L_{R,\mathrm{pop}}(\theta),
$$

where, for $m=|Q|\ge2$,

$$
\mathcal L_{\mathrm{shape}}(\theta;Q)
=\frac1m\|H_m(Y_Q-f_{\theta,0,Q})\|_2^2.
$$

It trains level-free query shape using source query labels only. All auxiliary losses are computed
from the same retained scalar head. The empirical objective
weights independent components equally.

If end-to-end episodic training is claimed, gradients must pass from query loss through
$c_\theta^{\mathrm{neural}}(S)$ into every declared trainable encoder and head. A staged procedure that freezes
the baseline, zero-shot potential, and adapter in sequence is a phased bilevel surrogate and must be
reported as such.

The checkpoint contains only shared parameters. The state $\widehat c_p(S)$ is computed at deployment
and discarded after the episode. Trainable target-ID embeddings, per-target heads, target-specific
normalization statistics, and ID/hash/family/nearest-target label memories are prohibited.

## 8. Query sections and certificates

For a fixed learned affine family, define

$$
\mathcal C_S
=\{c:\|c\|_2\le R, f_{\theta,c}(p,l_i,a\mid S)\in O_i\ \forall i\}.
$$

If this set is nonempty, its projection at query $q$ is a compact interval. Let its midpoint and
half-width be $m_q^{\mathrm{sec}}$ and $r_q^{\mathrm{sec}}$. If it is empty, the model must abstain.

The deployed neural meta-prediction generally differs from the section midpoint. Define

$$
r_q^{\mathrm{center}}
=|\widehat f_q^{\mathrm{meta}}-m_q^{\mathrm{sec}}|.
$$

A certificate around the neural meta-prediction must satisfy, on one simultaneous event for the same
latent target and query,

$$
r_q^{\mathrm{tot}}
\ge r_q^{\mathrm{center}}+r_q^{\mathrm{sec}}
+r_q^{\mathrm{repr}}+r_q^{\mathrm{trans}}+r_q^{\mathrm{obs}}.
$$

The section term alone does not cover representation misspecification, deployment shift, or
measurement noise. Acceptance and ranking may use the total certificate only after independent
calibration of both coverage and acceptance frequency.

## 9. Cold-target generalization interface

The QPSMP hypothesis class must be fixed and measurable. Token counts, support size, inputs, outputs,
and loss are bounded; masks are nonempty; $\lambda$ has a positive lower bound; and parameters are
compact or uniformly norm constrained. Learned encoders are part of the hypothesis class.

For $N_C$ IID source components, expected Rademacher complexity $\mathfrak R_{N_C}$, bounded
component loss $M$, and empirical regression regret $\gamma_{\mathrm{reg}}$, the deployment excess
risk obeys

$$
R_{\mathrm{dep}}(\widehat h)-\inf_{h\in\mathcal H}R_{\mathrm{dep}}(h)
\le
4\mathfrak R_{N_C}
+2M\sqrt{\frac{\log(2/\delta)}{2N_C}}
+\gamma_{\mathrm{reg}}+2\Delta_{\mathcal H}
$$

with probability at least $1-\delta$, under the conditions in the pure mathematical theory.

For a deployment-computable component descriptor $D$, the shift term may be decomposed as

$$
\Delta_{\mathcal H}
\le L_{\mathrm{mech}}W_1(\Pi_{\mathrm{src}}^D,\Pi_{\mathrm{dep}}^D)
+\Delta_{\mathrm{mech}},
$$

provided the specified source conditional-risk versions are uniformly Lipschitz on the union of the
two descriptor supports. The complete bound also contains the deployment approximation gap
$\varepsilon_{\mathrm{app}}^{\mathrm{dep}}$.

Three claims must remain distinct:

1. **Architectural generality:** the same ID-free checkpoint accepts an unseen sequence and support.
2. **Conditional theoretical generalization:** the component, complexity, and transport assumptions
   imply a risk bound.
3. **Empirical generalization:** an untouched confirmation split demonstrates the required effects.

The architecture guarantees the first property and provides an interface for the second. The third
requires prospective validation.

## 10. Admission protocol

### G0: Data and estimand validity

- Pair direction, inclusion, enumeration, and truncation are outcome independent.
- Forward and reverse pairs are both materialized when no canonical direction exists.
- Support and query measurement cells are disjoint.
- Confirmation labels do not enter representation or hyperparameter selection.
- Predeclared dependency edges do not cross hard splits.

### G1: Representation coverage

- Transformation descriptors have shared numerical support across source and confirmation domains.
- Nearest-neighbor distance, leverage, effective rank, and a planted positive control are reported.
- Exact transformation identity is not used as the transferable representation.

### G2: Cold protein-conditioned quotient

Set $c=0$ and use only the retained scalar potential. The correct model must improve over all four
predeclared controls:

1. the best outcome-redacted design-nuisance model;
2. a ligand-only model;
3. a matched wrong-protein model;
4. a shuffled-protein model.

All multiplicity-adjusted lower confidence bounds must exceed predeclared effect thresholds. Both
proteins in a cold-by-cold confirmation rectangle belong to source-disjoint homology components, and
training has not read either protein's affinity labels. A source protein on one side yields only a
source-to-cold contrast.

The independent unit for inference is never the rectangle row. Dependence must be represented by the
complete preregistered closure of shared protein, measurement cell, document, transformation, and
panel links, or by a preregistered dyadic/graph bootstrap with the corresponding nested blocks.

### G3a: Few-shot HyperSAR interaction adaptation

For $k\ge2$, using the same frozen zero-shot potential, the centered support update must improve
held-out delta risk, depend on the correct support-label binding, and outperform a matched
wrong-protein state. The $k=1$ arm verifies level calibration and the zero-SAR invariant but is not
part of the interaction-adaptation admission condition.

The few-shot endpoint always retains the source-learned query shape. Quotient-null evidence and
the $k=1$ case set only the SAR term to zero; they return
$f_{\theta,0}+\alpha_\theta(S)$ rather than a constant support mean. The zero-support crossed
potential remains a separate G2 estimand.
Wrong-protein or foreign-support corruption changes only the transient code supplied to the
recipient interaction trunk. The correct query protein, query ligand, zero-shot potential, and
learned level statistic remain fixed. The protocol must retain real-versus-permuted and
real-versus-foreign support-binding contrasts.

### G3b: Scalar DTA bridge

Define the three scalar arms explicitly:

$$
f_{\mathrm{no-int}}=b^{\mathrm{add}},
\qquad
f_{\mathrm{frozen}}=b^{\mathrm{add}}+s,
\qquad
f_{\mathrm{level}}=f_{\mathrm{frozen}}+\alpha_\theta(S),
\qquad
f_{\mathrm{full}}=f_{\theta,\Delta\Theta(t_S)}+\alpha_\eta(S).
$$

All arms share the same frozen nuisance parameterization. The level and full arms use the identical
support-level statistic. The HyperSAR estimand toggles only support-conditioned interaction
modulation while retaining the same zero-shot endpoint and learned level correction. Therefore
full-versus-level is the direct interaction-adaptation utility contrast.

### G4: Certificate calibration

On independent calibration and confirmation components, validate query acceptance coverage, total
certificate coverage, and accepted-pair ranking error. Scalar regression, confidence intervals,
rank correlation, and listwise metrics remain separate estimands.

### G5: Cold transport

Report component-level statistical complexity, descriptor transport distance, conditional-mechanism
sensitivity, and the deployment approximation gap. A sequence-identity split alone is not a
transport theorem.

The required order is

```text
G0 -> G1 -> G2 -> G3a -> G3b -> G4 -> G5 -> downstream integration
```

Failure of an effect gate means that the predeclared effect was not identified under the frozen
model, data domain, sample size, and threshold. It is not a function-class impossibility theorem.

## 11. Minimal implementation contract

```text
protein_tokens = ProteinEncoder(sequence, optional_legal_structure)
query_embedding = DrugEncoder(query_molecular_graph, assay_context)
support_embeddings = DrugEncoder(support_molecular_graphs, assay_context)

level_adjustment = LearnedPermutationInvariantLevelChannel(support_residuals)
task_code = AmortizedTargetConditioner(
    support_interactions,
    support_embeddings,
    centered_support_residuals,
)
query_code = SiameseRelativeConditioner(
    task_code, query_interaction, query_embedding)
adaptive_interaction = HyperSARLowRankModulation(query_code)
y_meta = ScalarPotential(protein_tokens, query_embedding, adaptive_interaction)

delta = scalar(endpoint_plus) - scalar(endpoint_minus)
rectangle = delta(task_a) - delta(task_b)

transport_diagnostics = valid_representation_transport_observation_bounds
```

Implementation tests must verify endpoint-difference identity, four-endpoint rectangle identity,
swap antisymmetry, support permutation invariance, baseline invariance to support binding permutation,
target-ID renaming invariance, and complete exclusion of query labels from encoding, adaptation, and
model selection.

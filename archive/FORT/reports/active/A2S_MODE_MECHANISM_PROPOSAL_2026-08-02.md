# A2S meta-adaptation redirection — from transport to target state

**Pre-implementation design. No code was written for this document.**

Date: 2026-08-02
Scope: source-only. `locked` role and the A2S recipient roster sealed and not requested.
Inputs: the measured record through `A2S_TRACE_Q2_MECHANISM_DECISION_2026-08-01.md`.
Labels: every substantive statement is **FACT**, **INFERENCE** or **HYPOTHESIS**.

---

# 1. Failure diagnosis

## 1.1 What each attempt actually established

| Attempt | What it learned | Why it failed as a meta-adaptation mechanism |
|---|---|---|
| A2S-CMAL | attention over support residuals × unbounded learned scale | **FACT:** support-dependent but not beneficial. Unbounded scale on a residual aggregate (C2); attention read residuals, so no clean assignment control |
| Anchor / shrunk anchor | episode-level level offset | **FACT:** real and identifiable, and **exactly** rank-null — measured at 0.0000 CI in all 45 policy × k × stratum cells of the Q1 gate. It is calibration, not adaptation |
| Global code / SVD / FiLM / SCAO | continuous target code on frozen features | **FACT:** `τ_z ≈ 0.185`, `σ ≈ 0.997` ⇒ `ρ₅ ≈ 0.147`; the identifiability certificate fired on 0.000 of episodes. A *continuous code estimated by residual regression* is not identifiable at `k ≤ 5` |
| MDK / BIR (kernel-ridge posterior) | closed-form local posterior | **FACT:** effective dof 0.99/2.89/4.74; at k=5 indistinguishable from pooled ridge. Closed-form interpolation, not a learned adaptation rule |
| TAMSK (episode kernel router) | mixture over Tanimoto scales | **FACT (this session):** a static learned mixture *ties* one global scalar (mixture collapses to γ=1 with weight 0.91). The kernel-scale degree of freedom is not where the signal is |
| A2S-TRACE | per-pair transport reliability | **FACT:** −0.0001 CI [−0.0006, +0.0005] over the analytic bar, with a positive control recovering +0.016–0.026. Null with power |

## 1.2 The single structural fact behind all of them

**FACT.** Every mechanism in the programme's history has produced a correction of the form

```
Δ_q  =  Σ_i  w(x_q, x_i, p, ·) · r_i
```

— a weighted sum of measured support residuals. Anchor is `w = 1/k`. Ridge/KRR/MDK is `w = k_qS(K+λI)^{-1}`. CMAL is learned attention. TAMSK is a mixture of such `w`. TRACE is a learned modulation of such `w`.

**FACT (Q1, measured).** Any operator of that form is worth **nothing** unless the query is a close chemical analogue of a support compound. Fixed Tanimoto KRR minus frozen base, paired component 95 % LCB on target-macro CI:

| nearest support Tanimoto | k=3 | k=5 |
|---|---:|---:|
| `< 0.20` | −0.003 | −0.007 |
| `0.20 – 0.35` | −0.003 | −0.006 |
| `0.35 – 0.55` | −0.003 | +0.006 … +0.012 |
| `≥ 0.55` | +0.023 … +0.036 | +0.031 … +0.048 |

identical in sign across all three declared support policies.

> **INFERENCE — the diagnosis.** The programme has spent six architectures inside a function class
> whose reach is bounded by chemical distance. `Δ_q = Σ_i w_qi r_i` can only move a query that has a
> near neighbour in a 3–5 compound support set. On an unseen target with a scaffold-cold library —
> the actual A2S task — most queries have no such neighbour, so the ceiling of the whole class is
> near zero there. **This is not an estimation problem and no better `w` fixes it.**

## 1.3 What the same measurements say is still available

**FACT.** The v2 high-data target oracle (a ridge fitted on the target's own abundant labels) improves
probe ranking with LCB +0.026 / +0.070 / +0.071 at k=1/3/5. Target-specific ranking structure exists
and is recoverable — with many labels.

**FACT (this session).** Hindsight ceilings on the same probe episodes, over fixed KRR at k=5:
episode-level transport scale +0.078 … +0.107 CI; episode-level support subset +0.059 … +0.065;
per-query support subset +0.353 … +0.366 (absolute CI 0.93–0.95).

**INFERENCE.** The action class is expressive; nothing label-free predicts which action is right; and
the largest structured, *episode-level* quantity — how much to trust support evidence at all — is
worth an order of magnitude more than everything the per-pair route delivered.

**INFERENCE — the reframing this document acts on.**
The transferable object is not *how to weight compounds*. It is **what kind of target this is**, and
the query correction must be a function of the query alone, conditioned on that state:

```
Δ_q  =  R_θ( z_t , x_q ) ,    z_t = A_θ(S_t)      ← no support compound appears in the correction
```

A correction of this form is **not distance-limited**. It can act on a scaffold-cold query, which is
precisely the regime where the entire previous class measures zero.

---

# 2. Why TRACE cannot be the final contribution

**FACT.** TRACE's learned object is `m_qi = 2σ(g_ψ(φ(x_q,x_i), c_p))` — a per-*pair* reliability
multiplier on a fixed kernel row. Its state has no target semantics: there is no `z_t`. Removing the
support set does not remove a *target representation*, it removes the evidence the weights multiply.

**FACT.** Measured against the analytic bar it adds −0.0001 CI [−0.0006, +0.0005] at k=3 and
−0.00003 [−0.0006, +0.0005] at k=5, on 74–76 components, with a positive control that recovers
+0.016 to +0.026 CI in an injected world. It is a null with power, not an unmeasured possibility.

**INFERENCE.** Even a *successful* TRACE would have been a better similarity estimator. It satisfies
the letter of "learned" but not the substance of "learning how to adapt to an unseen target."

**Disposition.** TRACE is retained in exactly two subordinate roles, neither of which is a claim:

1. **Evidence encoder.** Its label-free episode statistics (LOO consistency, centred residual–kernel
   alignment, support similarity mass, support spread, Gram conditioning) are a legitimate input
   `e_t` to a target-state inference network.
2. **Mandatory local baseline.** Globally-scaled Tanimoto KRR is the bar in the `t ≥ 0.55` stratum and
   must be reported beside any new result. Any future gain claimed *inside* that stratum is measured
   against it.

---

# 3. The identifiability budget, stated honestly

**FACT.** The support supplies `k` residuals attached to `k` compounds. After the level (rank-null,
identifiable, worth 0.0000 CI) is removed, at most `k − 1` centred contrasts carry rank-relevant
information: 0 at k=1, 2 at k=3, 4 at k=5 — before noise.

**FACT.** For a *continuous* code fitted by residual regression, the per-direction reliability is
`ρ_k = τ²/(τ² + σ²/k)`, and the measured `τ/σ ≈ 0.186` gives `ρ₅ ≈ 0.147`, i.e. ≈ 0.11 bits. That
route is closed.

> **INFERENCE — the loophole the τ/σ arithmetic leaves open, and the core of this proposal.**
> `ρ_k` governs *estimating a continuous quantity whose prior dispersion the data fixes*. It does not
> govern **choosing between a small number of hypotheses**. The discriminability of two hypotheses
> `ψ_a`, `ψ_b` from `k` labels is
> `D_k = Σ_{i≤k} ( ψ_a(x_i) − ψ_b(x_i) − Δ̄ )² / σ²`,
> which depends on how differently they *predict*, not on how dispersed fitted codes are. `τ` is
> handed to us by the data; **`D_k` is an object the outer loop can maximise.** A dictionary whose
> members differ by ≈ 1 pKi RMS on typical compounds has `D_5 ≈ 5` at `σ ≈ 1` — identifiable.
>
> This is the only place in the whole record where the outer loop can change the inner problem's
> conditioning rather than merely solve it better. A2S-IDA aimed at this and attacked the wrong
> object (a continuous basis, governed by `τ`). Discrete model selection is governed by `D_k`.

**Declared budget.** Rank-relevant target-state capacity: `0` bits at k=1, `≤ 1.5` bits at k=3
(≈ 3 modes), `≤ 2.5` bits at k=5 (≈ 6 modes). Any proposal exceeding this is rejected on sight.

---

# 4. Six candidate meta-adaptation mechanisms

Each is stated as: learned object → why source targets can teach it → why `k ≤ 5` can identify it →
distinction from prior art.

### C1 — Discrete Response-Mode Selection (**A2S-MODE**)

**Learns:** a small dictionary of target *response modes* `Ψ = {ψ_0 ≡ 0, ψ_1, …, ψ_M}`, each a
correction function of the **query alone**, plus an amortised prior over modes and a selection rule.
`z_t` = posterior over `M+1` modes.
**Transferable because:** abundant source targets share a small number of recurring SAR response
patterns (e.g. potency rising with lipophilicity/size, or with polarity, or flat); a mode is a
property of the target, not of a compound pair.
**Identifiable because:** it is model selection, governed by `D_k` (§3), which the training objective
maximises. With a shrunk level, k=1 carries almost no mode evidence and k=3/5 carry 1–2 bits.
**Not:** MoE/multi-task clustering — those fit weights on abundant task data; here the weights are
inferred from `k ≤ 5` labels of an *unseen* task, under a declared bit budget, with a null expert and
a separation objective on the dictionary itself.

### C2 — Behavioural Prototype Assignment (**PROTO**)

**Learns:** `P` prototypes over *source targets* (clusters in fitted-head space), each with a
prototype head; `z_t` = posterior assignment of the recipient to prototypes.
**Transferable because:** "which known target does this behave like" is exactly what a bank of
abundant source targets can teach.
**Identifiable because:** discrete assignment among `P ≤ 8`.
**Not:** protein-similarity transfer. **FACT:** sequence/KLIFS-group transfer was falsified in this
programme (own-group-cold G0 fails; recorded as *TR group not resolvable*). PROTO assigns by
**measured support behaviour**, not by sequence — a different claim.
**Relation to C1:** C2 is C1 with prototype-initialised, source-target-derived modes. They share one
inference machine and should be one experiment, not two.

### C3 — Conditional Ranking-Rule Policy (**RULE**)

**Learns:** a library of interpretable monotone ranking rules on ligand axes (size, logP, HBD/HBA,
aromatic rings, TPSA) plus a policy that selects/weights them from support evidence.
**Identifiable:** discrete choice among ≈ 5 rules.
**Weakness — INFERENCE:** it is C1 with a hand-fixed, rank-1-per-axis dictionary. Strictly less
capacity, but a valuable *control*: if RULE ties C1, the learned dictionary is not load-bearing.
**Disposition:** adopt as a mandatory nested control, not as the mechanism.

### C4 — Sparse SAR-Transformation Dictionary (**TRANSFORM**)

**Learns:** a dictionary of substructural transformation responses ("this exchange raises potency");
`z_t` = sparse set of atoms active for this target.
**Transferable because:** matched-pair SAR is the medicinal-chemistry unit of transferable knowledge.
**Identifiability risk — FACT:** it needs the query to be reachable from a support compound by a
dictionary transformation, so it inherits the *same distance limitation* that §1.2 identified as the
class-level failure. Measured nearest support–query Tanimoto in the passive construction is ≈ 0.223.
**Disposition:** reject for now on the measured geometry; revisit only inside the `t ≥ 0.55` stratum
where it would compete with a baseline that already wins there.

### C5 — Sequential Target-State Filter (**STATE**)

**Learns:** an amortised recursive update `z ← U_θ(z, (x_i, r_i))` over a low-dimensional learned
state manifold; prediction `f_0 + R_θ(z, x_q)`.
**Weakness — INFERENCE:** with a *continuous* `z` this is governed by `τ`, i.e. the closed route, and
it is a renamed CNP/MetaFun with a declared dimension. It becomes interesting only if the manifold is
discrete or heavily quantised — at which point it is C1 with extra machinery.
**Disposition:** reject as primary; the recursive form adds no identifiability.

### C6 — Risk-Calibrated Adaptation Program (**PROGRAM**)

**Learns:** a tiny program `z_t = (mode m_t, trust level, abstain flag)` trained with a
harm-asymmetric loss.
**Disposition:** not an independent mechanism — it is C1 plus the abstention layer. Adopt the null
mode and the harm-asymmetric term **into** C1. It is how criterion 9 (negative-transfer rate) is met.

### Screen

| | learned object is a target state | identifiable at k≤5 | acts on scaffold-cold queries | load-bearing in main path | structurally abstains | verdict |
|---|---|---|---|---|---|---|
| **C1 MODE** | ✅ posterior over modes | ✅ separation-driven | ✅ correction ignores support chemistry | ✅ | ✅ null mode + posterior-minus-prior | **SELECT** |
| C2 PROTO | ✅ | ✅ | ✅ | ✅ | ✅ | **merge into C1** as initialisation + ablation |
| C3 RULE | ✅ | ✅ | ✅ | ✅ | ✅ | **mandatory control** |
| C4 TRANSFORM | ✅ | ⚠ | ❌ distance-limited | ✅ | ✅ | defer |
| C5 STATE | ⚠ continuous | ❌ τ-governed | ✅ | ✅ | ⚠ | reject |
| C6 PROGRAM | — | — | — | — | ✅ | absorb into C1 |

---

# 5. A2S-MODE — formulation

## 5.1 Objects

- `f_0(p, x)` — frozen support-free base (unchanged, cross-fitted, already built).
- `g(x) ∈ R^d`, `d ≤ 32` — compact **label-free** ligand basis (fixed descriptors + a small learned
  projection). Capacity is declared, not tuned.
- **Mode dictionary** `Ψ = {ψ_0, …, ψ_M}` with `ψ_0 ≡ 0` (the null mode) and
  `ψ_m(x) = ⟨u_m, g(x)⟩`, `u_m ∈ R^d`, `M ≤ 5`.
- **Prior network** `π_θ(e_t) ∈ Δ^{M}` — label-free: protein embedding, support chemistry summary,
  `k`, and (optionally) TRACE's episode statistics. It never sees a label.
- `σ²` — meta-learned residual scale; `a` — level, with a shrinkage prior `a ~ N(0, s²)`.

## 5.2 Adaptation at meta-test (one pass, no fitting)

Support residuals `r_i = y_i − f_0(p, x_i)`. For each mode `m`, profile the level under its prior:

```
â_m   = ( s² / (s² + σ²/k) ) · mean_i [ r_i − ψ_m(x_i) ]                       (empirical Bayes level)
ℓ_m   = − (1/2σ²) Σ_i ( r_i − â_m − ψ_m(x_i) )²  − â_m²/(2s²)                  (mode log-evidence)
z_t   = softmax_m [ β · ℓ_m + log π_θ,m(e_t) ]                                  ← THE TARGET STATE
```

`z_t ∈ Δ^{M}` is the entire target-specific object. **Dimension: `M ≤ 5`, entropy budget ≤ 2.5 bits.**

## 5.3 Prediction

```
ŷ_q = f_0(p, x_q)  +  â                                   ← LEVEL channel (rank-null, RMSE only)
                   +  Σ_m ( z_{t,m} − π_θ,m(e_t) ) ψ_m(x_q)   ← RANK channel: posterior MINUS prior
```

**The correction contains no support compound.** It is a function of `x_q` and the inferred state.

Structural consequences, provable from the form:

| property | why |
|---|---|
| **Support removed ⇒ exact no-op** | no support ⇒ `z_t = π_θ` ⇒ the rank term is identically 0. The claim is literally *what the k labels taught* |
| **k=1 rank action ≈ 0** | with one point and a level, all modes fit nearly equally; `ℓ_m` differences collapse. Any material k=1 rank gain is a leakage alarm |
| **Bounded** | `\|Δ_q\| ≤ 2·max_m \|ψ_m(x_q)\|`, a quantity fixed at training time; no learned scale multiplies a residual aggregate (C2 excluded) |
| **Level/rank separation** | `â` is an episode constant on its own reported channel |
| **Not distance-limited** | `ψ_m(x_q)` is defined for every query, analogue or not |

**Deliberate departure, stated up front — INFERENCE.** Under a likelihood formulation `r_S ≡ 0` is
*informative*: it says the base is already right, so the posterior should move to `ψ_0`. TRACE's
bitwise `r ≡ 0 ⇒ Δ ≡ 0` guarantee is therefore **replaced** by two others: *support removal* is the
exact no-op, and residual-null is a **registered prediction** (`z_{t,0} → 1`), not a control.

## 5.4 Training objective (source episodes only)

```
L =  w_rank · L_smoothCI( ŷ, y ; within-episode pairs )        ← bounded surrogate, see below
   + w_level · L_level( â )
   + w_sep  · L_separation
   + w_harm · E[ relu( L_rank(adapted) − L_rank(base) ) ]      ← harm-asymmetric, criterion 9
   + w_ent  · E[ H(z_t) ]  (annealed)                          ← discourages permanently diffuse states
```

**The separation term is the innovation carrier**, the discrete analogue of what A2S-IDA tried to do
for a continuous basis:

```
L_separation = − E_{t, S ~ k} [  min_{m ≠ m'}  Σ_{i≤k} ( ψ_m(x_i) − ψ_{m'}(x_i) − Δ̄_{mm'} )² / σ²  ]
```

It maximises the *worst-case* pairwise discriminability of the dictionary **on a random k-subset**,
i.e. it shapes the dictionary so that `k` labels can actually tell its members apart. `D_k` is
reported per run alongside realised mode-selection accuracy.

**FACT (measured this session, carried forward).** The ranking term must use the **bounded
smoothed-CI surrogate** `σ(−Δŷ/τ)`, not the convex RankNet logistic. On this substrate the convex
surrogate's optimum sits far below the CI optimum, and a model trained on it will not find the right
correction magnitude.

## 5.5 Where the information comes from — and the falsifiable difference from the transport class

The declared pipeline is

```
chemical evidence e_t (label-free)  +  k support residuals r_S
                    │
                    ▼
        transferable adaptation state  z_t ∈ Δ^M      (≤ 2.5 bits, §3)
                    │
                    ▼
   small target-specific intervention  Δ_q = Σ_m (z_{t,m} − π_{θ,m}) ψ_m(x_q)
```

with similarity/reliability models (TRACE, atom-pair KRR, metric learning) entering **only** at the
first arrow, as encoders of `e_t`, and standing elsewhere **only** as baselines.

**FACT (Q1).** The *transport* class `Δ_q = Σ_i w_qi r_i` is admitted only when the query is close to
a support compound. Its information requirement is **support→query proximity**.

**INFERENCE — the structural difference, and it is testable.** A2S-MODE's information requirement is
different in kind. The mode log-evidence

```
ℓ_m − ℓ_{m'}  ∝  Σ_{i≤k} ( ψ_m(x_i) − ψ_{m'}(x_i) − Δ̄ ) · ( … )
```

is large when the **support set itself spans the axes that separate the modes**. It says nothing
about where the queries are. So:

> **Registered prediction P0.** MODE's realised gain correlates with **support diversity on `g`**
> (the centred spread of `ψ_m(x_i)` across the support) and is **approximately flat in
> support→query Tanimoto**. The transport class shows the opposite pattern — flat in support
> diversity, steeply increasing in support→query Tanimoto.
>
> If MODE's gain instead tracks support→query proximity, it is transport in disguise and the claim is
> withdrawn. This single 2-D readout (gain vs support diversity × support–query Tanimoto) separates
> the two mechanism families cleanly and costs one extra column in the existing records table.

This is why Gate A3 (§7) is the decisive gate rather than a nice-to-have: it asks for a gain exactly
where the transport class measures zero.

**Consequence for episode construction — INFERENCE.** If P0 holds, the *support policy* matters for a
new reason. Q1 varied the policy to change support→query proximity; MODE cares about support
**spread**. A diversity-aware support policy is therefore a legitimate, declarable design axis for
MODE (and a confound to control: any policy comparison must hold `k` and the query set fixed).

## 5.6 Comparison with prior art, and the exact increment

| Method | What adapts | Why A2S-MODE differs |
|---|---|---|
| MAML | all weights, by gradient steps on `k` points | continuous, high-dimensional, `τ`-governed, no abstention. MODE replaces gradient adaptation with **discrete model selection over a meta-learned, separation-shaped dictionary** |
| ANIL | last layer only, by gradient | same objection at lower dimension; still a continuous `d`-vector from `k ≤ 5` points |
| CNP / ANP / MetaDTA | a continuous context vector from an `(x,y)` set encoder | no declared state dimension, no identifiability statement possible, and — **FACT, measured in this programme** — a plain CNP's apparent gain survived residual-null, wrong-target and deranged support |
| AdaMBind | task-adaptive continuous update | same identifiability burden; no declared bit budget |
| Deep Kernel Transfer / ADKF-IFT | meta-learn a representation so a **target-specific GP solve** works | the estimand is target-specific and the action is similarity-driven ⇒ inherits the measured distance limit of §1.2. MODE has **no similarity term in the correction** |
| GP / KRR / ridge | closed-form interpolation of support labels | rank-null outside `t ≥ 0.55` (Q1). MODE's action is defined everywhere |
| Retrieval | nearest-neighbour transfer | same limit, plus candidate-set dependence |
| Fine-tuning | many-shot weight update | needs ≫ 5 labels; also violates the no-recipient-SGD contract |
| Mixture-of-experts / multi-task clustering | expert weights fitted on abundant task data | MODE infers the mixture from `k ≤ 5` labels of an **unseen** task, under a declared bit budget, with a null expert, and **trains the experts to be separable at that budget** |
| A2S-TRACE | per-pair transport reliability | no target state; distance-limited; measured null |

> **The increment in one sentence.** Previous work asked *how strongly should this support compound
> speak about this query*; A2S-MODE asks *which of a few learned target behaviours is this target
> exhibiting*, and it shapes those behaviours so that `k ≤ 5` labels can tell them apart.

---

# 6. Required definitions (deliverable item 7)

**Core innovation.** A meta-learned, **separation-shaped discrete response-mode dictionary** with an
amortised label-free prior, where the target state is the posterior-minus-prior over modes and the
query correction is similarity-free.

**Trainable module (in the main prediction path).** `{u_m}` (dictionary), the projection inside `g`,
`π_θ` (prior net), `β`, `σ`, `s`. The frozen base is not retrained.

**Information source.** `k` support residuals against the frozen base (the only label channel);
label-free episode evidence for the prior; abundant source targets for the dictionary.

**Training objective.** §5.4 — bounded smoothed-CI ranking + separated level + **separation
regulariser** + harm-asymmetric penalty + annealed entropy.

**Ablations (mapping to the user's Requirement 2 and criteria 1–9).**

| # | Ablation | Expected if the claim is true |
|---|---|---|
| A1 | remove adaptation module (`z_t ≡ π_θ`) | recovers the frozen base **exactly**; gain → 0 |
| A2 | replace with closed-form KRR / ridge / **globally-scaled KRR** / TRACE | MODE ≥ these **outside** `t ≥ 0.55`; inside it, MODE is complementary, not necessarily better |
| A3 | random support (drawn from a different target's compounds) | gain → 0 |
| A4 | wrong-target support, norm-matched | gain → 0 |
| A5 | label permutation / magnitude-matched residual derangement | gain → 0 (posterior scrambles) |
| A6 | support removal (k = 0) | **structurally exact** no-op |
| A7 | `M = 1` (single global head) | isolates "modes" from "a better global base" — the key confound |
| A8 | null mode removed | negative-transfer rate rises; total gain may not |
| A9 | separation term removed | mode-selection accuracy falls toward chance; this is the innovation's own falsifier |
| A10 | RULE control (fixed interpretable dictionary) | if it ties, the learned dictionary is not load-bearing |
| A11 | protein zero / protein shuffle in `π_θ` | sizes the protein claim honestly; if ≈ 0, the method is not called protein-conditioned |
| A12 | prototype-shuffled state (assign a random mode) | gain → 0 |

**Reported alongside:** target-macro CI and NDCG@10 with paired component bootstrap (≥ 2,000 draws,
homology component as the unit), RMSE/MAE, risk–coverage, **conditional harm rate**, negative-transfer
rate, per-mode usage, realised `D_k`, mode-selection accuracy vs chance, and every number split by
`k ∈ {1,3,5}` **and** by relation stratum.

---

# 7. Pre-registered gates that must run **before** implementation

These are measurements on the existing substrate, not the mechanism. Each can kill the route cheaply.
The programme's own discipline (`no mechanism may be trained before the admissible object is
measured`) applies here exactly as it did to Q1.

**Gate A0 — is there anything to select?**
Fit a per-target head on the compact basis `g` using each probe target's *abundant* labels; measure
its ranking gain over the frozen base, **per relation stratum**. If the gain is confined to
`t ≥ 0.55`, the whole premise of §1.3 is wrong and MODE offers nothing new.
*Prior expectation (v2 high-data oracle LCB +0.026/+0.070/+0.071) says this passes; it has never been
measured stratum-by-stratum.*

**Gate A1 — mode sufficiency.**
Cluster source-target heads into `M ∈ {2,3,4,6}`; on **held-out** probe targets measure the gain of
the *oracle best mode* against (i) the frozen base and (ii) `M = 1`. Pass requires the `M > 1`
oracle to beat `M = 1` with a paired component LCB above the 0.005 MDE.
*This is the confound that killed the group-transfer route before: a "mode" that is really just a
better global head.*

**Gate A2 — k-shot identifiability (decisive).**
With the modes from A1 frozen, select by `k`-shot evidence and measure (a) selection accuracy vs
chance `1/M` at k = 1/3/5, and (b) realised gain of k-shot selection as a fraction of the oracle-mode
gain. Pass requires accuracy materially above chance at k = 3 and 5, and a realised gain LCB above
MDE. **If A2 fails, A2S-MODE is dead and the honest deliverable is the bit-budget bound.**

**Gate A3 — complementarity (the reason to prefer this route).**
Restrict to the `t < 0.35` strata where the entire transport class measures zero. Pass requires a
positive k-shot MODE gain there. This is the only place a new mechanism can add something the
existing analytic bar does not already have. Report the P0 readout at the same time: realised gain
as a 2-D function of **support diversity on `g`** × **support→query Tanimoto**. A gain that rises
with proximity rather than with diversity is transport in disguise and fails the gate even if the
headline number is positive.

**Gate A4 — synthetic positive control.**
As in Q2: inject a world with true discrete modes and verify the pipeline recovers a known fraction.
**No null may be reported from a gate that has not passed A4.**

**Stop rule.** A1 and A2 both pass → implement A2S-MODE and run the ablation ladder. A1 or A2 fails →
report the measured mode-sufficiency ceiling and the `k`-shot bit budget as the deliverable, and do
not build a seventh architecture.

---

# 8. Maximum scientific risk, stated plainly

1. **HYPOTHESIS at risk.** That source targets partition into a *small* number of *separable*
   response modes at all. If per-target heads form one diffuse cloud, A1 fails and no discrete state
   exists to infer. This is the single most likely failure and it is testable in one pass.
2. **The `M = 1` confound.** A mode dictionary can improve results purely by being a better global
   ligand model. A7/A1(ii) exist to catch exactly this, and the claim dies if it is not caught.
3. **Power.** 74–76 components in the admitted stratum, ~110 overall; MDE80 ≈ 0.005 CI. A true effect
   below that is undetectable here regardless of mechanism quality, and the deliverable would then be
   the required sample size.
4. **Criterion 7 conflict — flagged for your decision.** The requirement asks for improvement at
   k=1. Under this design the k=1 **rank** channel is structurally near-silent (§5.3); k=1 can only
   improve **RMSE**, through the level channel, which is already a replicated result. I consider the
   k=1 rank silence a feature (it is a free leakage detector), but it means criterion 7 can be met at
   k=3 and k=5 only. Relaxing it requires weakening the level shrinkage, which trades the detector
   away.
5. **Honest scope.** Even a full pass gives *stratum-and-policy-qualified* evidence on source probe
   components. Confirmation still requires freezing the protocol and opening `locked` once, and the
   recipient roster stays sealed after that.

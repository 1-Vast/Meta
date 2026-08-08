# The Differential Reference Panel: a novel `z_bio` module

Status: **PROPOSAL. Not a registered stage, not an artifact, not hashed.**
Third companion to `IDENTIFICATION_ROADMAP_AND_Z_ADMISSION.md` and
`SOLUTION_MENU_LITERATURE_INFORMED.md`. 2026-08-08.

Derived from a direct reading of `theory/FINAL_FROZEN_THEORY/`. Changes no
frozen object. Nothing executed; no label read.

---

## Part I — Three things the frozen theory says that the project has not used

### I.1 The theory explicitly disclaims the functional that was measured

`chapters/03_SCOPE_AND_REFERENCE/07_SCOPE_LIMITATIONS.md`, "Not claimed":

> pairwise, listwise, or metric ranking; coherent joint-order learning;
> derivation of ranking from affinity regression

Eight stages — P1C, P1R1, P1R2A, P1R2B0, P1R2B1, P0, H0A, H0C — used within-task
concordance as the primary readout. The frozen theory states in its own scope
chapter that it provides no theorem for that quantity. E-AFF-R0 later proved
empirically that the metric is exactly blind to the task-level channel. **The
theory said so first.** This is not a reinterpretation of any verdict; it is an
observation that the readout was outside the guaranteed scope from the beginning.

### I.2 `z` must be *sufficient*, not merely informative

From `01_FOUNDATIONS.md`, the conditional base risk is

```text
L_0(z, beta) = E[ L(beta, Y) | zeta = z ]
```

and (S-CONT) requires this to have an everywhere-defined version continuous in
`z`. The learning target `g*_mu(z) = argmin_p J_mu(z,p)` depends on `z` **only
through `L_0(z,·)`**.

The admission criterion for a biological statistic is therefore sharper than the
project's current phrasing. It is not "does `z_bio` correlate with affinity" or
"does it beat a baseline". It is: **does conditioning on `z_bio` pin down the
conditional law of `Y`?** Two tasks with the same `z_bio` must have the same
affinity law. That is a conditional-distribution requirement, and it is the
reason a location/law-valued readout is the theory-conformant one.

### I.3 The theory contains a computable dimension budget for `z`, and nobody has computed it

From `04_META_LEARNING_FORMULATION.md`, `D_N = (m+1) nu_N` where `nu_N` is the
mesh node count; at resolution `r` in `d_z` dimensions, `nu_N = (r+1)^{d_z}`.
The consistency schedule in `06_CALIBRATION_AND_GENERALIZATION.md` requires

```text
D_N log(Lambda N) / N  ->  0.
```

With `m = 7` and the D0/D1 corpus at `N ≈ 3,817` governed tasks, demanding
`D_N ≲ N` gives

```text
8 * (r+1)^{d_z}  ≲  3,817     =>     (r+1)^{d_z} ≲ 477
     r = 1  =>  d_z <=  8
     r = 2  =>  d_z <=  5
```

**The frozen theory's own rate says `z` should be roughly five to eight
dimensional at the available data scale.** Two consequences follow immediately:

- `config.py`'s `d_z = 28` is far outside the budget. It is documented as a
  declared engineering choice and a placeholder, and it should be re-registered
  downward rather than treated as a target to fill.
- The seven-dimensional `z_bio` used in L0 and L0R was, by this criterion,
  **exactly the right size**. That choice was never justified from the theory;
  it should have been, and it is now.

Any proposed `z_bio` that cannot fit in ~5–8 coordinates is inconsistent with
the theory's own generalization guarantee. This is a hard design constraint, and
it is the single most useful thing the theory folder offers a module designer.

---

## Part II — The diagnosis: the project has been attempting the harder of two problems

Two distinct questions have been conflated throughout the ledger.

| | Question | Instrument | Status |
|---|---|---|---|
| **Q1** | Can `delta[t,l]` be **measured** for a new target from its `k` support labels? | support labels | **never asked** |
| **Q2** | Can `delta[t,l]` be **predicted** for a new target from protein features alone? | ESM + geometry | asked ~20 times, always null |

Every stage from P1C to L0R attacked Q2. But MetaSieve is a **few-shot system**:
its deployment contract gives `k = 5` measured affinities for the new target, and
the frozen theory is explicitly *support-conditioned* — `z = z(S,Q,gamma)`, with
`S` in the statistic.

Inductive matrix completion theory says why the ordering matters. Recovering a
low-rank interaction from side information has sample complexity scaling with
the **intrinsic feature dimension rather than the ambient matrix dimension** —
but only when the latent factors lie in (or near) the span of the supplied
features. That span condition is the precondition for Q2, and **it has never
been tested**, because testing it requires knowing the latent factors, which
requires solving Q1 first.

Q1 is also strictly easier: measuring beats predicting. And Q1 is the
product-relevant question for a few-shot system.

**The proposal below is an instrument for Q1**, from which Q2 becomes a plain
regression with a directly observed target instead of a 1.67% variance component
chased end-to-end through a frozen operator.

---

## Part III — The Differential Reference Panel (DRP)

### III.1 The construction

Fix, label-blind and before scoring, a **reference panel** `R` of source targets
and the source interaction estimate on them.

**Offline, from source data only.**

1. Fit main effects on the source corpus: `mu[s]` per assay stratum (the `kappa`
   channel), `alpha[t]`, `beta[l]`. Cross-fitted by closure component.
2. Form the source interaction residual `delta_hat[t,l]` on the panel targets
   times shared ligands.
3. Take its leading `r` singular directions. This is the multiplicative half of
   an AMMI decomposition, or equivalently a factor-analytic FA(`r`) fit:
   ```text
   delta_hat[t,l]  ~  sum_{q=1..r}  lambda_q  u_q(t)  v_q(l).
   ```
   Freeze the ligand loadings `v_1..v_r` and hash them.

**Online, for a new target `t` with support ligands `L_S`, `|L_S| = k`.**

4. Form all `C(k,2)` within-target support differences
   `d_t(i,j) = y(t,l_i) - y(t,l_j)`.
5. Subtract the reference panel's mean difference on the same ligand pair,
   `d_R(i,j)`. The result is a **double difference**:
   ```text
   DD_t(i,j) = [y(t,l_i) - y(t,l_j)] - [ y_R(l_i) - y_R(l_j) ].
   ```
6. Project the `DD_t` vector onto the frozen `v_q` restricted to `L_S`, giving
   **relative interaction coordinates** `u_hat_1(t) .. u_hat_r(t)` in closed form.
7. Emit
   ```text
   z_bio  =  [ tau_hat ,  u_hat_1 .. u_hat_r ,  c ]
   ```
   with `tau_hat` the level coordinate (support mean relative to `b_pop[kappa]`)
   and `c` the identifiability coordinate defined in III.4.

Prediction for a query ligand is closed-form:
`delta_hat(t, l_q) = sum_q u_hat_q(t) v_q(l_q)`.

### III.2 Why the double difference is the right object

Under the model `y[s,t,l] = mu[s] + alpha[t] + beta[l] + delta[t,l] + eps`, the
step-5 subtraction cancels `mu[s]` (same stratum), `alpha[t]`, `alpha` of every
reference, and both `beta[l_i]` and `beta[l_j]`. **Every additive term is gone
by construction.** What remains is pure interaction contrast plus noise.

This matters because the ledger's dominant failure mode is main effects
masquerading as mechanism: `77.58%` ligand main, `20.76%` protein main, `1.67%`
interaction (P1R2A). A statistic that is *algebraically* free of the first two
cannot be fooled by them. Compare H0C, where a support-fitted nuisance was added
to both arms and the derangement contrast was structurally voided.

### III.3 The information amplification

This is the quantitative case, and it is the reason to expect a different
outcome from the same corpus.

F-01 treated the support as five noisy point estimates of a target's level, and
found `k=5` support means estimated level at RMSE `0.363`–`0.389`, no better
than the between-target spread. That is the correct verdict for that use.

The same five labels yield `C(5,2) = 10` within-target differences, each
compared against a reference panel of `R` targets. The support is not five
scalars; it is a **ten-dimensional differential fingerprint measured against a
frozen panel**.

Noise accounting, using the project's own `sigma_assay = 0.47971`:

```text
Var(d_t(i,j))        = 2 sigma^2                 SD = 0.678
Var(d_R(i,j))        ≈ 2 sigma^2 / R             negligible for R >~ 50
Var(DD_t(i,j))       ≈ 2 sigma^2 (1 + 1/R)       SD ≈ 0.68
```

Against X0's registered interaction size (interaction RMS `= 0.5 sigma`), the
signal in one pair difference has SD `≈ sqrt(2) * 0.24 = 0.34`, so

```text
per-pair SNR   ≈ 0.34 / 0.68  =  0.5
aggregate SNR  ≈ 0.5 * sqrt(C(k,2))
```

| `k` | pairs | aggregate SNR |
|---:|---:|---:|
| 5 | 10 | **1.6** |
| 10 | 45 | 3.4 |
| 20 | 190 | **6.9** |
| 40 | 780 | 14.0 |

Two conclusions, both actionable. At the production contract `k = 5` the design
is marginal and should not be the first test. At `k = 20` — **which is exactly
the support size H0A, H0C, L0 and L0R already used** — the SNR is roughly 7, and
that is a regime in which a real effect should be plainly visible. The first DRP
experiment should therefore run at `k = 20`, where the corpus already supports
it, and treat `k = 5` viability as a separate downstream question.

### III.4 Abstention becomes computable rather than declarative

`THEORY_BIOLOGY_INTEGRATION.md` §4.3 lists "support coverage — projection of the
query onto the support row space" as coordinate 7 and correctly flags it as an
identifiability coordinate rather than a chemistry one. Under DRP it acquires a
formula.

Let `V_S` be the frozen ligand-loading matrix restricted to the support ligands.
Step 6 is a least-squares projection onto `V_S`, so define

```text
c  =  sigma_min(V_S) / sigma_max(V_S)      (bounded to [0,1])
```

If the support ligands fail to excite the interaction directions, `V_S` is
ill-conditioned, `u_hat` is unidentified, `c -> 0`, and the correct output is
mass on `beta_0` — the population band for that assay stratum. Because the
operator emits a **set** of laws (`02_TARGET_AND_OPERATOR.md` §1), a wide band is
a valid honest output, not a failure. Abstention is a point of `Delta_m`.

This also makes support-set *selection* a first-class design variable: choose
support ligands to maximize `sigma_min(V_S)`, which is a classical D- or
E-optimal design problem. Nothing in the ledger has ever chosen supports on
purpose.

### III.5 Every frozen constraint is satisfied, and two open problems close

| Frozen requirement | How DRP meets it |
|---|---|
| bounded | coordinates squashed to `[0,1]`; `c` bounded by construction |
| finite-dimensional, `d_z ≈ 5–8` (§I.3) | `1 + r + 1` with `r = 2–3` ⇒ `d_z = 4–5` |
| deterministic | closed-form projection, no sampling |
| permutation-invariant in `S` | `u_hat` is a least-squares fit on the unordered `DD` multiset |
| free of query labels | uses support labels only, which the theory admits via `z(S,Q,gamma)` |
| `d_adapt <= k` | `d_adapt = r + 1 ≈ 3–4` against `k = 5` — **the `m = 7 > k = 5` problem dissolves**, because adaptation happens in the coordinate estimate, not in `Delta_7` |
| sign of the affinity increment | fixed by the ordered anchor ladder (dominance gap `0.0`, verified), not estimated |

Two things worth noting. First, **there is no inner-loop optimizer.** Steps 4–7
are linear algebra. The failure mode that consumed E0, E0S, E0R0 and E0R1 —
a non-converged gradient descent reported as a scientific negative — is
structurally impossible here. Second, the sufficiency argument of §I.2 becomes
statable: if `delta` is rank-`r`, then `u_hat(t)` identifies the target's
interaction behaviour up to the panel span, so `z_bio` is sufficient for the
interaction component of the conditional law. That is the form of claim the
theory actually requires.

### III.6 Intellectual lineage, and what is new

**Relative representations** (Moschella et al., ICLR 2023 oral) propose "the
latent similarity between each sample and a fixed set of anchors as an
alternative data representation", enforcing invariance "without any additional
training". DRP is that idea moved out of embedding space into **affinity space**,
with one strengthening: the relation used is a double difference, so the
representation is invariant to additive target and ligand effects *algebraically*
rather than approximately.

**Connectivity Map** is the biological precedent — characterize a perturbation by
its differential signature against a reference panel, and match by signature
similarity. DRP is the affinity-space analogue for targets rather than
perturbations.

**Difference-in-differences** supplies the identification language: the reference
panel is the control group, additivity is the parallel-trends assumption, and
`delta != 0` is precisely a violation of it.

**AMMI / factor-analytic mixed models** supply the low-rank structure and the
empirical prior that one to three interaction components dominate — which is what
makes `r = 2–3`, and hence the theory's dimension budget, plausible.

**What is new, as far as I could establish.** Kinase selectivity profiling
characterizes *compounds* against target panels. DRP inverts this: it
characterizes a *target* by its differential response against a panel of
reference *targets*, computed through shared ligands, and uses that as a few-shot
task statistic. I searched for a direct precedent and did not find one — but
absence of search hits is weak evidence, and a proper novelty check is required
before any claim is made in writing.

There is also a specific relationship to `arXiv:2510.14419`. Its Proposition 2
shows that a drug- and target-permutation-equivariant algorithm scores *exactly*
`0.5` on interaction for off-training-set targets, and offers side information as
the remedy. DRP is a **second, complementary remedy**: it breaks the permutation
symmetry using the new target's own support labels rather than its features.
That route sits outside the paper's IDIT/IDOT/ODIT/ODOT partition, which is
inductive; DRP is transductive/few-shot. Whether the proposition extends to that
setting is a question worth putting to those authors.

---

## Part IV — Second module: Mondrian-conformal band calibration

The frozen operator emits a credal set `K(beta)` with an **asymptotic** guarantee
and a fixed `2h` floor, under (S-IID) exchangeability. It carries no
finite-sample coverage statement.

Conformal prediction supplies exactly that, under the exchangeability assumption
the theory already makes. Fisch et al. (ICML 2021), "Few-shot Conformal
Prediction with Auxiliary Tasks", casts conformal prediction "as a meta-learning
paradigm over exchangeable collections of auxiliary tasks", obtains tighter sets
while keeping marginal guarantees, and validates on computational chemistry for
drug discovery. The structural match to the frozen formulation — exchangeable
tasks, set-valued output, few-shot — is close.

Proposal: calibrate `b_pop` and the emitted band **Mondrian-conformally,
stratified by `kappa`**. This does three things at once:

1. gives the band per-stratum finite-sample validity;
2. activates the `kappa` channel, which `THEORY_BIOLOGY_INTEGRATION.md` §4.2
   notes has never been used and which `run_eaff_l0r.py:219` bypasses with a
   single pooled band;
3. makes abstention *calibrated* — a wide band then means a quantified "I don't
   know" rather than a declared one.

None of this touches `V`, `h`, `Delta_m`, `mu`, the band polytope or the
operator. It changes how `b_pop` is estimated, which the theory lists as a
deployment component.

---

## Part V — Staged experiments

| Stage | Question | Controls | Freeze | Cost |
|---|---|---|---|---|
| **DRP-0** | Does the corpus support panel double differences? Label-blind census of ligand overlap between candidate targets and a reference panel. | none needed | [E] | days |
| **DRP-A** | Does *measured* `u_hat` improve prediction at `k = 20`? | support-label permutation (F-09 style); reference-panel shuffle; coupling-null | [R] | one panel |
| **DRP-B** | Can protein features **predict** `u_hat`? Plain regression of `u_hat(t)` on ESM + P1B geometry, reporting `R^2` with closure-component bootstrap. | protein derangement | [R] | cheap once A exists |
| **DRP-C** | Full few-shot Gate under the frozen operator with conformal calibration, at `k = 20` then `k = 5`. | the full five-arm ladder | [R] | one panel |

**DRP-0 inverts a known confound into an enabling condition.** The F-11
root-cause assessment recorded that "both benchmarks have complete ligand overlap
across target splits", treating it as an advantage that contaminated the ligand
baseline. Complete ligand overlap is exactly what a reference-panel double
difference needs. What made absolute prediction easy to fake makes differential
measurement possible. On this criterion the ranking is: kinase profiling panels
(complete overlap by construction) > KIBA/DAVIS > ChEMBL (sparse overlap) — which
independently reproduces the corpus recommendation from the solution menu, by a
different argument.

**Licensing.** DRP-A licenses few-shot claims only. DRP-B is the inductive
matrix completion span test and is what licenses any zero-shot claim. Neither
alone authorizes `z` admission; the Z0–Z4 chain in the roadmap still applies in
full.

---

## Part VI — Honest risks

1. **DRP creates no information.** If `Var(delta) = 0`, `u_hat` is pure noise and
   the module will faithfully report nothing. Its virtue is that it *measures*
   rather than *predicts* `delta`, making it the most sensitive available test of
   `Var(delta) > 0` — informative even when it fails.
2. **`k = 5` is marginal** at SNR ≈ 1.6 (§III.3). Do not lead with it. If DRP
   works at `k = 20` and fails at `k = 5`, the honest conclusion is that the
   production contract needs more supports, which is a product decision rather
   than a scientific failure.
3. **Reference-panel choice is gameable.** The panel, the rank `r`, and the
   ligand loadings `v_q` must be frozen and hashed before any scoring, selected
   label-blind or on source folds disjoint from evaluation.
4. **Stratum matching is required.** The `mu[s]` cancellation in §III.2 assumes
   the new target and the panel are compared within one assay stratum. Across
   strata the cancellation fails and Landrum & Riniker's cross-assay noise
   returns. This is a hard eligibility constraint on panel construction.
5. **The DD noise estimate assumes independence** across the four cells. Shared
   ligands and shared references induce correlation; the effective sample size is
   below the nominal `C(k,2)`, and the table in §III.3 is therefore an upper
   bound. A cluster bootstrap over ligands and references is required.
6. **Low-rank `delta` is an assumption**, imported from AMMI practice in a
   different domain. If the true interaction is high-rank or sparse-and-spiky
   (a few activity cliffs rather than smooth loadings), the leading-`r` projection
   will miss it. A rank-sweep diagnostic on the source panel should precede
   DRP-A, and it is label-blind on source folds.

---

## Part VII — Why this is worth one panel

The project's last twenty stages share a shape: predict an interaction from
protein features, evaluate through a metric that cannot see it, and stop at a
null. DRP changes all three terms at once. It measures interaction instead of
predicting it; it is algebraically blind to the main effects that dominated every
previous statistic; and it is closed-form, so it cannot fail the way E0 through
E0R1 failed.

It also fits the frozen theory more tightly than anything currently in `model/`:
`d_z = 4–5` against a derived budget of 5–8, `d_adapt = 3–4` against `k = 5`,
sufficiency statable under a low-rank assumption, sign fixed by the verified
anchor ladder, and abstention computable from a condition number.

If it produces a null at `k = 20` with SNR ≈ 7 and clean controls, that null is
strong evidence — far stronger than any null in the ledger — that
`Var(delta) ≈ 0` in accessible data, and the project can close Claim B honestly
and write it up.

---

## Sources

- Frozen theory, read directly: `theory/FINAL_FROZEN_THEORY/chapters/{01_FOUNDATIONS, 02_TARGET_AND_OPERATOR, 04_META_LEARNING_FORMULATION, 06_CALIBRATION_AND_GENERALIZATION, 07_SCOPE_LIMITATIONS}.md`
- Relative representations enable zero-shot latent space communication (ICLR 2023 oral) — https://arxiv.org/abs/2209.15430
- Improving Relative Representations with Learned Anchors and Whitened Inner Products — https://arxiv.org/pdf/2605.30596
- Few-shot Conformal Prediction with Auxiliary Tasks (ICML 2021) — https://arxiv.org/abs/2102.08898
- Interaction Concordance Index — https://arxiv.org/abs/2510.14419
- Sample efficient inductive matrix completion with noise and inexact side information — https://arxiv.org/abs/2605.17189
- Speedup Matrix Completion with Side Information (NeurIPS 2013) — https://papers.nips.cc/paper/2013/file/e58cc5ca94270acaceed13bc82dfedf7-Paper.pdf
- Fast and Sample Efficient Inductive Matrix Completion via Multi-Phase Procrustes Flow (ICML 2018) — http://proceedings.mlr.press/v80/zhang18b/zhang18b.pdf
- AMMI decomposition of interaction — https://pmc.ncbi.nlm.nih.gov/articles/PMC6483959/
- Jarquín et al., reaction norm model — https://pubmed.ncbi.nlm.nih.gov/24337101/
- Kramer et al., Experimental Uncertainty of Heterogeneous Public Ki Data — https://pubs.acs.org/doi/10.1021/jm300131x
- Landrum & Riniker, Combining IC50 or Ki Values from Different Sources — https://pubs.acs.org/doi/10.1021/acs.jcim.4c00049

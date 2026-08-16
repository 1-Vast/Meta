# Biological knowledge as a load-bearing inductive bias for strict dual-cold DTA

**Date:** 2026-07-27
**Type:** literature review + diagnosis + architecture decision. **No experiment was run, no model trained,
no data read, no API called.** `sealed_test_consumed=false`; `confirmation_labels_read=true` (pre-existing).
**Scope:** model-side only. The LEXOR data-layer program (`task.md` Part 6) is owned by a separate task and
is not modified here.

---

## 0. Summary of the decision

| question | answer |
| --- | --- |
| Why did previous biological-knowledge attempts fail? | Three stacked causes, ranked: (C3) the identifying measurement variation was destroyed by database curation; (C1) the **row (protein) subspace is misspecified**, which provably annihilates the cold-start sample-complexity advantage no matter how `B` is parameterised; (C2) the protein coordinate is *pooled*, while the specificity signal is *sparse and positional*. |
| Keep `b(d)` + biological residual? | **Yes, unchanged.** `b(d)` is the only signal that reproduces provenance-family-disjoint (+0.700). |
| Replace the protein encoder? | **Yes — this is the single change.** Every failure since 2026-07-25 is a different parameterisation of `B` on top of the *same* pooled ESM-2 `u_t`. |
| Add privileged-information / distillation? | **No.** Rejected on the program's own information-headroom gate (`Δ_info +0.0154` vs threshold `0.045`) and on the closed physical-structure teacher. |
| Add causal / invariant objectives? | **No as a loss; yes as a split axis.** RECRO L0 proved the "environments" in this corpus are 91.6% fake; IRM on fake environments is vacuous. |
| Final proposal | **PARC** — Pocket-Anchored Row Coordinate — plus **γ-interpolation**, an estimator whose output *is* the load-bearingness of the biological coordinate, with a CI. Exactly 2 innovations. |
| Can it be confirmed today? | **No.** M0/M1 are train-only mathematical development. The predictive gate stays blocked behind a powered independent panel. This is stated up front, not discovered later. |

---

## 1. What the ledger actually says (facts, not narrative)

Five facts constrain any proposal. Each is cited to `history.md`.

**F1 — `b(d)` reproduces; `g(t,d)` does not.**
RECRO L0: raw within-target potency ordering reproduces across provenance-disjoint document families at
**+0.700**. The target-specific residual falls from `+0.3345 [0.219,0.446]` (all documents) to
`+0.0901 [-0.0557,+0.2358]` (document-family-disjoint), a 73% drop with `LCB < 0`; the matched
wrong-target contrast is `residual_correct − residual_wrong = −0.0744`. Verdict
`RECRO_SIGNAL_EXPLAINED_BY_PROVENANCE`.

**F2 — the interaction is real, low-rank, and *global*; it is not protein-indexed.**
On the dense Metz panel the interaction residual SD is `0.66292 pK` (48.9% of variance, distributed —
top-1% ligands carry 6.01% of the energy), reproduces on held-out cells (`PA4 LCB95 +0.2431`), and is
ESM-aligned (`PA5 p = 0.000488`, random-protein control `p = 0.2077`).
Then HQ-GBMA Stage D asked whether the protein predicts *which* low-dimensional subspace a target's own
interaction direction lies in, and got the decisive ordering:

```
global shared subspace  0.340 [0.163, 0.548]
protein(random)         0.314
protein(shuffle)        0.237
pooled-only             0.201
protein(TRUE ESM)       0.108 [-0.143, 0.324]      true − global = −0.232  (LCB −0.083)
```

A **shared** rank-6 basis contains 34% of held targets' noise-corrected interaction energy; the
ESM-conditioned map is *significantly worse than a protein-free shared basis*. PD-M independently found
3/8 projected directions both feature-explainable and component-transferable, with measured spectral
effective rank **22–26**, and with the ligand-side null easier to clear than the target-side null in every
direction.

**F3 — the protein coordinate is the variable that has never been changed.**
Pooled/frozen ESM-2 fails or is beaten by its own destruction controls in *five structurally independent*
mechanisms: `BM1_RR_FAIL_STOP` (random protein +0.0035), `CHEMBL_CFRI_GATE_Z_FAIL_STOP` (random exceeds
CFRI by 0.0083), `PANEL_GATE_PC_FAIL_STOP` (shuffled protein scores *higher*),
`PROTEIN_CONDITIONED_PRIOR_NOT_LOAD_BEARING` (both destroyed-protein controls at or above true),
`HQGBMA_STAGE_D_FAIL_STOP` (above). SCGD adds a sixth (protein-free SCGD beat protein SCGD).
The one time the coordinate itself was changed — C2 CAPIT/ASPIRE-P0, the KLIFS 85-residue **aligned pocket
composition** — the ordering inverted, over 302 strict homology components:

```
aligned pocket        0.4539
group centroid        0.4253      pocket − group    +0.0286 [+0.0116, +0.0459]   (threshold 0.030)
pooled ESM            0.3983      pocket − ESM      +0.0556 [+0.0403, +0.0714]
random same-group     0.3670      pocket − random   +0.0869 [+0.0665, +0.1068]
within-group shuffle  0.3562      pocket − shuffle  +0.0977 [+0.0761, +0.1188]
```

That is the **largest protein-specific margin anywhere in the program**, and it passed every destruction
control. Its two disqualifying caveats are recorded and binding: it was an intentionally *favourable
shared-ligand oracle* (ligand-warm, therefore not dual-cold), and it missed its own materiality threshold
by 0.0014. It is a coordinate diagnostic, not a result.

**F4 — few-shot support adaptation is the weaker half, not the rescue.**
`GATE_R_FAIL_FEW_SHOT_CLAIM_RETRACTED` (the +0.093 gain was an evaluation-protocol artefact),
`PANEL_GATE_PC_FAIL_STOP` (`−0.0007` vs B0; adapts the mean, not the order), `SCGD_FAIL_STOP`,
`QACO_FAIL_STOP`, `SI0_COVARIANCE_ALIGNMENT_FAIL_STOP` (+0.0271 < +0.030), and decisively
`RB_DR_QMAPD_ORACLE_INFORMATION_FAIL_STOP`: with an oracle teacher given up to 16 same-target labels,
`Δ_info = +0.0154 [-0.0155, +0.0464]` against threshold `max(0.03, 0.0452)`. Meanwhile the best-ranking arm
anywhere in the program is `I0`, the **zero-shot** unconstrained interaction (`+0.0376 [-0.0033, +0.0818]`
over B0 on Metz), which has no posterior, no cross-fitting and no support labels.

**F5 — the binding constraint is measurement design.**
`NO_OPEN_POWERED_INDEPENDENT_PANEL`. An admissible panel needs ~100 independent homology components **and**
~40 scaffold-diverse query ligands per target after firewalls. Metz has the shape but is spent and its
labels are rounded to 0.1 pK (81.7% within-target ties); Davis `MDE80 0.1596`; Reinecke `MDE80 +0.0668`
(median 5 query ligands); SPD median 14 compounds/assay; BindingDB-native 38 targets.

---

## 2. Diagnosis: why biological knowledge has not been load-bearing

### C1 — Inexact side information destroys the cold-start advantage. This is a theorem, not a hunch.

Strict dual-cold DTA *is* inductive matrix completion (IMC): predict `y(t,d)` for an unseen row and an
unseen column from row features `u_t` and column features `v_d`, through `y ≈ b(d) + u_t^T B v_d`. The IMC
literature settles what side information buys and when it stops buying it:

* With **exact / perfectly predictive** side information, sample complexity collapses from linear in the
  matrix dimension to `O(log n)`, and cold-start on both axes becomes possible at all
  (Xu et al. 2013; Jain & Dhillon 2013).
* With **noisy or partial** side information ("dirty IMC", Chiang, Hsieh & Dhillon, NeurIPS 2015), the
  degrees of freedom remain **linear in n**, the advantage is lost, and the target matrix is *not*
  recoverable — you need per-row observations, which dual-cold by construction forbids.
* Yang & Ma (arXiv:2605.17189, 2026) close part of this gap for the noisy case and show the estimation
  error degrades **with the level of subspace misspecification**, and propose a *penalized interpolation
  between IMC and ordinary matrix completion* trading sample efficiency against robustness to imperfect
  features.

This is the unification the ledger has been missing. CFRI (neural joint head), the Bayesian precision gate
(BM0/PC), the hierarchical covariance factor (HIER), the Grassmann subspace map (HQ-GBMA), the deep kernel
(SCGD), the antisymmetric operator (QACO) and the nuclear-norm shrinkers (ORRC/OSA) are **six different
parameterisations of `B` sitting on top of one and the same `u_t`**. If `u_t` is a misspecified row
subspace, dirty-IMC says none of them can succeed, and Stage D's `true − global = −0.232` says exactly
that: conditioning on a misspecified `u_t` is *worse than not conditioning at all*, because the map spends
its capacity fitting training-target idiosyncrasy.

**The failure is in the row feature map, not in `B`'s parameterisation.** Six negatives, one cause.

### C2 — The signal is sparse and positional; the coordinate is pooled.

Mean-pooling a 1280-d residue field over 300–1000 residues is a low-pass filter over exactly the residues
that determine specificity. The nearest sister field converged on this independently: in enzyme–substrate
specificity, "substrate scope is often determined by a select few residues within the active site", so a
whole-enzyme encoder "is likely to minimize or stifle their signal"
([arXiv:2607.05084](https://arxiv.org/html/2607.05084v1), July 2026). The same field also supplies the
sharpest warning about what "generalisation" means: ESP
([Nat. Commun. 14, 2787, 2023](https://www.nature.com/articles/s41467-023-38347-2)) reaches 88% accuracy on
enzymes at <40% identity to training — i.e. it *does* generalise on the protein axis — yet on the
**substrate** axis at the same identity band it collapses to `MCC ≈ 0.01`, statistically
indistinguishable from random. Dual-cold demands both axes at once; the published field has demonstrated
one.

FORT's own C2 result is this same finding in the DTA setting, and it is the only place in the ledger where
a protein coordinate cleared shuffle, random and pooled-ESM controls simultaneously.

### C3 — The identifying variation was removed by curation, and it gates C1/C2.

RECRO L0 is the program's most important measurement: 79.75% of co-measured cross-document cells are exact
duplicates; 984 ChEMBL documents collapse into 463 provenance families; 91.6% of nominal cross-document
comparisons are *intra*-family. Papyrus F0-P: one aggregated row per `(connectivity, target)` — zero
document-replicated cells. So the independent measurement environments that would identify `g(t,d)`
causally are largely absent from aggregated public data.

This is why C3 must be ranked first: **a perfect protein coordinate still cannot be confirmed** on any
currently registered open substrate at the 0.03 floor. Any honest plan therefore separates *estimating*
the coordinate's contribution (possible now, train-only) from *claiming* predictive improvement (blocked).

### Independent corroboration from 2026 literature

The program's private negatives are now public findings, which materially raises the value of the negative
result and lowers the risk that FORT is measuring a local artefact:

* **HonestAffinity** ([arXiv:2606.03422](https://arxiv.org/abs/2606.03422), June 2026) isolates exactly two
  priors — frozen ESM-2 650M embeddings and a learned pocket-position marker — under a leak-aware protocol
  on LP-PDBBind. Result: a **split-conditioned reversal**. Both priors improve on canonical/familiar splits
  and *reduce* Pearson R on every strict no-leak tier; the `Pocket-NoESM` variant wins on all strict tiers.
  Their recommendation — report *paired canonical and leak-proof ablations* — becomes a mandatory reporting
  rule below.
* **Target mirroring**
  ([bioRxiv 2026.06.29.735309](https://www.biorxiv.org/content/10.64898/2026.06.29.735309v1)): homologous
  proteins with *low* sequence identity still show highly correlated binding profiles; >6,000 such ChEMBL-36
  assay pairs; leakage persists at identity thresholds as low as **0.2**. Sequence-identity splitting is
  therefore insufficient, and a ligand-only baseline reaches `r = 0.66` on FEP+ — the exact signature FORT
  measures as "B0 is the only reproducible signal".
* **Graber et al.** ([Nat. Mach. Intell. 2025](https://www.nature.com/articles/s42256-025-01124-5)): after
  leak-free re-splitting, benchmark performance drops substantially, and some models "perform comparably
  well on CASF datasets after omitting all protein or ligand information from their input data."

### The five research questions, answered

1. **How should incomplete scientific knowledge enter a neural model?** Von Rueden et al.
   ([IEEE TKDE 2021 / arXiv:1903.12394](https://arxiv.org/abs/1903.12394)) enumerate four entry points:
   training data, **hypothesis set**, learning algorithm (regulariser), and final hypothesis (validation).
   FORT has tried the regulariser (mechanistic MPF field: `MECHANISM_PRIOR_NOT_IDENTIFIED`), the data
   (GO/structural pretraining: `PRETRAIN_SIGNAL_REAL_DOWNSTREAM_FAIL_STOP`), and extra branches (atom–residue
   field, six-channel interaction field, contact graphs — all removable or harmful). It has **never changed
   the hypothesis set through the coordinate system**. Dirty-IMC says that is the only entry point that can
   change the cold-start rate.
2. **How to prevent a model from ignoring the knowledge?** Structurally — make it the sole path to the
   quantity of interest (FORT already does this correctly: `A0` provably cannot reorder; `k=0` is exact).
   But the opposite failure is now the dominant one: *forcing* a knowledge path that is harmful
   (BridgeFIRE S1-R physics gate `−0.256 [LCB −0.395]`; HIER protein prior anti-load-bearing; HonestAffinity's
   strict-tier reversal). Rule adopted: **never hard-wire; gate every knowledge path by a learned scalar that
   is permitted to reach zero, and report that scalar as a primary result.**
3. **How to verify the knowledge is actually used?** FORT's battery (shuffle / random / matched-wrong /
   pooled-composition / cross-target / label-permutation / exact `k=0`) is already at or above the published
   standard. One upgrade is forced by Stage D: **"beats random" is the wrong null.** The correct null for a
   protein-conditioned module is a *capacity-matched, protein-free shared model*. Two further controls are
   new here (§4.4): pocket-set shuffle and structure-token shuffle, which separate "the pocket" from "the
   protein" and "local structure" from "residue identity".
4. **How can scientific priors improve OOD?** Only by reducing the *effective dimension of a genuinely
   predictive row subspace*. That is literally the IMC condition. Priors that add capacity without reducing
   subspace misspecification cannot improve OOD and empirically degrade it.
5. **How should imperfect biological knowledge be represented?** As a **shrinkage target with an estimated
   weight**, not as a hard structure — Yang & Ma's penalized interpolation, i.e. a hierarchical model whose
   prior mean is the shared global object and whose protein-conditioned deviation is an explicitly
   penalised, ablatable increment.

---

## 3. Candidate comparison (all seven required alternatives)

Ranked by admissibility under the ledger. "Compat." = compatibility with the frozen contract.

| # | candidate | why it may work | why it fails / is closed here | data required | complexity | compat. | novelty | verdict |
|---|---|---|---|---|---|---|---|---|
| 1 | **Structure-compressed protein coordinate + shared low-rank residual (PARC)** | Directly attacks subspace misspecification, the only quantity dirty-IMC says controls the cold-start rate. In-program evidence: aligned pocket beats pooled ESM `+0.0556 [+0.0403,+0.0714]` with all destruction controls passing. Cross-domain: active-site conditioning is the convergent conclusion of the enzyme-specificity field. | Risk: AF apo pockets miss cryptic sites; HonestAffinity shows pocket markers can *hurt* on strict tiers. Both are declared as arms, not assumed away. | Local KLIFS (353/358 genes) for the kinase instantiation; AFDB + a 3Di tokeniser for the multi-family generalisation (dependency **not yet audited**) | Low (~10⁴ params on top of frozen features) | High — reuses `M_X^W` projection, φ(d), and EC-Helmert-EB verbatim | Medium-high: the *estimand* (load-bearingness of a biological coordinate) is the novel object | **ADOPT** |
| 2 | Bayesian biological prior `p(θ | knowledge)` | Principled; already implemented | **Structurally dead.** ORRC-EB v2 §1 *derives* that a positive coordinate-wise target-conditioned precision gives `sign(m_j) = sign(b_j)` for every positive `v` — it can shrink, never rotate or flip. Confirmed empirically 3× (`PC`, `HIER`, `BM1-RR`). | — | — | — | none | **REJECT (proved non-identifiable)** |
| 3 | Privileged information (structure at train, sequence at test) | LUPI has real precedent for protein affinity ([PMC6238365](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6238365/); [arXiv:2601.03704](https://arxiv.org/abs/2601.03704)) | Two independent kills. (i) FORT already ran the *information-headroom precondition*: `Δ_info +0.0154 [-0.0155,+0.0464]` vs threshold `0.0452`, and its own §7 rule says "do NOT implement privileged completion distillation". (ii) The structure teacher is closed: `BRIDGEFIRE_S1_FAIL` shows that even at the **native crystal-pose upper bound** the joint model (0.257) loses to the chemistry-marginal control (0.394) and to ligand heavy-atom count alone (0.429). Distilling a teacher that fails its own load-bearing test cannot help. Note also that AFDB coverage makes structure *universally available*, so it is not privileged in the first place. | — | — | — | low | **REJECT** |
| 4 | Knowledge distillation (structure teacher → sequence student) | Same literature as #3 | Same as #3; additionally `MECHANISM_PASS_UTILITY_FAIL_STOP` (residue-teacher pocket distillation: KL 2.54→0.074, RMSE +0.065/+0.048/+0.041) is the exact experiment already run and failed. Even the 2026 LUPI paper reports the student does not reach the structure teacher and attributes it to dataset scale — FORT's scale is smaller. | — | — | — | low | **REJECT** |
| 5 | Causal / invariant learning (IRM-family) to remove assay & provenance shortcuts | Correct diagnosis of the confound | IRM needs ≥2 *genuinely independent* environments per unit. RECRO L0: 91.6% of nominal cross-document comparisons are intra-family — the environments are fake, so the invariance penalty is vacuous. IRM is also known to lose to plain ERM in weakly-diversified high-dimensional regimes. | would require the LEXOR/prospective layer | — | — | — | **REJECT as loss; ADOPT as split axis** (provenance-family firewall + new binding-profile-correlation firewall) |
| 6 | Energy-based / physics-inspired | Mechanistically appealing | `PHYSICAL-STRUCTURE PROGRAM CLOSED` — falsified across four framings (deployment docking ρ≈0.005; holo vs AF paired Δ −0.019; native-pose standalone; native-pose increment `−0.256 [LCB −0.395]`). Reopening requires a *new* registered mechanism justification, which none of this supplies. | — | — | — | — | **REJECT (closed route)** |
| 7 | Mechanistic regularisation | Interpretable | MPF 7-basis signed interaction field: not identifiable from scalar affinity supervision (`MPF_VIABLE_BUT_MECHANISTIC_FIELD_NOT_YET_LOAD_BEARING`); STRATA-MX R2 structural prior `+0.0002 [-0.0001,0.0004]`, identical to random gate and shuffled pocket. | — | — | — | — | **REJECT** |
| 8 | Graph / equivariant architectures | SOTA elsewhere | Subsumed by #6's native-pose upper bound, and data-blocked (no local complexes; AFDB apo only). Also, `A0`/`I0` evidence says the bottleneck is not expressiveness — the *unconstrained* head is already the best arm. | — | — | — | — | **REJECT** |

---

## 4. The proposal: PARC + γ-interpolation

### 4.1 Mathematical formulation

Keep the frozen contract and the exact ORRC observed-edge projection. Write the model on the projected
residual, with main effects profiled out by the weighted projector `M_X^W`:

```
y(t,d) = b(d) + g(t,d) + eps
g(t,d) = phi(d)^T w_t
w_t    = w_bar + gamma * Delta(u_t)                      (1)
```

* `phi(d) in R^r` — the **shared** ligand-side interaction basis, `r` chosen by nested folds (PD-M measured
  effective rank 22–26, so `r` is *not* fixed at 8; the old ceiling was inherited from the provably
  non-identifiable adapter dimension and is dropped).
* `w_bar in R^r` — the **shared global** coefficient. This is the arm that *won* Stage D (0.340 vs 0.108).
  It is the model's prior mean and its protein-free fallback.
* `Delta(u_t) = W_2 sigma(W_1 u_t)`, with `W_1 in R^{h x d_u}`, `W_2 in R^{r x h}`, `h ≈ 16`. Deliberately
  **lower capacity than the shared basis** — Stage D's failure mode was capacity spent on training-target
  idiosyncrasy.
* `gamma in [0,1]` — a single learned scalar, initialised at 0, penalised toward 0.

`gamma` is the penalized interpolation of Yang & Ma (2026): `gamma = 0` is ordinary shared low-rank
completion (robust, no side information); `gamma = 1` is full inductive matrix completion (sample-efficient
*iff* `u_t` is exact). Its estimate `gamma_hat` with a component-bootstrap CI **is the primary reported
quantity**: it is a direct, interpretable, preregisterable measurement of how load-bearing the biological
coordinate is, and it is informative whether it is positive or zero.

Estimation uses the existing **EC-Helmert-EB** machinery, repurposed rather than kept as a separate
innovation: the per-target empirical coefficient `w_hat_t` comes with a sandwich covariance `V_t`, and
`w_hat_t ~ N(w_bar + gamma*Delta(u_t), Sigma_0 + V_t)`. Subtracting `V_t` is what makes `gamma_hat`
unbiased; treating cross-fitted coefficients as noise-free is exactly the error blueprint v3 identified,
and it would bias `gamma_hat` upward — i.e. it would *manufacture* the positive result we are testing for.

### 4.2 The innovation itself — PARC, the Pocket-Anchored Row Coordinate

`u_t` is built label-free, with no ligand, no pose, no complex, no docking, no taxonomy label:

1. **Pocket residue set `S_t`** (`|S_t| ≈ 30–60`), obtained label-free. Two registered instantiations:
   * *Kinase instantiation (available now, zero new dependency):* the KLIFS 85-residue aligned pocket,
     already local and already used by `research/kirhub_pocket_oracle.py`.
   * *Family-agnostic instantiation (required for any multi-family substrate):* geometric cavity detection
     on the AlphaFold monomer. **Dependency not yet audited** — the local `dataset/structure/alphafold`
     cache was removed during the 2026-07-25 cleanup, and no Foldseek/3Di tooling is installed in the `drug`
     env. This is an explicit E-gate item, not an assumption.
2. **Per-residue encoding**, position-preserving: (i) the frozen **ESM-2 residue** vector — the residue axis
   the program has always had and always pooled away; (ii) a discrete **local-structure token** in the
   Foldseek 3Di 20-state alphabet, i.e. the SaProt structure-aware vocabulary
   ([bioRxiv 2023.10.01.560349](https://www.biorxiv.org/content/10.1101/2023.10.01.560349v5.full)), which is
   precisely "structure compressed into a sequence-shaped protein representation" and is *not* a physical
   interaction channel.
3. **Reduction to `u_t in R^{d_u}`**, `d_u ≈ 32–64`: a low-rank projection of the pocket residue field plus
   the 3Di composition histogram over `S_t`. The *sparsity* is the knowledge; the pooling inside the pocket
   is kept deliberately dumb so that any gain is attributable to the restriction, not to a new pooling
   mechanism.

Deployment: AFDB covers essentially every UniProt accession, so this is sequence-only-at-deployment in the
operational sense. A sequence-predicted-3Di fallback arm is registered and **reported separately, never
summed** — measured, not distilled.

### 4.3 Information flow, and where biological knowledge enters

```
sequence ──► frozen ESM-2 residues ─┐
                                    ├─► restrict to S_t ─► u_t (d_u ≈ 32–64) ─► Delta(·) ─┐
AF structure ──► 3Di tokens ────────┘                                                     │
                                                                          w_bar ──► + ────┴─► w_t
ligand ──► Morgan+descriptors ──► b(d)  and  phi(d) ──────────────────────────────────────────► y
```

Knowledge enters at **exactly one place: the hypothesis set**, as a restriction of the row coordinate to a
geometrically-defined pocket in a structure-aware alphabet. It does **not** enter the loss, the prior over
`theta`, or as an additional branch — the three entry points that have already failed here.

### 4.4 Ablation and destruction battery

Primary contrast is **not** "beats random". It is:

| # | control | what it isolates | requirement |
| --- | --- | --- | --- |
| A1 | `gamma = 0` shared-global model | the correct null (the Stage-D winner) | `PARC − shared` LCB95 > 0 |
| A2 | **pocket-set shuffle** — random size-matched residue set from the *same* protein | the pocket restriction, holding length, composition marginals and ESM scale fixed | LCB95 > 0 — **new to this program** |
| A3 | **3Di shuffle within `S_t`** — permute structure tokens, keep residue identity | the structural half vs the sequence half | reported; may legitimately be null |
| A4 | matched wrong target (size- and family-matched) | target identity | LCB95 > 0 |
| A5 | random features at matched `d_u` | capacity | LCB95 > 0 |
| A6 | cross-target coefficient swap (`w_t` from another target) | that `w_t` is target-specific at all | LCB95 > 0 |
| A7 | capacity-matched protein-free control (same params spent on `phi`/`w_bar`) | that the gain is not capacity | LCB95 > 0 |
| A8 | **synthetic `gamma`-recovery** — simulate with known `gamma*`, check `gamma_hat → gamma*` | that a null `gamma_hat` means "no signal", not "broken estimator" | recovery within CI — **new; protects the negative result** |
| A9 | knowledge-conformity readout (von Rueden output stage) | that predicted `w_t` similarity tracks *pocket* similarity more than whole-sequence similarity | reported; if it tracks whole-sequence similarity the restriction is inert |
| A10 | **paired canonical / strict-firewall reporting** | HonestAffinity's split-conditioned reversal | both reported side by side, always |

A8 is the control that makes a null publishable. A2 is the control that makes a positive result meaningful.
Neither exists in the current battery.

### 4.5 Training objective and inference

Train-only, cross-fitted on held homology components, on the exact ORRC projected residual:

```
min over (phi, w_bar, W_1, W_2, gamma)
    sum over observed edges  w_e * ( M_X^W y )_e  -  phi(d_e)^T w_{t_e} )^2
  + lam_B * ||shared low-rank block||_*        (frozen convex nuclear norm, eps = 1e-8 lam_B for uniqueness)
  + lam_D * gamma * ||Delta||_F^2              (the interpolation penalty; lam_D by nested folds)
```

Within-target pairwise ranking is retained as the auxiliary objective — it is the program's one validated
reusable technique (BridgeFIRE: 0.045 → 0.257, scale-invariant, and the natural loss when a target's
absolute baseline is unknowable). No support labels, no episodes, no meta-learning on the critical path.

Inference is a single forward pass: `y_hat = b(d) + phi(d)^T (w_bar + gamma_hat * Delta(u_t))`.
Compatibility with the Bayesian posterior contract is preserved but *not exercised*: `phi(d)` is exactly the
design matrix a later support posterior would use, and with zero support labels the `k=0` contract is
satisfied vacuously.

### 4.6 Why this should improve dual-cold, and why it avoids each prior failure

| prior failure | why PARC is not that |
| --- | --- |
| atom–residue field, six-channel interaction field, contact graphs | no interaction channel at all; `u_t` is protein-only |
| docking / BFEO / BridgeFIRE / Gate-P | no pose, no complex, no docking score, no ligand-conditioned geometry |
| protein-conditioned *precision* (BM0, PC) | ORRC v2 §1 proves precision can only shrink; PARC moves the coefficient **mean**, so it can rotate and flip sign |
| protein-conditioned *covariance* (HIER) and *subspace* (HQ-GBMA Stage D) | those conditioned the subspace and lost to a shared basis; PARC **keeps the subspace shared** — it adopts Stage D's winner as its null and its prior mean |
| pooled ESM everywhere | sparse, positional, structure-aware coordinate |
| few-shot posterior (Gate R, PC, SCGD, QACO, SI0, O1) | off the critical path |
| GO / family / taxonomy embeddings | forbidden and unused; no taxonomy label enters `u_t` |
| unconstrained head absorbing the base (the reason CFRI was constrained) | main effects removed exactly by `M_X^W`, not by a soft orthogonality penalty |

---

## 5. Comparison with the existing blueprint (ORRC-EB v3)

| component | v3 status | under this proposal |
| --- | --- | --- |
| Exact observed-edge projection `M_X^W` | retained, KKT-audited (`1.17e-13`) | **retained verbatim** |
| Shared low-rank `B` under fixed nuclear norm | train-only reference | **retained** — becomes `phi`, `w_bar` |
| OSA-ORRC monotone singular-value shrinker | `OSA_ORRC_ARCHITECTURE_FAIL_REVIEW` (`−0.0037 [-0.0124,+0.0045]`, boundary selection) | **retired.** PD-M2A showed the remaining headroom is not in `B`'s estimator |
| EC-Helmert-EB error-corrected prior | preregistered innovation | **repurposed, not retired** — it becomes the unbiased estimator for `gamma`, which is its highest-value use |
| Rank ceiling `r <= 8` | recorded obsolete (measured 22–26) | dropped; `r` selected by nested folds |
| `panel_davis` as confirmation source | recorded obsolete (underpowered) | unchanged; stays sealed, `consumed=false` |
| Protein coordinate | **never a variable** | **the single innovation** |
| Innovation budget | 2 (OSA-ORRC, EC-Helmert-EB) | 2 (PARC coordinate; `gamma`-interpolation estimator + identifiability battery) — no increase |

The relationship is succession, not competition: ORRC-EB spent its budget improving the estimator of `B`
and measured that path flat; PARC says the residual headroom is in `u_t` and reuses ORRC's projection and
EB machinery unchanged.

---

## 6. Staged plan, and the honest limit

**M0 — coordinate-misspecification audit (train-only, non-gating, ~1 GPU-hour).**
Re-run HQ-GBMA Stage D with **one changed input**. Same estimator, same folds, same error-corrected
containment, same `r=6`, same seed. Arms: shared-global / pooled-ESM / **PARC (KLIFS-aligned pocket)** /
pocket-set shuffle / 3Di shuffle / random / matched-wrong. This is the cheapest discriminating test in the
entire design and it answers the program's central question directly: *is the protein uninformative, or was
the coordinate wrong?* If PARC does not beat pooled ESM **and** does not close the gap to the shared-global
arm, the route dies for one GPU-hour and the program gains a clean, quantitative, publishable answer.

**M1 — `gamma` estimation (train-only, non-gating).** Fit (1) with EC-Helmert-EB, report `gamma_hat` with a
component-bootstrap CI, the full A1–A10 battery, and the A8 synthetic recovery check. Deliverable: a
*measurement* of biological-coordinate load-bearingness, not a performance claim.

**M2 — predictive gate. BLOCKED.** Requires a registered, power-verified independent panel
(~100 independent components × ~40 scaffold-diverse query ligands after firewalls) — i.e. `LEXOR L4/L5`
or a prospective factorial panel. No result from M0 or M1 may be described as predictive evidence or used
to authorise M2, F1–F4, or any Mamba comparison.

**New firewall axis, mandatory from now on.** Sequence identity is registered as **insufficient** (target
mirroring: leakage persists at identity 0.2). Every split must add a **binding-profile-correlation
firewall**: two targets fall in the same holdout block if their measured cross-ligand binding profiles
correlate above a frozen threshold, regardless of sequence identity.

**Declared expected outcome.** The most likely M0 result is that PARC beats pooled ESM (C2 predicts a
`+0.05`-scale margin) but still does not beat the shared-global basis. That is a *useful* outcome: it would
localise the residual failure to the coefficient map rather than the coordinate, and — combined with
HonestAffinity's independent split-conditioned reversal and the dirty-IMC theorem — would let the program
state the strongest form of its negative result: **strict dual-cold protein-conditioned SAR reordering is
not identifiable from any currently available protein coordinate, and this is a subspace-misspecification
result with a theorem behind it, not a tuning failure.**

---

## Sources

- [Informed Machine Learning – A Taxonomy and Survey of Integrating Knowledge into Learning Systems (von Rueden et al., arXiv:1903.12394 / IEEE TKDE 2021)](https://arxiv.org/abs/1903.12394)
- [Sample-efficient inductive matrix completion with noise and inexact side information (Yang & Ma, arXiv:2605.17189)](https://arxiv.org/abs/2605.17189)
- [Matrix Completion with Noisy Side Information (Chiang, Hsieh & Dhillon, NeurIPS 2015)](https://www.cs.utexas.edu/~inderjit/public_papers/dirtyIMC_nips15.pdf)
- [Fine-grained Generalisation Analysis of Inductive Matrix Completion (NeurIPS 2021)](https://proceedings.neurips.cc/paper/2021/file/d6428eecbe0f7dff83fc607c5044b2b9-Paper.pdf)
- [SaProt: Protein Language Modeling with Structure-aware Vocabulary (bioRxiv / ICLR)](https://www.biorxiv.org/content/10.1101/2023.10.01.560349v5.full)
- [HonestAffinity: Leak-Aware Evaluation of Protein and Pocket Priors for Binding Affinity Prediction (arXiv:2606.03422)](https://arxiv.org/abs/2606.03422)
- [Identifying and Addressing Systematic Data Leakage in Protein-Ligand Affinity Benchmarks — "target mirroring" (bioRxiv 2026.06.29.735309)](https://www.biorxiv.org/content/10.64898/2026.06.29.735309v1)
- [Resolving data bias improves generalization in binding affinity prediction (Graber et al., Nat. Mach. Intell. 2025)](https://www.nature.com/articles/s42256-025-01124-5)
- [Learning to Generalize: Deep Models and Robust Benchmarks for Drug–Target Affinity Prediction (ChemRxiv 2025)](https://chemrxiv.org/doi/full/10.26434/chemrxiv-2025-gmrdb)
- [Rethinking Benchmarks and Models for Enzyme Specificity Prediction (arXiv:2607.05084)](https://arxiv.org/html/2607.05084v1)
- [A general model to predict small molecule substrates of enzymes — ESP (Nat. Commun. 2023)](https://www.nature.com/articles/s41467-023-38347-2)
- [Learning protein binding affinity using privileged information (Abbasi et al., PMC6238365)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6238365/)
- [Investigating Knowledge Distillation Through Neural Networks for Protein Binding Affinity Prediction (arXiv:2601.03704)](https://arxiv.org/abs/2601.03704)
- [Learning with Privileged Knowledge Distillation for Improved Peptide–Protein Docking (ACS Omega 2025)](https://pubs.acs.org/doi/10.1021/acsomega.5c00967)
- [PocketDTA: multimodal architecture using 3D structural data of target binding pockets (Bioinformatics 2024)](https://academic.oup.com/bioinformatics/article/40/10/btae594/7811139)
- [Sequence-based drug-target binding site pre-training enables cryptic pocket detection (J. Cheminform. 2026)](https://link.springer.com/article/10.1186/s13321-026-01227-0)
- [Out-of-distribution Generalization for Total Variation based Invariant Risk Minimization (arXiv:2502.19665)](https://arxiv.org/html/2502.19665v2)
- [Generalization Beyond Benchmarks: Evaluating Learnable Protein-Ligand Scoring Functions on Unseen Targets (arXiv:2512.05386)](https://arxiv.org/html/2512.05386v1)

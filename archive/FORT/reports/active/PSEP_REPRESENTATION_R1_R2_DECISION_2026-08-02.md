# PSEP — invariant representation learning: decision report (R1, R2)

Date 2026-08-02 · Runners `research/psep_invariance.py`, `research/psep_representation.py`
Artifacts `reports/active/psep_{invariance,representation}_2026-08-02.json` (+ parquets)
Role read: **`discover` only**. `validate` and `confirm` never opened.

**Verdict on the registered gate: `NO_ARM_BEATS_STANDARD_MULTITASK_DTA`** — the
required margin (lower95 > 0.005 over standard multitask) is not met.
**But one mechanism is identified, derived, and confirmed**, and the entire
invariance family is falsified. Both are stated precisely below.

---

## 1. The formalisation, measured

    y = f(x,p) + g(c) + eps

Variance decomposition of the residual after a target-agnostic linear base, on
263 318 discover rows (one-way, groups with >= 5 rows):

| factor | groups | R² | group-mean SD |
|---|---:|---:|---:|
| **assay** | 9 158 | **0.581** | 1.092 |
| **document** | 7 838 | **0.544** | 1.050 |
| ligand identity | 1 911 | 0.364 | 0.669 |
| target unit | 823 | 0.248 | 0.747 |
| accession | 616 | 0.214 | 0.729 |
| homology component | 508 | 0.206 | 0.710 |
| **endpoint** | 3 | **0.005** | 0.240 |

Nested: **unit × document 0.642** vs **unit alone 0.261**.

`g(c)` is roughly **2× everything about which protein you are looking at**. That
is the quantitative reason every target-conditional mechanism in this programme
drowned. Endpoint contributes essentially nothing (R² 0.005), retroactively
validating the decision to pool pKi/pIC50/pKd with `unit = accession|endpoint`.

## 2. The ceiling — measured, not assumed

Intra-document replicate noise, from 83 236 aggregates whose replicates sit
inside one document (median SD; sigma corrected for the half-normal bias of an
N=2 estimate, `sigma = median / 0.6745`):

| endpoint | within-doc label SD | sigma | **CI ceiling** | headroom over base 0.533 |
|---|---:|---:|---:|---:|
| pIC50 | 0.589 | 0.279 | **0.843** | **+0.310** |
| pKi | 0.595 | 0.440 | **0.735** | **+0.202** |
| pKd | 0.642 | 1.772 | **none** | noise exceeds signal |

Two consequences: **pKd is unlearnable within documents** and should be excluded
from headlines; **pKi is ~2× noisier intra-assay than pIC50**, inverting the usual
assumption that Ki is the cleaner endpoint.

Every effect in this programme should be read against **+0.20 … +0.31**, not
against 1.0.

## 3. Results

### R1 — six arms, matched capacity, 76 meta_test components

| arm | gain vs base | vs ERM | context R² |
|---|---:|---:|---:|
| erm (MSE, selected on val **MSE**) | −0.0025 | — | 0.636 |
| fixed_effects | +0.0007 | +0.0032 [+0.0010,+0.0054] | 0.635 |
| centred | +0.0170 | +0.0195 [−0.0042,+0.0434] | 0.326 |
| irm | −0.0062 | −0.0037 | 0.447 |
| dro | +0.0084 | +0.0109 [−0.0056,+0.0274] | 0.503 |
| adversarial | +0.0095 | +0.0120 [−0.0003,+0.0239] | **0.825** |

**The adversarial arm produced the *highest* document dependence** (0.825 vs ERM
0.636) while being trained by gradient reversal specifically to suppress it. Its
CI gain therefore did not come from its claimed mechanism. Had only the CI column
been reported, this would have read as "context-adversarial learning helps".

### R2 — 2×3 factorial, 5-fold component cross-fitting, **373 components**

Every component tested exactly once; all arms share folds, encoder, capacity,
budget, and are **model-selected on within-document concordance** (so no arm is
penalised for being tuned against a different quantity — this alone moved
`mse_raw` from −0.0025 in R1 to +0.0084).

| arm | gain vs base | vs mse_raw |
|---|---:|---:|
| mse_centred | +0.0069 [−0.0041,+0.0170] | −0.0015 [−0.0119,+0.0086] |
| mse_fixed | +0.0082 [−0.0025,+0.0185] | −0.0002 [−0.0026,+0.0020] |
| mse_raw | +0.0084 [−0.0023,+0.0188] | — |
| **rank_raw** | **+0.0192 [+0.0087,+0.0295]** | **+0.0108 [+0.0009,+0.0207]** |
| **rank_centred** | **+0.0192 [+0.0087,+0.0295]** | **+0.0108 [+0.0009,+0.0207]** |

**Decomposition: objective +0.0108 · context −0.0015 · interaction +0.0015.**

## 4. The mechanism: the objective *is* the invariance

`rank_centred` is numerically identical to `rank_raw` — same CI to four decimals,
same epoch count, every fold. This is a mathematical identity, not a coincidence:
the ranking loss is computed on **within-document pairs**, and

    (y_i − m_d) − (y_j − m_d) = y_i − y_j

so the document mean is invisible to every gradient. **`g(c)` is structurally
unidentifiable in a within-document pairwise objective.**

That single fact explains the entire invariance null across R1 and R2. IRM,
GroupDRO, adversaries, fixed effects and centring are all solving a problem the
objective has already solved.

It also explains why MSE-based nuisance removal is *not* a substitute:

| approach | what it removes | result |
|---|---|---:|
| `mse_centred` | nuisance **and** real across-document signal — documents genuinely differ in potency | +0.0000 |
| `mse_fixed` | nuisance via free intercepts, but the target still carries it while fitting | +0.0019 |
| **`rank_raw`** | nuisance from the **comparison**, leaving `f` unconstrained in level | **+0.0108** |

Centring removes the nuisance from the *target* and destroys real signal with it.
Pairwise within-document ranking removes it from the *comparison* only. That
distinction is the mechanism, it is derivable in advance, and the three MSE arms
are its ablation.

## 5. Gate, honestly scored

| # | Requirement | Result |
|---|---|---|
| 1 | improves held-out provenance-separated targets | **yes** — +0.0192 [+0.0087,+0.0295] over the frozen base, 373 components |
| 2 | removes document dependence | **yes, structurally** (identity above) — not via a learned penalty |
| 3 | outperforms standard multitask DTA | **+0.0108 [+0.0009,+0.0207]** — excludes zero, **but the lower bound is below the registered 0.005 MDE** |
| 4 | survives document split | **yes** — all evaluation is on scaffold+document+assay-separated rows |
| 5 | clear ablation | **yes** — three MSE arms; and `rank_centred ≡ rank_raw` is a bit-identical no-op ablation |

**Requirement 3 is not met at the registered margin.** The correct statement is
that the ranking objective is the only intervention that significantly beats the
frozen base, and it beats standard multitask MSE by +0.0108 with an interval
excluding zero but not clearing the pre-registered threshold.

## 6. What this corrects

The "+0.0274 target-agnostic representation improvement" reported after the
operator gate was **substantially an objective effect, not a representation
effect**. The same feature space under MSE gives +0.0084 (n.s.); under the
ranking surrogate, +0.0192. Attributing it to nonlinear capacity was wrong; the
uncontrolled variable was the loss.

## 7. Position and remaining risk

Against measured headroom (+0.20 pKi / +0.31 pIC50), +0.0192 is **6–10 % of what
is obtainable**. This is a real, derived, ablated mechanism — and a small one.

- **Single seed.** Cross-fitting controls fold noise (fold-to-fold spread ~0.016,
  comparable to the effect); seed variance is not yet measured. Multi-seed is the
  first thing required before any claim.
- **The margin is thin.** +0.0108 [+0.0009, +0.0207] would not survive a stricter
  MDE. More components (`validate`/`confirm` are sealed and unused) or more seeds
  are the honest routes to resolving it.
- **pKd should be dropped** from any headline (§2).
- **`centred` RMSE in R1 is invalid** (restored absolute scale using evaluation
  document means); only its ranking metric was used.
- **`context_r2` is an in-sample probe** and should be re-measured cross-fitted if
  document dependence ever becomes load-bearing.

**STOP conditions in force.** No further invariance mechanism may be proposed on
this substrate — the objective already quotients `g(c)` exactly, so there is
nothing left for one to remove. Any future gain must come from the function class
or from data, not from nuisance handling.

---

## 8. R3 — multi-seed replication (post-hoc addendum, same day)

`research/psep_seeds.py`, 5 seeds x {mse_raw, rank_raw} x 5 folds, folds held
**fixed across seeds** (`component_fold` keys on the frozen programme seed, not
the training seed) so the per-component paired contrast is comparable across
seeds and only training randomness varies.

| seed | rank − mse |
|---|---:|
| 20260802 | +0.0108 [+0.0009,+0.0207] |
| 20260803 | +0.0035 [−0.0061,+0.0128] |
| 20260804 | +0.0169 [+0.0066,+0.0266] |
| 20260805 | +0.0033 [−0.0076,+0.0139] |
| 20260806 | +0.0140 [+0.0033,+0.0251] |

**Positive in 5/5 seeds** (sign stable; seed-contrast mean +0.0097, SD 0.0061).
Seed-averaged, component bootstrap, n=373: **rank − mse = +0.0097 [+0.0013,
+0.0180]** — tighter than the single-seed R2 estimate on both ends, but the lower
bound remains **below the pre-registered 0.005 MDE**.

**Verdict: `OBJECTIVE_EFFECT_REPLICATES_BUT_BELOW_MDE`.** The mechanism is real,
reproducible across training randomness, mathematically derived, and passes its
own ablation — and it is not large enough to clear the threshold set before this
result was seen. Per standing programme discipline, this is not rescued by
further seeds, a relaxed MDE, or a different metric. ~3–5 % of measured headroom
(+0.20 pKi / +0.31 pIC50).

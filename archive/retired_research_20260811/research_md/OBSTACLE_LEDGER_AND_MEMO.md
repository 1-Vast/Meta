# Meta-fewshot obstacle ledger and theory-to-model memo

Updated 2026-08-10, after Phase 0 (`FEWSHOT_EPISODE_DATA_NOT_IDENTIFIABLE`,
commit `3bc5ed5`).

## 1. Current-state diagnosis

| dimension | state |
|---|---|
| execution authority | `task.md`; episodic stage authorized, Phase 0 failed closed |
| frozen theory | `A(F,z)=K(B(z)F(z))` unchanged; simplex, ridge, band, mesh untouched |
| biological coordinate | frozen 288D T-BASIS, structurally validated, **not affinity-admitted** |
| meta-learning evidence | none; no model preregistered or trained |
| earliest unresolved obstacle | **O1, evaluation identifiability** |

Phase 0 measured, label-blind and with zero affinity reads: leakage exactly
zero on target, ligand, scaffold, document and protein-homology-40; source
supply ample (442 targets, 220 at `k=5`); evaluation depth insufficient
(68 targets, 16 at `k=5`, `MDE_d = 0.622` against a declared `0.600`).

## 2. Obstacle ledger

### O1 — Evaluation identifiability — **RESOLVED AS A CORPUS ARTIFACT**

```text
Obstacle          held-out panel too shallow to identify the few-shot effect
Causal hypotheses H1 the data genuinely lack independent target depth
                  H2 depth exists but was filtered out by a corpus built for a
                     different estimand
                  H3 depth requires a new source (PDSP non-kinase stratum)
Minimal discriminator
                  label-blind recount of the SAME governed projection under
                  few-shot rules only (single chain, Ki, >= k+3 ligands),
                  dropping the cycle-positive quotient requirement
Result            H2 confirmed; H1 rejected; H3 not needed yet
Next decision     preregister a few-shot-shaped corpus rebuild with full closure
```

Evidence:

| quantity | quotient-shaped corpus | few-shot rules only |
|---|---:|---:|
| single-chain Ki rows | 12,457 cells | **25,072** |
| distinct Ki targets | 510 proteins | **910** |
| targets usable at `k=1` | — | 584 |
| targets usable at `k=2` | — | 499 |
| targets usable at `k=3` | — | 459 |
| targets usable at `k=5` | 236 total (220 src + 16 eval) | **394** |
| `>=8` ligands and `>=2` documents | — | **218** |

`MDE_d` by held-out panel size: `16 -> 0.622`, `30 -> 0.454`, `50 -> 0.352`,
`100 -> 0.249`, `150 -> 0.203`.

**Root cause.** The CQ corpus required *cycle-positive* panels, because the
crossed-rectangle quotient estimand needs closed rectangles. The few-shot
estimand needs only per-target `(ligand, affinity)` depth. Inheriting the
quotient filter discarded roughly half the usable few-shot material and, worse,
preferentially discarded exactly the targets that make good independent
evaluation tasks: those appearing in one or few panels without completing
rectangles.

**This is not threshold-moving and not re-cutting dependent data.** It is
building the corpus for the estimand actually being tested, which is what
O1 prescribes as label-blind expansion of independent target panels.

### Honest limits of this result

1. The 394 targets have **not** passed protein-40 homology, scaffold and
   document conflict closure. Closure will reduce independent components, and
   the giant-component pathology (85.86% in one component) may persist.
2. Target depth is demonstrated; **independent component depth is not**. Only a
   full closure run can decide whether the expanded panel yields enough
   independent evaluation units.
3. The conservative figure is `218` targets with `>=8` ligands and `>=2`
   documents. If closure preserves even half of those as independent evaluation
   tasks, `MDE_d` lands near `0.20-0.25`.

### O2 — Biological information — **UNTESTED**

`phi = ` frozen 288D T-BASIS. `CQ_TBASIS_LINEAR_AFFINITY_WITNESS_NOT_OBSERVED`
rejected one *population-shared* linear direction (explained fraction
`0.000709`); it did not test protein-specificity under a target-conditioned
model. Correct-vs-wrong-protein and correct-vs-ligand-only admission remain
open and must be measured inside the episodic design, not before it.

### O3 — Coefficient-sharing failure — **UNTESTED, the primary hypothesis**

`d = 0` versus `d = 1..5` has never been run. The CQ-R2 failure is *consistent
with* target-specific heterogeneity averaging to zero under a shared `w`, but
is not evidence for it.

### O4/O5/O6 — Support association, query observability, statistic sufficiency

All untested. They are downstream of O1 and O3 and cannot be diagnosed until an
adequately powered evaluation panel exists.

## 3. Theory-to-model memo

### What the theory generates, and what it does not

The scientific spine is:

```text
source tasks -> shared function family -> support section -> query-valid state
             -> law-valued output
```

The load-bearing mathematical fact for the section is elementary and must be
labelled as such:

```text
a_t = M_S^T (M_S M_S^T + lambda I)^{-1} r_S   =>   a_t in row(M_S)
=>  dim(identifiable target-specific adaptation) <= rank(M_S) <= k
```

**Provenance, stated exactly:**

| claim | source |
|---|---|
| `a_t in row(M_S)`, rank bound `<= k` | linear algebra, not a theorem of `FINAL_FROZEN_THEORY` |
| positive ridge `lambda > 0` gives existence/uniqueness | strong convexity, consistent with the frozen archive's ridge role but applied here to the *section*, not to `g_mu^star` |
| `A(F,z) = K(B(z)F(z))`, simplex, band, mesh | `FINAL_FROZEN_THEORY`, unchanged |
| `d <= 5`, `k <= 5`, coverage thresholds | engineering choices constrained by support size |
| target-coefficient heterogeneity exists | **empirical hypothesis, untested** |

The frozen archive explicitly does not provide ranking guarantees, so no
ranking metric may be claimed as theorem-backed.

### Why this is not MAML, not ridge regression, and not an uncertainty head

- **Not MAML**: there is no free inner loop. Task freedom is bounded a priori by
  `rank(M_S) <= k`, so the adaptation cannot exceed what the support observes.
- **Not ordinary ridge**: the *family* `U` is learned across source tasks so
  that `k <= 5` observations suffice to section it. Ridge on a fixed basis has
  no cross-task outer objective.
- **Not an uncertainty head**: rank, conditioning and coverage are properties of
  the same section geometry `P_S = M_S^T (M_S M_S^T)^+ M_S`, not a bolted-on
  predictor. `c_q = ||P_S m_q||^2 / (||m_q||^2 + eps)` decides where adaptation
  is admissible, and off-coverage queries receive zero correction rather than a
  confident extrapolation.

### Traceability

| object | role | source | Gate |
|---|---|---|---|
| `phi(P,L)` 288D | observable coordinate | engineering input | biological admission (O2) |
| `U`, `d<=5` | shared family geometry | meta-learned | family Gate (O3) |
| `a_t` section | task identification | linear algebra | support controls (O4) |
| rank/coverage | observability certificate | section geometry | coverage Gate (O5) |
| `z` (4-5 bounded) | admitted statistic | frozen interface | conditional sufficiency (O6) |
| `K(B(z)F(z))` | output | frozen theory | unchanged |

## 4. Next decision

Branch selected by evidence: **expand the evaluation data**, by rebuilding the
corpus for the few-shot estimand rather than inheriting the quotient filter.

The next stage must be separately preregistered and must:

1. rebuild from the same governed projection under few-shot admission rules
   (single chain, Ki only, canonical ligand identity, scaffold present);
2. apply the *same* strict conflict closure already used by CQ — document,
   protein-40 homology and Murcko scaffold union — with no relaxation;
3. report independent component counts and the giant-component share **before**
   any split is frozen, and fail closed if component depth is still
   insufficient;
4. only then re-run the Phase 0 power gate with thresholds unchanged
   (`>=30` evaluation targets at `k=5`, `MDE_d <= 0.600`).

No model is authorized until that gate passes. O2-O6 remain untested, and
target-coefficient heterogeneity is neither supported nor refuted by anything
run so far.

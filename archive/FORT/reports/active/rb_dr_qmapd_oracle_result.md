# RB-DR-QMAPD — Stage O1 oracle result (2026-07-26)

Preregistration `reports/active/rb_dr_qmapd_oracle_preregistration.md`; design
`rb_dr_qmapd_mathematical_design.md`; feasibility `rb_dr_qmapd_power_audit.md`. Frozen Teacher A
(validated exact-Cholesky Bayesian residual posterior), single seed 1729, CUDA (`drug`). Train + spent
development rows only; no confirmation/Davis/sealed label. Runner `research/rb_dr_qmapd_o1.py`; data
`reports/active/rb_dr_qmapd_oracle_result.json`. This is `ARCHITECTURE_MECHANISM_RESULT_ONLY`.

## Decisive question

Does additional same-target evidence (k → K) contain a reproducible, scaffold-generalising SAR ranking
signal beyond the same teacher at k shots?

## Result — primary pair (k=4, K=16), 76 independent components, median 24 query ligands

| quantity | value | grouped 95% CI | threshold |
|---|---:|---|---|
| ρ(B0) | 0.3563 | — | — |
| ρ(T_4) | 0.3638 | — | — |
| ρ(T_16) | 0.3794 | — | — |
| **Δ_info = ρ(T_16)−ρ(T_4)** | **+0.0154** | [−0.0155, +0.0464] | ≥ max(0.03, MDE80 0.0452) |
| Δ_arch = ρ(T_4)−ρ(B0) | +0.0108 | [−0.0158, +0.0405] | — |
| **Δ_total = ρ(T_16)−ρ(B0)** | **+0.0262** | [−0.0079, +0.0621] | ≥ max(0.05, MDE80 0.0498) |

Identity check `Δ_total − Δ_info − Δ_arch = 0.0` (exact, identical rows). Empirical paired MDE80_info
`0.0452` (paired SD 0.1371/component), MDE80_total `0.0498`. RMSE(T_16) ≤ 1.02·min(RMSE) holds.

Secondary robustness pair (k=4, K=32), 40 components, median 17 query ligands (larger jump):
ρ(B0)=0.4025, ρ(T_4)=0.3868, ρ(T_32)=0.4211; Δ_info=+0.0333 [−0.0129, +0.0768]; Δ_arch=−0.0132
(k=4 support *hurts* on this ligand-rich subpopulation, recovered by k=32); Δ_total=+0.0186.

## Verdict

```
RB_DR_QMAPD_ORACLE_INFORMATION_FAIL_STOP
```

Frozen criteria 1–4 fail: Δ_info (+0.0154) is below its threshold (0.0452) with LCB below zero, and
Δ_total (+0.0262) is below 0.05 with LCB below zero. A single failed mandatory criterion stops the
program. (Criterion 5 is also technically unmet only because `recent_half` is undefined on a dateless
panel; that is not the decisive failure — 1–4 are.)

## Honest interpretation (mechanism, not null)

The signal is **positive and monotone in direction** — more same-target support consistently improves
within-target ranking (T_16 > T_4 > B0 primary; T_32 − T_4 = +0.0333 for the largest jump) — and the
oracle decomposition is exact. So this is not a flat null: there is a *small* amount of extractable
same-target SAR headroom. But it is **below the practically meaningful minimum effect and not
statistically resolvable** on this substrate: the per-component paired SD (0.137) with 40–76
independent components gives MDE80 ≈ 0.045–0.050, larger than the observed effect. Consistent with the
program's recurring finding, the binding constraint is component-level **power / label resolution**
(Metz pKi is rounded to 0.1 pK with 81.7% within-target duplicate values), not the student
architecture.

Per Section 7: Δ_arch is weakly positive while Δ_info fails, so the correct action is to **report the
k-shot baseline and not implement privileged completion distillation.** The student operator, teachers
B/C, completion designs, the Rao–Blackwellized law, the design-robust objective, and the one-seed gate
are **not built** (they are authorised only after O1 passes). No threshold, teacher, feature, or
hyperparameter was changed after the result; no data/model was added as a rescue. Confirmation, Davis,
and sealed labels were not read.

## What would change the answer (not executed here)

Only a more powered substrate can resolve a ~0.015–0.033 headroom: more independent endpoint-consistent
homology components (an additional open dense panel), or a continuous endpoint with finer label
resolution than 0.1-pK-rounded pKi. That is a measurement/data step, not an architecture step, and it
is out of scope for this fail-stopped gate.

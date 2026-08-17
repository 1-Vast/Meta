# Stage R14 screening: the regression-compatible ranking term is inert — claim withdrawn

> **Correction (2026-08-16, repaired and verified).** The dispersion artifacts
> in this directory originally carried a
> `"seal": "physical: QPSMPData include_meta_test=False"` claim that was false —
> the flag was never passed and the class default admitted the sealed split. No
> number here is affected, and that is now demonstrated: re-running the audit
> under the repaired fail-closed seal reproduces **105 of 105 numeric A0 fields
> bit-identically**. The seal blocks have been corrected in place with no
> numeric field touched. Full scope and verification:
> `tools/research/a2_readiness_v2/GOVERNANCE_INCIDENT.md`.

Numerical authority: `DISPERSION_R14_3seed.json` (+ `.rows.jsonl`),
`CONTRAST_R14_k0.json`, `ALIGNMENT_OPERATING_POINT.json`, and the nine per-run
`RESULT.json` files. Population: double-cold `meta_val`, 41 targets / 19
components, three seeds, 1200 steps, matched to the incumbent configuration.
`meta_test` never read.

Gates were preregistered in
`../stageR14_diagnostics_20260816/PREREGISTRATION.md` before any run.

## Arms

One changed variable, the *form* of the existing within-target ranking term
at its existing weight 0.5. The incumbent already carries RankNet at that
weight, so it is the matched misaligned control rather than a ranking-free
one.

| arm | ranking term |
|---|---|
| A0frozen | the retained R3R4 incumbent checkpoints (RankNet @ 0.5) |
| A0repro | the same configuration, **retrained** here (RankNet @ 0.5) |
| R1listce | regression-compatible ListCE @ 0.5 — the candidate |
| R3norank | `ranking_loss_weight = 0.0` — the **necessity control** |

## Outcome (k=0, three-seed, equal-component weighting)

| arm | MSE | **ordering floor** | `r` | shape | calib | CI | sd_p |
|---|---:|---:|---:|---:|---:|---:|---:|
| A0frozen | 2.1488 | 0.6922 | **0.213** | 0.9130 | 1.2358 | **0.580** | 0.132 |
| A0repro | 2.0911 | **0.6850** | 0.162 | 0.9032 | 1.1879 | 0.552 | 0.147 |
| R1listce | **2.0700** | 0.6931 | 0.186 | 0.9005 | **1.1695** | 0.544 | 0.124 |
| R3norank | 2.0924 | 0.6951 | 0.178 | 0.9049 | 1.1875 | 0.543 | 0.121 |

Paired component bootstrap against A0repro (positive = better):

| contrast | ordering floor | `r` |
|---|---|---|
| A0frozen − A0repro *(same config!)* | −0.0072 [−0.0716, +0.0648] | +0.0503 [−0.0461, +0.1509] |
| R1listce − A0repro | −0.0081 [−0.0300, +0.0137] | +0.0231 [−0.0300, +0.0750] |
| R3norank − A0repro | −0.0101 [−0.0250, +0.0038] | +0.0153 [−0.0242, +0.0605] |

## Gate verdicts

| gate | requirement | measured | verdict |
|---|---|---|---|
| **O1** (primary) | R1's k=0 `r`/ordering floor beats A0 with a positive component lower bound | floor 0.6931, **worse** than both A0repro (0.6850) and A0frozen (0.6922); `r` 0.186 < 0.213; interval crosses zero | **FAIL** |
| O2 | calibration no worse than 1.236 + 0.02 | 1.1695 | pass |
| **O3** | k=0 MSE improves **and** CI ≥ 0.570 | MSE 2.0700 (better) but CI 0.544 | **FAIL** |
| **O4** (necessity) | R1's `r` exceeds R3's with a positive lower bound | floor 0.6931 vs 0.6951 (Δ 0.002); `r` 0.186 vs 0.178 — both far inside the same-config noise band | **FAIL** |
| O5 | direction holds in 3/3 seeds | not reached | — |

**Three gates fail, including the primary and the necessity control. Under the
preregistered rule the core-innovation claim is withdrawn.** The family is
stopped here rather than advanced to a formal run.

## Why it failed, measured rather than inferred

The preregistration listed "the alignment identity does not hold" as the first
failure mode. It holds — that was verified before implementation. The failure
is the one nobody preregistered, and it is a property of the construction
itself.

`ALIGNMENT_OPERATING_POINT.json` measures how much gradient each term supplies
**at the model's actual operating point** (`r ≈ 0.2`, `sd_p/sd_y ≈ 0.2`, the
values Phase 2 measured), not at the optimum:

| term | mean ‖∂L/∂p‖ at the operating point |
|---|---:|
| squared error | 3.93e-01 |
| hinge margin | 2.28e-01 |
| RankNet | 1.51e-01 |
| **regression-compatible ListCE** | **6.53e-03** |

ListCE supplies **4.3% of RankNet's gradient and 1.7% of the regression
term's**. At weight 0.5 it is inert, and R3norank — simply deleting the
ranking term — reproduces it to within 0.002 on every metric. Two independent
lines of evidence agree.

The mechanism is the alignment property itself. The gradient is

    ∂L/∂p_k  ∝  w_k/T(p_k) − Σ_j w_j / Σ_j T(p_j)

so its scale goes as `1/T(p)`. Exactly the `1/T(p)` factor that makes the term
vanish at `T(p) ∝ w` also damps it everywhere else, and the damping is worst
precisely where this model lives — under-dispersed predictions make `T(p)`
nearly constant across the panel, so the bracket nearly cancels term by term.

**The general lesson: exact regression-compatibility and useful ranking
pressure are in tension.** A ranking term that is stationary at the regression
optimum and scale-free is necessarily weak wherever the model is
under-dispersed. This is a real constraint on the RCR-style construction that
its original setting (binary relevance, well-dispersed scores) does not
expose, and it was not visible in the pre-implementation alignment check,
which only probed `s = y`.

## A second finding: the frontier is close to retraining noise

`A0repro` is the incumbent configuration retrained here. Against the frozen
checkpoints it differs by **0.058 in k=0 MSE (2.0911 vs 2.1488) and 0.051 in
`r` (0.162 vs 0.213)** — same architecture, same configuration, same seeds,
different run. CUDA training here is documented as not bitwise deterministic.

That noise is comparable to the entire measured frontier: B3's best-on-record
2.055 sits **0.094** below A0's 2.149, only 1.6× the same-config spread. The
k=0 Pareto frontier in `report/BOUNDARY_20260816.md` should be read with that
in mind — its three points are separated by roughly one to two retraining
standard deviations, and none of its MSE differences was ever resolved by a
component bootstrap.

This also means **A0frozen's recorded 2.149 is one draw, not a constant**, and
comparisons of a newly trained arm against it inherit that variance.

## What is closed, and what is not

**Closed.** The loss-form axis for this family. R9 (cliff dose), R10
(variance), R12 (margin form) and now R14 (regression-compatible listwise)
have each varied the within-target ranking objective and none moved the
ordering floor. Combined with the Phase 2 finding that every ranking-primary
arm is worse than the regression-dominant incumbent, the evidence is that
**within-target ranking terms are not the lever for `r` on this data.**

**Not closed, and not tested here.** Whether a *larger* ListCE weight or a
shift closer to the label range would supply usable gradient. The shift was
fixed at 2.0 pK in advance and the preregistration forbade a post-hoc sweep;
tuning it now against a failed run would be exactly the selection bias Stage
R0 measured at 0.468 MSE. If that direction is ever reopened it needs its own
preregistration, and it starts from a weak prior given the four failed
loss-axis stages.

**The next evidence-supported hypothesis** is the one Phase 2 left standing:
`r` is bounded by what the ligand representation carries about within-target
potency ordering, not by the objective. Testing that requires a
representation-side probe — how much within-target ordering is extractable
from the trunk's ligand embeddings under a `meta_train`-only fit — which is a
separate, unstarted line and a diagnostic, not a model change.

`meta_test` labels were used for no fitting, selection or reported metric. A historical process-isolation incident remains open: some processes parsed and indexed the split.

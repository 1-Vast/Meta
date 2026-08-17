# Stage U — protein × local chemical transformation interaction: **NEGATIVE at U0**

Authorities: `PREREGISTRATION.md` (SHA-256
`fdc0a830aa92882d07b9aea50f22a4c72fc6d93f92c55a3be6bc15cd6a645c11`, recorded
inside `U0_CENSUS.json` and `U0_DECISION.json`), `U0_RELIABILITY.json`,
`U0_CENSUS.json`, `U0_DECISION.json`, `U0_OBSERVATIONS.jsonl.gz`,
`U0_PROVENANCE_meta_train.jsonl.gz`.
Structural verification: `tests/test_structural.py`, **22 passed**
(`RUN_SLOW=1`).

**Verdict: the frozen U0 admission gate fails. Stage U stops before U1 and
before any neural model.** No threshold was moved after the result; the
preregistration is unchanged. `model/` and production `scripts/` were not
modified by this stage, and the sealed confirmation split was not mounted.

---

## 1. What was frozen before any result

The estimand is the crossed double difference

    tau      = (shared core, R_a -> R_b, attachment environment,
                stereochemistry, charge change)
    delta_y(t,tau) = mu_tau + delta(t,tau) + noise
    D(tau,t1,t2)   = delta_y(t1,tau) - delta_y(t2,tau)

with `D` formed only within one core-inclusive exact transformation key. U0 was
preregistered to audit label aggregation, same-panel delta-pK reliability,
RDKit-based true MMP construction, the bipartite evidence graph, effective rank
and deployment coverage — and then to apply a frozen admission gate. U1
(interaction variance versus measurement noise) and U2 (the
transformation-conditioned local protein-region operator) were gated on U0.

## 2. U0 reliability — measured, with its limits stated

Both pre-aggregation artifacts hash-match the corpus manifest; **5,983 / 5,983
(100%)** of `meta_train` source rows were recovered. The corpus aggregation rule
was verified verbatim: within-panel median, then equal-panel median across
panels, endpoint exact positive uncensored Ki, `pKi = 9 - log10(Ki[nM])`.

| level | groups | identifiable | result |
|---|---:|---|---|
| L1 same panel, same protocol | 99 (98 disagreeing) | yes | residual variance 0.429 pK² (sd 0.655) |
| L2 same panel, different protocol | 0 | **no** | not identifiable; no number invented |
| L3 across panels | 133 (54 disagreeing) | yes | residual variance 0.182 pK² (pooled, duplication-deflated) |

Derived same-panel delta-pK reliability: point **0.858 pK²**, clustered
bootstrap over repeated-measure groups 95% CI **[0.688, 1.040]** pK².
Cross-panel: **1.221 pK²**. The repeated-measure subset is 4.1% of cells and
selected; these are supervision-reliability estimates, **not a benchmark MSE
floor**, and are labeled as such in the artifact.

## 3. U0 census — the graph is rich on one axis and dominated on another

Built with `rdkit.Chem.rdMMPA.FragmentMol` (Hussain–Rea, single cut), isomeric
SMILES, core-inclusive SHA-256 exact keys and a coarse key that strips stereo
and reduces attachment context. The bank contains **44,233** total
observations; the primary same-panel fit bank has **37,945** observations /
**243 targets** / **187 components** / **30,463 exact keys** (30,461 coarse).

| diagnostic | fit same-panel | internal same-panel |
|---|---:|---:|
| observations | 37,945 | 4,589 |
| targets / components | 243 / 187 | 34 / 25 |
| exact keys | 30,463 | 4,046 |
| keys spanning >=3 targets and >=3 components | **1,001** | 0 |
| top-1 key share | 0.05% | 0.07% |
| top-10 key share | 0.43% | 0.50% |
| graph connected components | 140 | 23 |
| largest component share | 36.6% | 30.1% |
| incidence numerical rank | 215 | 34 |
| incidence stable rank | 69.1 | 18.8 |
| incidence condition number | 336.8 | 34.7 |
| exact-key reuse fit→internal | 11.53% | — |
| core reuse fit→internal | 26.82% | — |

These are **empirical sufficient-richness diagnostics**, not a proof of
persistent excitation.

Deployment coverage on the frozen nested episode banks (query labels not read):
`C_k` exact = **0.226 / 0.362 / 0.442 / 0.526** at k=1/2/3/5. MMP is therefore
at best a training signal and a partial inference mechanism, never a universal
reference-based mechanism.

## 4. The frozen admission gate — 9 of 11 checks pass, 2 fail

| check | threshold | measured | verdict |
|---|---:|---:|---|
| same-panel fit observations | >= 2,000 | 37,945 | PASS |
| fit targets | >= 50 | 243 | PASS |
| exact keys spanning >=3 targets and >=3 components | >= 30 | 1,001 | PASS |
| internal same-panel observations | >= 300 | 4,589 | PASS |
| internal components | >= 10 | 25 | PASS |
| top-1 key observation share | <= 0.05 | 0.0005 | PASS |
| top-10 key observation share | <= 0.20 | 0.0043 | PASS |
| **top-1 target observation share** | **<= 0.25** | **0.2963** | **FAIL** |
| top-5 target observation share | <= 0.75 | 0.4726 | PASS |
| **top-1 component observation share** | **<= 0.25** | **0.2963** | **FAIL** |
| top-5 component observation share | <= 0.75 | 0.5116 | PASS |

The two failures are the same fact seen twice: one target, in a one-target
CD-HIT40 component, contributes **11,243 / 37,945 (29.63%)** of same-panel fit
MMP observations from 86 ligands. The next target contributes 1,954 (5.1%).
Transformation-key concentration is not the problem — the top key is 0.05% of
observations — but target/component concentration is, and it is exactly the
failure mode the degree-concentration gate exists to forbid: a few
high-throughput targets would carry the fit bank.

## 5. Decision

Per the frozen stop rule (“if any threshold fails: stop immediately, write a
negative report, train no neural model”):

- **U1 was not run.**
- **U2 was not implemented or trained.**
- **No neural model, checkpoint or trained arm exists in this stage.**
- **No threshold was moved and no unpreregistered module was added.**

## 6. Claim boundary — stated plainly

This is a **negative admission decision**, not a biological falsification. It
establishes that, under the frozen Stage U protocol, the core-inclusive
transformation graph on the governed BindingDB-Ki `meta_train` partition does
not satisfy the preregistered sufficient-balance condition, so the
protein × transformation interaction was never estimated. It does **not**
establish that such interaction signal is absent, and it must not be described
as evidence of absence.

A concurrent repository audit (`tools/research/stageV_core_mmp/`) independently
classified this Stage U design as an adaptive correction of Stage T and
superseded it with a corrected preregistration; that document is left untouched
here. The Stage U files remain read-only evidence of the frozen attempt.

## 7. Verification

* `RUN_SLOW=1 pytest tools/research/stageU_mmp_interaction/tests -q`:
  **22 passed**, including MMP canonical direction/inverse, attachment/stereo/
  charge preservation, physical meta-test seal, within-target same-core bank
  checks, no cross-population leakage, PYTHONHASHSEED-independence of the bank
  digest, no `hash()`, and the assertion that no neural model was trained before
  the U0 gate.
* Environment: conda env `drug`, Python 3.11.15, torch 2.6.0+cu124,
  CUDA available, RDKit 2023.09.6, numpy 1.26.4, scipy 1.17.1.
* Git commit at execution: `5bb373668969985ada9f1e208fc3fcfe886b3123`.

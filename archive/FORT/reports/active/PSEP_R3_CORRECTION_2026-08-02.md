# CORRECTION — R3 reclassified · programme status reset

Date 2026-08-02 · supersedes the R3 addendum in
`PSEP_REPRESENTATION_R1_R2_DECISION_2026-08-02.md` §8

**Prior claim (WITHDRAWN):** "the within-document ranking objective *is* the
invariance mechanism … real, reproducible, mathematically derived, ablated."

**Corrected status: `OBJECTIVE_EFFECT_UNRESOLVED_DUE_TO_CROSS_TASK_PAIRING`.**

**Programme status: `NO CORE MECHANISM IDENTIFIED`.**

All R1/R2/R3 artifacts are preserved unmodified. Nothing is reinterpreted in
place; this document records what they do and do not support.

---

## 1. Independent audit of the review's claims

Every claim was re-derived from the artifacts before being accepted. All are
confirmed; the central one is **worse** than reported.

### 1.1 The ranking loss is predominantly cross-target — CONFIRMED, understated

Training groups are keyed on `(fold, document)` only
(`psep_representation.py:246`, `psep_seeds.py:93`). A ChEMBL document routinely
reports one compound series against several targets, so a "within-document" pair
is not a within-target pair.

| measured | value |
|---|---:|
| training document groups | 12 129 |
| total training pairs | 249 122 844 |
| **cross-target-unit pairs** | **226 483 925 = 90.9 %** |
| cross-accession pairs | 226 336 091 = 90.9 % |
| largest group | 12 492 rows spanning **33 units** |
| **share of all pairs held by the top 10 groups** | **97.0 %** |

The review said 90.9 %; it did not note that **97 % of every gradient came from
ten documents**. The model is ligand-only, so cross-target pairs cannot teach
target-specific reordering. The measured effect is consistent with global ligand
ranking and document-composition learning, and is dominated by a handful of large
papers.

### 1.2 Assay nuisance is not cancelled — CONFIRMED

Evaluation pairs require a shared document, not a shared assay
(`psep_transfer.py:79`). Measured: **171 705 / 1 551 246 = 11.1 %** of evaluation
within-document pairs cross assays (review: 12.3 %; the gap is pair-capping).

The identity `(y_i − m_d) − (y_j − m_d) = y_i − y_j` removes only a nuisance
*shared by both measurements*. Where assays differ, an assay-specific offset
survives. The correct invariant context is at least `(unit, document, assay)`.

### 1.3 The sealed-role firewall is declarative, not operational — CONFIRMED

`psep_substrate.py:343` assigns roles; `psep_substrate.py:361` then fits the
frozen base on **all** activities. The base therefore saw `discover`,
`validate` and `confirm` rows (cross-fitted by component, but not by role).

- `rank − mse` is largely unaffected — the base cancels from that contrast.
- **Every "gain vs base" number in this programme is not role-pure.**
- `validate` and `confirm` are **development-contaminated for base-dependent
  claims** and are no longer pristine confirmation sets.
- `"validate_read": false` in the payloads is a hard-coded literal, not a proof.

### 1.4 Seed uncertainty omitted — CONFIRMED

`psep_seeds.py:148` averages the five observed seeds per component, then
bootstraps components. The reported `[+0.0013, +0.0180]` is *conditional on those
seeds*. A two-way seed × component bootstrap gives **+0.0097 [−0.00033,
+0.01938]** — crossing zero. Sign test on 5/5 positive: p = 0.0625. Endpoint
intervals all cross zero.

### 1.5 Optimisation budget confounded — CONFIRMED

Measured per fold: MSE **36–42** optimiser steps/epoch vs ranking **150** — a
**3.6–4.2× per-epoch** advantage before early stopping (~1.7× after). Parameter
count was matched; optimisation and sampling budget were not. Pair exposure also
grows quadratically with document size, compounding §1.1.

### 1.6 `rank_centred ≡ rank_raw` is a sanity identity, not a causal ablation — CONFIRMED

It proves centring cancels inside document pairs. It does **not** show that
nuisance cancellation *caused* the gain: the ranking arm simultaneously changes
label geometry, pair sampling weights, optimiser steps, and magnitude handling.

Prior art also weakens the novelty framing: **ActFound** (Nat. Mach. Intell.
2024) already combines within-assay pairwise differences with meta-learning for
incompatible assay scales.

### 1.7 R3 is not meta-learning — CONFIRMED

No support set, no `k`, no target-conditioned state, no protein input, no
adaptation operator. Predictions do not change when target support changes. Even
a fully replicated effect would be a **support-free training objective**, not a
few-shot adaptation mechanism.

---

## 2. What the artifacts still support

- Fixed folds reproduce **exactly** (seed 20260802 reproduced R2 bit-for-bit).
- All 541 unit × seed cells complete; components are the bootstrap unit.
- Five seed estimates directionally consistent.
- The registered MDE was never relaxed, and the failing verdict was reported.
- The **negative** results (adaptation, operators, routing) rest on wrong-support
  / wrong-target / document-split controls that are *not* affected by §1.1–§1.5,
  because those controls compare arms under identical pairing and budget.

## 3. What is withdrawn

- "The objective is the invariance mechanism" — the derivation is sound in the
  abstract but **does not describe what was actually trained** (§1.1).
- "+0.0097 real and reproducible" — downgraded to unresolved (§1.4).
- Any base-referenced magnitude as role-pure (§1.3).

## 4. Directed next steps (adopted)

1. **Freeze R3.** Reclassified; artifacts preserved.
2. **PSEP v2, role-safe.** Physically separate role files; fit preprocessing and
   the OOF base on `discover` only; tests that fail on any validate/confirm
   affinity read. Treat current validate/confirm as contaminated for
   base-dependent claims; reserve a *new* external confirmation source.
3. **One corrected discovery seed.** Pair keys `(component, unit, document,
   assay)`; assert zero cross-unit and cross-assay pairs; cap pairs per context
   and weight contexts equally. Equal-budget arms: pointwise MSE ·
   within-context pairwise-difference MSE · within-context sign-ranking ·
   context-shuffled and cross-target-pair controls.
4. **Strict gates**: lower95 > 0.005, survives removal of the largest documents,
   consistent pIC50/pKi direction, beats budget-matched MSE *and*
   pairwise-difference MSE. Null ⇒ stop the branch.
5. **Then** multi-seed with a **two-way** seed × component bootstrap.
6. **Then** exactly one episodic operator on within-assay pairwise sufficient
   statistics, with the k-specific / k=1-no-op / wrong-support / query-dependence
   / ablation requirements.
7. **`validate` stays closed** until that gate passes.

## 5. Lesson

I derived the pairwise-cancellation identity correctly and then failed to check
that the code implemented it. The identity holds for within-*document-and-target*
pairs; the implementation grouped by document alone, so 90.9 % of pairs violated
its premise. **A derivation is not evidence that the thing derived is the thing
that ran** — the check is one `groupby` on the training pairs, and it was never
run until review forced it.

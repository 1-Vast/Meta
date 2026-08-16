# PSEP — abundant-to-scarce transfer: decision report

Date 2026-08-02 · Runners `research/psep_{transfer,source_select,oracle_stability}.py`
Artifacts `reports/active/psep_{transfer,source_select,oracle_stability}_2026-08-02.json`
Role read: **`discover` only**. `validate` and `confirm` never opened.

**Verdict: `HEADROOM_IS_A_DOCUMENT_LOTTERY_NOT_TARGET_STRUCTURE`.**
The abundant-to-scarce *routing/selection* reframe fails, for a reason that
unifies every negative result in this programme.

---

## 1. Information-theoretic diagnosis: why k-shot adaptation failed

The target-specific object is a vector in `R^266` whose task-population
covariance has **participation ratio 115** (M2/M3). Each support label supplies
one noisy linear functional of it, contaminated by a per-context offset carrying
**68 % of residual variance** — a nuisance about twice the signal.

k labels identify at most a k-dimensional projection of a ~115-dimensional
object. Measured recovery tracks that bound: k=5 obtains 12 % of the full-head
gain, k=20 obtains 62 %.

**This motivated the reframe, and the counting argument was sound.** Choosing
among `N ≈ 550` pre-fitted source heads costs `log2(550) ≈ 9.1 bits`; five noisy
scalars can carry that, while estimating a 115-dimensional vector needs far more.
Selection is *information-feasible* where estimation is not.

**The argument was necessary but not sufficient.** Bits being available says
nothing about whether the channel is informative. It is not (§3, §4).

---

## 2. T1 — bounding the whole routing family

Leave-one-homology-component-out, 538 units / 379 components, scored on
pair-half B with every selective arm choosing on half A.

| arm | vs base | vs uniform | neg-transfer |
|---|---:|---:|---:|
| uniform (global multitask floor) | +0.0001 | +0.0000 | 0.41 |
| random_source | −0.0041 | −0.0042 | 0.53 |
| **universal_best** (best mean transfer on *other* components) | **+0.0118** [+0.0059,+0.0179] | +0.0117 | 0.41 |
| protein_top1 | +0.0027 | +0.0026 | 0.47 |
| chemistry_top1 | −0.0005 | −0.0006 | 0.49 |
| chemistry_softmax (τ=0.02) | +0.0058 [+0.0003,+0.0113] | +0.0057 | 0.46 |
| own head (~140 labels) | +0.0134 | +0.0133 | 0.39 |
| oracle_split (honest ceiling) | +0.1222 | +0.1221 | 0.03 |
| oracle_insample (**biased**) | +0.1698 | +0.1698 | 0.00 |

**Selection bias, measured and removed.** The first pass of this gate reported a
ceiling of +0.165 — ten times the target's own label-fitted head. That was the
maximum of ~550 noisy estimates scored on the pairs that chose them; the bias of a
max over N draws (`≈ σ√(2 ln N)`) reproduced essentially the whole effect.
Split-sample selection removed **+0.0476** of it.

**No observable predicts which source transfers.** Spearman between similarity and
realised transfer: chemistry **−0.048** [−0.062,−0.034] (wrong sign), protein
**+0.002** [−0.003,+0.007] (exactly zero).

---

## 3. T2 — can k labels select the source?

Selection by centred support MSE over pre-fitted source heads. Gains vs uniform.

| k | support_select | wrong_support | **select − wrong** | support-fit rank ρ |
|---|---:|---:|---:|---:|
| 3 | −0.0123 | −0.0127 | **+0.0004** [−0.0029,+0.0038] | +0.023 |
| 5 | −0.0113 | −0.0118 | **+0.0005** [−0.0038,+0.0051] | +0.029 |
| 10 | −0.0083 | −0.0076 | **−0.0007** [−0.0072,+0.0060] | +0.039 |
| 20 | −0.0047 | −0.0061 | **+0.0013** [−0.0056,+0.0089] | +0.049 |

**Selection driven by a *wrong target's* support performs identically at every k.**
The criterion carries no target-specific information. This is the same failure
signature as the `cnp` operator: the mechanism runs, produces numbers, and is
indifferent to whose data it receives.

Support fit barely predicts transfer (ρ ≈ 0.03), and hard selection on a
near-uncorrelated criterion over 550 candidates is *worse than not selecting*
(−0.0113 at k=5); top-5 averaging recovers most of the loss (−0.0029); plain
uniform beats both. Winner's curse, exactly as expected when the ranking signal
is noise.

**Correction to this table's `ridge_k` column (not shown above).** In the run,
`ridge_k` was granted the best of three ridge values scored *per episode on the
evaluation half* — eval-label selection, the same bias class removed from the
oracle. Its printed values (+0.0121 … +0.0302) are therefore inflated; the honest
figure from M2 is **+0.0019 at k=5 over base**. No conclusion here depends on it,
because `support_select` is negative and target-indifferent regardless.

---

## 4. T3 — the decisive test: is the ceiling real?

T1's pair-halves share query *rows and documents*, so a head chosen on half A is
chosen partly for fitting those particular molecules in those particular
documents. T3 splits the query rows **by document**: choose the source on
documents half 1, score on documents half 2 — same target, both halves disjoint
from support.

157 units / 127 components (units need ≥ 4 query documents).

| quantity | value |
|---|---:|
| oracle_insample_half2 (**biased**) | **+0.1378** [+0.1255,+0.1512] |
| **oracle_doc_split** (honest, exploitable) | **−0.0053** [−0.0192,+0.0077] |
| document-specific component | **+0.1431** |
| rank correlation of source quality across documents | **+0.0518** [+0.0325,+0.0725] |
| top-10 overlap across documents | 0.036 (chance ≈ 0.018) |
| **own head, same document split** | **+0.0244** [+0.0113,+0.0383] |

**The entire apparent headroom is document-specific.** Selecting the best source
on one set of the target's documents and applying it to another set of the *same
target's* documents yields **nothing** (−0.0053, interval spanning zero). Source
quality does not replicate across documents: rank correlation +0.05, top-10
overlap 2× chance.

**The harness validates itself on the same split.** The target's *own* head
**does** replicate across documents (+0.0244 [+0.0113,+0.0383]). So the
measurement can detect cross-document transfer when it exists — and finds none
for source selection. This is a structural control, not a statistical one.

With ~550 candidate directions, some head aligns with whatever idiosyncratic
ordering a given document happens to carry. That is the whole of the +0.12–0.14.
**No router, learned or fixed, can reach it.**

---

## 5. Core mechanism gate

| # | Condition | Result |
|---|---|---|
| 1 | source population improves unseen target prediction | **only unconditionally** (universal_best +0.0118) |
| 2 | learned routing beats fixed similarity | **N/A — fixed similarity is already ≈0, and the target is a lottery** |
| 3 | negative transfer decreases | **no** — support_select raises it (0.59 vs 0.39 for uniform) |
| 4 | benefit remains under provenance separation | **no** — 100 % of the headroom is document-specific |
| 5 | mechanism contribution removed by ablation | **fails** — wrong-target support gives identical results |

**No routing/selection mechanism is admitted.** Per the standing STOP rules this
is not rescued by a learned router, a mixture-of-experts, or an abstention head:
the object they would route to does not survive a document split.

---

## 6. The unifying finding across the whole programme

Every apparent **target-conditional** effect on this substrate is
**document-specific**; every effect that survives provenance separation is
**unconditional**.

| effect | apparent | after provenance separation |
|---|---:|---:|
| chemical head, conventional split (D0) | +0.0756 | +0.0029 (93 % collapse) |
| document-mean oracle (D0) | +0.0860 | **exactly 0** (structural) |
| best-source routing (T1→T3) | +0.1378 | **−0.0053** (100 % collapse) |
| CNP operator support-conditioning | +0.0274 | **0** (wrong support identical) |
| — | | |
| target's own head, ~140 labels (M4/T3) | | **+0.0134 … +0.0244** ✓ |
| universal_best source head (T1) | | **+0.0118** ✓ |
| target-agnostic nonlinear head (operator gate) | | **+0.0274** ✓ |

The two survivors are both *unconditional*. Nothing that requires knowing *which
target you are looking at* survives except the target's own labelled head — which
needs ~140 labels, not 5.

---

## 7. The correct meta-learning problem supported by the data

Not `support → adaptation`. Not `source population → routed transfer`.

> **Population-level multitask learning of a target-agnostic chemistry model,
> evaluated under simultaneous scaffold/document/assay separation.**

The evidence: a nonlinear target-agnostic head trained across source components
gains **+0.0274 [+0.0067]** on held-out components — larger than the entire
target-specific object (+0.0154), and ~14× what k=5 adaptation delivers.

The scientific contribution is the **measurement apparatus**, which is
transferable and which the field lacks: an 828-component provenance-separated
substrate plus five controls that each independently caught a false positive in
this programme — the document-mean oracle, the matched-capacity target-agnostic
head, the wrong-target/wrong-support derangement, `universal_best`, and the
document-split oracle.

---

## 8. Remaining risks and STOP conditions

- **`own` head at +0.0244 in T3 vs +0.0134 in T1.** T3's subset (157 units with
  ≥4 query documents) is deeper than average, so the object is larger in
  well-measured targets. This does **not** rescue routing (measured on the same
  subset) but it does suggest depth, not method, governs the effect.
- **Single seed** throughout discovery, by programme discipline. None of the
  conclusions turn on seed: the failing controls (wrong support identical,
  document-split oracle null) fail structurally.
- **pIC50-weighted corpus** (314/379 components). Endpoint splits in the JSONs
  show the same pattern.
- **`validate` and `confirm` remain sealed.** The one positive worth confirming
  (+0.0274 target-agnostic head) has **not** been confirmed on a sealed role.
  That is a deliberate hold: burning a sealed role is a one-way action and needs
  an explicit decision.

**STOP conditions now in force:**
1. The `support → adaptation` track is closed (operator gate, 3 families).
2. The `source routing/selection` track is closed (T1–T3).
3. No further few-shot mechanism may be proposed on this substrate without first
   passing a **document-split** replication of whatever object it claims to use.

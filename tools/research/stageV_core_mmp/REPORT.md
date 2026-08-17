# Stage V — core-inclusive MMP × protein interaction: **STOPPED before training**

**Verdict: the corrected Phase-1 test does NOT establish a transferable
protein-conditioned interaction signal. It also does not refute one. Three
independent frozen gates fail, and all three are statements about the
dataset's ability to identify the requested estimand — not about biology.**

Authorities: `PREREGISTRATION.md` (SHA-256 `c567f660…5844d4`),
`STAGE0_FORENSICS.json`, `STAGE_U_GOVERNANCE_AUDIT.md`, `V0_V1_RESULT.json`,
`ENVIRONMENT.json`, `tools/research/stageT_mmp/CORRECTION_20260817_CORE_KEY.md`.
Structural and leakage tests: `tests/`, **31 passed** (`RUN_SLOW=1`, 174 s).

No neural model was built or trained. Nothing promoted. The sealed split was
never mounted; the development-validation split was never read.

---

## 1. Stage 0 — the Stage T defect, quantified

Stage T's key was `f"{attachment_context}|{R_a}>>{R_b}"` — **no shared core** —
and it median-pooled multiple cores into one target effect while its descriptor
never read the core. The requested estimand required the core.

The decisive measurement holds the **protein fixed**: within one target, one
Stage T key realised on ≥2 cores.

| statistic | value |
|---|---:|
| across-core \|Δ delta_y\| median | **0.269 pK** |
| mean / p95 / max | 0.415 / **1.268** / 3.401 |
| Stage T fit `D` rows whose targets had **disjoint** core sets | **8,010 / 19,851 (40.4%)** |
| internal | **219 / 759 (28.9%)** |
| observations under a multi-core key (fit) | 7,130 (18.8%), max 57 cores per key |

Stage T's `D` truth has sd 0.804 pK, so the uncancelled nuisance is the same
order as the signal. **Stage T's "cancels `mu_tau` exactly" claim is false for
40% of its training rows**, and the defect supplies a mechanism for its own
inverted controls: fitting a target that is 40% contaminated is substantially
fitting noise, which is why its label-shuffled arm generalised better and the
zero predictor won.

**Consequences, separated as required:**

* **Valid Stage T result** — the concrete coarsened-key, pooled-protein
  discriminator is rejected. Its measurements stand.
* **Invalid Stage T claim** — "protein-conditioned SAR latent space is formally
  closed under the current BindingDB protocol." **Withdrawn.** It tested one
  coarsened representation and cannot close a family it never instantiated. The
  figure "1,112 rich keys" must not be reused; the core-inclusive count is 1,001.
* **Genuinely unmeasured** — the core-inclusive estimand with any protein
  operator; a local (non-pooled, non-target-ID) protein-region operator;
  stereochemical edits (1 fit / 0 internal) and charge-changing edits
  (326 fit / 2 internal).

## 2. Stage U — audited, stopped, superseded

Full audit in `STAGE_U_GOVERNANCE_AUDIT.md`. Timeline: Stage T's gate metrics
were read at **17:08**, Stage U's preregistration was frozen at **17:12**. Stage
U is therefore an **adaptive correction informed by Stage T, not an independent
confirmation**, and every claim inheriting from it carries that label.

Stage U's chemistry was right — core-inclusive key, core-consuming descriptor,
an interaction-variance gate before training, a local region operator, paired
substitution as the primary control. But four **load-bearing requirements were
missing** after metrics had been read:

1. no residue-token permutation control — fatal for a stage whose premise is
   *ordered* region tokens;
2. no capacity-matched random protein representation;
3. `fit_unsampled` was a random 10% of fit rows, retaining the same targets and
   keys, so it cannot detect target-key memorisation — the exact construction
   ruled out;
4. no identical-initialization / identical-minibatch-order requirement.

Plus two defects: gate 1 compared Pearson against a constant-zero predictor
(undefined), and the candidate had no nested zero predictor, so uncontrolled
output amplitude could lose to `A_zero` for reasons unrelated to protein
information.

Per the governing rule, Stage U is **stopped**, its preregistration is **not
edited**, and Stage V inherits **every Stage U threshold verbatim** (they were
frozen before any core-inclusive number existed) while adding the four missing
controls and repairing the two defects. **No threshold was loosened.**

## 3. V0 — core-inclusive census: **FAIL on concentration**

| check | measured | threshold | |
|---|---:|---:|---|
| same-panel fit observations | 37,945 | ≥ 2,000 | PASS |
| fit targets | 243 | ≥ 50 | PASS |
| rich exact keys (≥3 targets, ≥3 components) | **1,001** | ≥ 30 | PASS |
| internal same-panel observations | 4,589 | ≥ 300 | PASS |
| internal components | 25 | ≥ 10 | PASS |
| top-1 key share | 0.0005 | ≤ 0.05 | PASS |
| top-10 key share | 0.0043 | ≤ 0.20 | PASS |
| **top-1 target share** | **0.2963** | ≤ 0.25 | **FAIL** |
| top-5 target share | 0.4726 | ≤ 0.75 | PASS |
| **top-1 component share** | **0.2963** | ≤ 0.25 | **FAIL** |
| top-5 component share | 0.5116 | ≤ 0.75 | PASS |

One target, alone in its component, contributes **29.6%** of all same-panel fit
MMP observations. This independently reproduces the concurrent Stage U U0 stop.
It is an admission/balance failure, not evidence about biology.

### V0b — evaluability (rule inherited from Stage U: < 100 rows ⇒ not evaluable)

| surface | rows | components | keys | EIU | evaluable |
|---|---:|---:|---:|---:|---|
| fit (training) | 12,740 | 99 | 4,651 | 99 | yes |
| internal, **repeated keys** (Stage U's designated primary) | **32** | **4** | 30 | **4** | **NO** |
| internal, transformation-cold | 514 | 7 | 510 | 7 | yes |
| internal, all | 546 | 10 | 540 | 10 | yes |

**Internal rich exact keys (≥3 targets and ≥3 components): 0.**

Requiring a *complete* chemical context that repeats across targets **inside the
withheld protein components** leaves 32 rows over 4 components. The primary
evaluation surface for the corrected estimand does not exist on this dataset.

## 4. V1 — interaction variance vs supervision noise: **FAIL, resolved**

The question a neural model would answer only matters if there is between-target
variance to explain. Within one *complete* transformation, across targets:

| quantity | value | 95% interval |
|---|---:|---|
| `MS_effect` (pooled between-target mean square) | **0.4517** | [0.2104, 0.7509] |
| `sigma2_noise` (T0 same-panel difference variance) | 0.8576 | [0.6864, 1.0272] |
| **`theta = MS_effect - sigma2_noise`** | **−0.4059** | **[−0.6889, −0.0577]** |
| ratio `MS_effect / sigma2_noise` | **0.527** | |
| keys with ≥2 targets / effects / components | 4,651 / 12,133 / 99 | |

**The interval lies entirely below zero.** The internal partition agrees
descriptively: `MS_effect` 0.2610, ratio 0.304.

Read plainly: when two different proteins undergo the **identical** chemical
transformation on the **identical** scaffold, the spread in ΔpK is about half
what repeated measurements of a single protein–ligand pair already disagree by.

**The claim that survives without depending on the noise estimate at all:**
`MS_effect = 0.4517 pK²` is an *upper bound* on protein-attributable interaction
variance, because it still contains measurement noise. So
`Var(delta(t,tau)) ≤ 0.4517 pK²` (sd ≤ 0.67 pK) before subtracting any noise.

**Honest caveat on the noise reference.** T0 recorded that the technical-versus-
condition split is **not identifiable** (the assay protocol hash never splits a
panel). The L1 groups therefore mix technical replication with a paper reporting
two conditions, which likely **inflates** `sigma2_noise` and makes `theta` more
negative than a pure-measurement reference would. Two readings are admissible:
(a) there is no detectable interaction variance; (b) the noise reference is
overestimated. The frozen gate fails under either, and the upper bound above
holds under both.

### V1 sensitivity — the negative is not an artefact, but it is weaker where it matters

Gates had already fired; these are descriptive and were verified independently
of the concurrently recorded figures.

| subset | keys | effects | `MS_effect` | `theta` | 95% interval |
|---|---:|---:|---:|---:|---|
| all fit keys | 4,651 | 12,133 | 0.4517 | −0.4059 | [−0.6907, −0.0531] **resolved** |
| dominant target removed | 4,651 | 12,133 | **0.4517** | −0.4059 | [−0.6907, −0.0531] **resolved** |
| keys spanning ≥2 **components** | 2,616 | 6,892 | **0.6939** | −0.1637 | [−0.5085, **+0.1218**] *unresolved* |

Two things follow, and the second matters more than the first.

1. **The V0 concentration failure does not drive V1.** The dominant target
   (29.63% of observations) contributes only singleton keys, so removing it
   changes `MS_effect` not at all while top-1 target share falls to 7.32%. The
   two gate failures are independent.
2. **On the subset that actually bears on the question — transformations
   repeated across *different protein components* — the interval crosses zero.**
   `theta` = −0.164 [−0.509, +0.122]. So the strong resolved-negative statement
   holds on the pooled fit set but **not** on the cross-component restriction,
   where the evidence is *inconclusive rather than decisive*. `MS_effect` there
   is 0.6939, bounding protein-attributable interaction variance at sd ≤ 0.83 pK
   before any noise subtraction.

This is why the verdict below is **"not estimable"** and not "no interaction
exists". The pooled negative is partly carried by within-component key
repetitions, which are the least informative rows for a protein-cold question.

### V1 with a direct pair-level noise reference (post-hoc, descriptive)

`PAIR_LEVEL_NOISE_AUDIT.json` measures the supervision of the estimand itself:
repeated same-panel deltas of the *same* MMP pair, from the raw provenance.
Only **88 / 42,534** same-panel MMP pairs have both cells measured in more than
one shared panel, so the subset is reference-heavy and selected; **40 of the 88
are zero-range** (one physical measurement curated twice, not independent
replication).

| noise reference | point | pair-cluster 95% CI |
|---|---:|---:|
| pair-level, all 88 repeated pairs | 0.166 pK² | [0.098, 0.241] |
| pair-level, 48 disagreeing pairs only | 0.303 pK² | [0.200, 0.427] |
| preregistered T0 cell-level `2*L1` | 0.858 pK² | [0.686, 1.040] |

Re-running V1 against each reference (two-way key × component bootstrap):

| V1 contrast | T0 reference | pair-level all | pair-level disagreeing |
|---|---:|---:|---:|
| all keys `theta` [95% CI] | −0.406 [−0.704, −0.073] | +0.285 [−0.132, +0.434] | +0.148 [−0.488, +0.219] |
| cross-component keys `theta` [95% CI] | −0.164 [−0.504, +0.110] | +0.528 [+0.014, +0.606] | +0.391 [−0.327, +0.368] |

The only positive lower bound comes from the downward-biased all-group
reference. Under the conservative disagreeing-only reference no contrast has a
positive lower bound. The interaction variance is therefore **not identifiable
above the defensible noise envelope** on this corpus. This refines — but does
not overturn — the caveat above: the frozen V1 gate used the preregistered
reference and failed; the direct pair-level evidence says the noise reference
is likely lower, while the cross-component signal remains unresolved rather
than resolved.

### V0c — can the development-validation split supply the missing surface? No.

`METAVAL_STRUCTURE_CENSUS.json` (structure-only; no pK accessed or bound by the
module) counts the same-panel MMP relation on the development-validation split.
It carries **7,209 observations / 41 targets / 19 components / 4,968 exact
keys / 115 rich keys / 2,757 potential D rows**, so it is a plausible
transformation-cold surface. But the double-cold split forbids shared ligand
identities and scaffolds, and the measured exact-key overlap with `meta_train`
is **0**. It therefore cannot supply the repeated-key protein-cold surface
either. The Stage V primary surface is unsuppliable on both development
splits.

## 5. Decision, scoped exactly

Per the frozen stop rules:

* **V0 fails** → close only the **estimability of the exact-MMP route on this
  dataset**. A single target carrying 29.6% of the evidence is a support
  problem.
* **V0b leaves the primary surface not evaluable** (32 rows, 4 components; 0
  internal rich keys) → the requested estimand is **not identifiable** here.
* **V1 fails** → **stop before neural training.** The frozen gate required the
  2.5th percentile of `theta` to exceed zero; it is −0.689. On the pooled fit
  set the negative is resolved; on the cross-component subset it is negative but
  unresolved. Either way the gate fails and no larger network is a rescue for an
  estimand whose signal is bounded at sd ≤ 0.83 pK *before* any noise
  subtraction.

**Therefore: no V2 model was built.** The corrected Phase-1 test could not be
carried to the model stage, because the dataset cannot identify the estimand the
test requires.

### What may and may not be claimed

**May:** on the governed BindingDB-Ki double-cold corpus, with same-panel Ki
supervision and single-cut MMP transformations carrying full chemical context,
the protein × transformation interaction variance is not identifiable above
supervision noise, and no evaluation surface with adequate support exists on the
withheld protein components.

**May not:** that protein-conditioned interaction representations are
impossible, that protein biology does not modulate SAR, or that Phase 1 is
closed. **Insufficient support is not biological absence.** Stage T's global
closure claim remains withdrawn and is not reinstated by this stage.

## 6. Limitations

1. **No model was trained**, so nothing here speaks to any operator's capacity.
2. `sigma2_noise` rests on a small, selected, partly duplicate-inflated repeated-
   measure subset (T0), with the technical/condition split non-identifiable. It
   is **not** an MSE floor and is not quoted as one.
3. The V1 bootstrap resamples keys, components and noise groups; it does not
   capture training stochasticity, because there is no training.
4. Single-cut MMP with a full core is the strictest reasonable definition of
   "same transformation". A looser but still principled definition (for example
   core-similarity classes rather than core identity) would trade exactness for
   support and was **not** tested — it is a different estimand and would need
   its own preregistration.
5. Stereochemical and charge-changing edits remain unmeasured (1 and 2 internal
   observations).
6. This stage is an adaptive correction informed by Stage T. It is not an
   independent confirmation of anything.

## 7. What a future test would need

Not a bigger model. Either (a) a corpus where a complete transformation recurs
across many protein components — which BindingDB-Ki double-cold does not
provide, 0 internal rich keys — or (b) a preregistered, looser-but-principled
transformation equivalence class with its own cancellation analysis, since
loosening the key reintroduces exactly the `mu_tau` residual Stage 0 measured
at median 0.269 pK.

## 8. Verification

* `tests/`: **31 passed** (`RUN_SLOW=1`, 174 s) — 19 in `test_structural.py`
  plus 12 in `test_stage_v.py` from the concurrent Stage V work. The 19 include
  regression
  pins that Stage T's key omitted the core and that Stage V's key and descriptor
  include it; core-inclusive determinism; canonical direction and inverse
  negation; stereochemistry and charge preservation; double-difference algebra;
  no cross-target, cross-panel or fit/internal contamination; the physical
  meta-test seal; bank stability across `PYTHONHASHSEED` 0/1/12345; and
  parsed-AST checks for `hash()`, unseeded RNGs, label paths into the key, and
  the development-validation split name.
* Environment and artifact SHA-256 digests: `ENVIRONMENT.json`
  (conda `drug`, Python 3.11.15, torch 2.6.0+cu124, RDKit 2023.09.6,
  numpy 1.26.4, scipy 1.17.1, git `5bb3736`).

## 9. Commands

```bash
python -m tools.research.stageV_core_mmp.stage0_forensics
python -m tools.research.stageV_core_mmp.v0_census
RUN_SLOW=1 python -m pytest tools/research/stageV_core_mmp/tests -q
```

---

## 10. Phase 1 final decision pointer

`PHASE1_FINAL_DECISION.md` / `PHASE1_FINAL_DECISION.json` consolidate this stage
with Stage S, Stage P, Stage T and Stage U into the bounded Phase 1 verdict:
**BOUNDED NEGATIVE under the current BindingDB-Ki double-cold protocol** —
tested mechanisms closed as not estimable / negative; no biological
impossibility claim; MSA/coevolution externally blocked; Davis/KIBA
promotion-gated.

## 11. Synthetic calibration of the V1 statistic (post-hoc, descriptive)

`V1_SYNTHETIC_CALIBRATION.json` uses the exact chi-square null distribution on
the real 4,651-key / 12,133-effect graph (df = 7,482) to ask what interaction
size the observed `MS_effect = 0.4517` would imply under each noise reference.

* Under the preregistered T0 reference (0.858) a zero-interaction model
  predicts MS ≈ 0.858; the observed 0.4517 has empirical P ≥ 1.0 under that
  null — T0 is an upper-bound reference for the MMP-delta estimand.
* Under the direct pair-level references, if all excess above noise were
  signal, the implied interaction sd is **0.39 pK** (disagreeing-only
  reference 0.303) to **0.53 pK** (all repeated-pair reference 0.166).
* This is consistent with, and bounded by, the unresolved cross-component V1
  interval; it is a latent-effect size under an assumption, not evidence of a
  detectable signal, and it cannot reopen the frozen gate.

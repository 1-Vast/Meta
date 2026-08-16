# PSEP operator gate — results

Date 2026-08-02 · Runner `research/psep_operator.py` · Device CUDA (RTX 4060)
Artifacts `reports/active/psep_operator_2026-08-02.json`, `..._records_....parquet`
Pre-registration: `PSEP_OPERATOR_PREREGISTRATION_2026-08-02.md` (written before results were read)

Reported on **`meta_test`: 78 components / 114 units**, never seen in training or
model selection. `validate` and `confirm` roles never opened.

## Verdict: `NO_OPERATOR_PASSES_ADMISSION_GATE`

All three learned operators fail. But the *manner* of failure is informative and
produces one clean positive result, so this is not a null run.

---

## 1. Harness validity — all structural checks pass

These were pre-registered as correctness tests, not criteria.

| Check | Prediction | Measured |
|---|---|---|
| `intercept` on within-document CI, every k | exactly 0 (metric is intercept-invariant) | **+0.0000** |
| `attention` (pure transport) at k=1 | exactly 0 (convex weights = 1 ⇒ constant) | **+0.0000** |
| `ridge`, `krr` at k=1 | exactly 0 (same argument) | **+0.0000** |

The metric is doing exactly what the theory says it does. Note the contrast that
makes this worth stating: `intercept` improves **RMSE by −0.109** at k=5 while
moving within-document ranking by **exactly zero**. Reporting RMSE alone would
have scored calibration as a few-shot adaptation success.

---

## 2. The result table

Within-document CI gain over the support-free base, component bootstrap (n=78).

| arm | k | correct | random target | protein-hard | chemistry-matched | label-permuted |
|---|---:|---:|---:|---:|---:|---:|
| **cnp** | 1 | **+0.0271** | +0.0274 | +0.0270 | +0.0270 | +0.0272 |
| **cnp** | 3 | **+0.0273** | +0.0272 | +0.0268 | +0.0266 | +0.0272 |
| **cnp** | 5 | **+0.0274** | +0.0272 | +0.0270 | +0.0268 | +0.0273 |
| r2d2 | 1 | +0.0166 | +0.0166 | +0.0166 | +0.0166 | +0.0166 |
| r2d2 | 3 | +0.0119 | +0.0048 | +0.0121 | +0.0100 | +0.0088 |
| r2d2 | 5 | +0.0051 | −0.0011 | +0.0060 | +0.0004 | −0.0017 |
| attention | 3 | −0.0011 | −0.0057 | +0.0005 | −0.0000 | −0.0015 |
| attention | 5 | −0.0005 | −0.0064 | −0.0019 | +0.0001 | −0.0014 |
| krr | 5 | +0.0021 | −0.0015 | +0.0008 | −0.0005 | −0.0013 |
| ridge | 5 | +0.0005 | −0.0023 | −0.0000 | −0.0020 | −0.0052 |
| intercept | any | +0.0000 | +0.0000 | +0.0000 | +0.0000 | +0.0000 |

Pre-registered bar (criterion 1): lower bound > +0.005 at k=5. Nothing reaches it
with a support effect.

---

## 3. What actually happened, arm by arm

### `cnp` — the largest number in the study, and it is a shortcut

+0.0274 at k=5 is bigger than the *entire* target-specific adaptation object
(+0.0154, M4) and would have looked like a spectacular few-shot result. It is not.
**Permuting the support labels changes it by 0.0001. Random-target support changes
it by −0.0002.** The operator learned to route around the support set entirely and
compute a function of the query alone.

This is exactly the degenerate solution predicted for a decoder of the form
`g([x_q, r_q])`: nothing forces it to read `r_q`, and a target-agnostic chemistry
head is an easier optimisation target than an adaptation operator. Without the
wrong-support controls this would have been published as a mechanism.

### `r2d2` — collapses to its meta-learned prior, and support makes it worse

At k=1 the ridge correction is identically zero (centred single point), so
+0.0166 is purely the meta-learned prior `w0` — again a global head. As k grows,
the adaptation term **degrades** it: +0.0166 → +0.0119 → +0.0051.

There is a faint genuine support effect (correct +0.0051 vs random −0.0011 at
k=5), but it is not target-specific: **protein-hard support scores +0.0060,
*higher* than correct support.** The nearest non-homologous protein is as good as
the true one, so what little the support contributes is family-level context, not
target identity. The interval on correct (−0.0076 lower bound) does not exclude
zero regardless.

### `attention` — pure transport learns nothing

Exactly 0 at k=1 as predicted, and −0.001 to −0.000 at k=3/5, with correct support
indistinguishable from chemistry-matched wrong support. This independently
reproduces the earlier TRACE null on a different corpus, a different feature basis
and a learned rather than fixed kernel: **nothing label-free predicts which
support pair to trust.**

### `krr`, `ridge` — closed-form baselines behave as M2 predicted

krr at k=5: +0.0021 correct vs −0.0015 random. The correct-minus-random gap
(+0.0036) is the only support effect in the study whose sign is consistent across
all four controls, but its lower bound (−0.0014) does not exclude zero at 78
components.

---

## 4. The clean positive: the base model, not the adaptation, is the binding constraint

The `cnp` arm is a failed operator but a successful measurement. A nonlinear
target-agnostic chemistry head, trained on `meta_train` components and evaluated
on 78 held-out components under full scaffold/document/assay separation, improves
within-document ranking by **+0.0274 [+0.0067, …]** over the linear ridge base.

Set that against everything else measured in this programme:

| quantity | within-document CI |
|---|---:|
| **target-agnostic nonlinear chemistry head** (this run, held-out components) | **+0.0274** |
| full target-specific adaptation object, ~140 labels (M4) | +0.0154 |
| target-specific object at k=20 (M2) | +0.0096 |
| target-specific object at k=5 (M2) | +0.0019 |
| best learned operator's *support-attributable* effect (this run) | ≈ +0.003, n.s. |

**Improving target-agnostic chemistry is worth ~1.8× more than perfectly solving
target adaptation, and ~14× more than what k=5 adaptation delivers.** On this
substrate the few-shot adaptation problem is not where the value is.

RMSE caveat, stated plainly: the `cnp` head is trained on a ranking surrogate and
is badly mis-calibrated in absolute terms (RMSE +0.595 *worse*). The claim is
about ranking only, and a deployable version would need a calibration head.

---

## 5. Admission gate, item by item

| # | Condition | Result |
|---|---|---|
| 1 | improves over support-free model | cnp yes (+0.0274) — but not via support |
| 2 | improves over calibration | cnp yes; others no |
| 3 | improves over ridge/kernel adaptation | cnp yes; others no |
| 4 | **correct support > wrong support** | **ALL ARMS FAIL** — max gap +0.0036, n.s. |
| 5 | gain is query-dependent | cnp yes, but query-only |
| 6 | removing the mechanism removes the gain | **fails** — removing support changes nothing |
| 7 | survives held-out targets | yes (78 held-out components) |
| 8 | effect exceeds uncertainty | only for the support-free part |
| 9 | not explained by ligand similarity alone | untested — no arm got far enough |
| 10 | new inductive bias | no |

**Condition 4 is the one that matters and every arm fails it.**

---

## 6. Why this was predictable from the earlier measurements, and what it costs

M2 measured closed-form ridge against the support-free base at
k=3 **+0.0007**, k=5 **+0.0019**, k=10 **+0.0056**, k=20 **+0.0096**. The
pre-registered bar asked an operator to deliver at k=5 what ridge needs k≈10–20 to
deliver. Nothing did. The label-budget curve did not move.

This is consistent with the identifiability structure already established: the
task family has participation ratio 115 in 266 dimensions, so k ≤ 5 spans a
vanishing fraction of it, and meta-learning the *operator* cannot manufacture
information the support does not contain.

---

## 7. STOP conditions — invoked

Per the pre-registration, this is not to be rescued by more parameters, longer
training, extra losses or additional modules. Three declarations:

1. **The support-conditioned operator track is closed on this substrate.** Three
   operator families (transport/attention, CNP cross-attention, ANIL/R2D2
   closed-form) all fail condition 4.
2. **`cnp`'s +0.0274 must never be reported as a few-shot result.** It is a
   target-agnostic head. Any future use must carry the label-permutation control
   in the same table.
3. **The remaining live question is not adaptation but the base model.** The
   measured ordering says the paper-sized opportunity is target-agnostic
   chemistry under provenance-separated evaluation, where a simple nonlinear head
   already buys +0.0274 over the linear base on held-out components.

## 8. Risks to this conclusion

- **Single seed.** Per programme discipline, multi-seed was reserved for
  mechanisms that pass the gate. None did. A seed sweep would not change
  condition 4, which fails by construction (label permutation is a no-op).
- **Frozen base.** All operators correct a fixed ridge base. An operator trained
  jointly with the encoder might behave differently — but `cnp` and `r2d2` both
  had free encoders and both converged to support-independence, which is evidence
  the failure is not the frozen base.
- **k ≤ 5 only.** The gate targeted the few-shot regime by design. M2 already
  shows the object is recoverable at k ≥ 20 by closed-form ridge with no operator.
- **pIC50-weighted corpus.** 314 of 379 components are pIC50; endpoint-split
  results are in the JSON and show the same pattern.

# A2S Gate T0 — is there a transferable adaptation object?

Date: 2026-08-02 · **Revision 2 (supersedes revision 1 of the same date)**
Branch: `research/a2s-transfer-object-20260802`
Runner: `research/a2s_transfer_object_gate.py` · Tests: `tests/test_a2s_transfer_object_gate.py` (20 passed; `tests/` 302 passed)
Roles opened: source `fit` and `probe` only. `locked` and the recipient roster were not requested.

> ## Status: EXPLORATORY EVIDENCE ONLY
> The source `probe` outcome was consumed once by PIRS and, per the ledger, may
> not be reused for model selection. This gate evaluates on `probe`. Nothing here
> may drive a go/no-go decision; it is diagnostic only. A confirmatory version
> must use nested `fit`-only development tasks.

**FACT — decision: `NO_TRANSFERABLE_CHEMICAL_HEADROOM_OBJECT_IS_MEASUREMENT_CONTEXT`.
All four gates fail.**

| Gate | Question | Revision 1 | **Revision 2 (correct)** |
|---|---|---|---|
| T0A | is the headroom chemical or measurement context? | PASS | **FAIL** (marginal) |
| T0B | does a source head transfer? | PASS | **FAIL** |
| T0C | can a label-free proposal shrink the space? | FAIL | FAIL |
| T0D | can k ≤ 5 select the object? | FAIL | FAIL |

---

# Corrections to revision 1

Revision 1 reported `TRANSFERABLE_AT_FULL_SUPPORT_BUT_NOT_IDENTIFIABLE_AT_K5`
and called discrete transfer "the programme's first positive transfer result".
**That claim is retracted.** An external review identified the defects; I
reproduced every one of them against my own artifacts.

**C1 — the gate was not reproducible, and its verdict flips.** `build_basis` used
`torch.svd_lowrank`, which draws a random projection and accepts no generator.
Two calls in one process returned bases differing by 9.13 in max absolute value.
Replaying the unchanged command returned a different verdict. **This affects every
gate built on that basis** — A0–A4, G1–G4, R0, HOTSPOT — not only T0. Fixed by an
exact symmetric eigendecomposition of the 1024×1024 fit-role Gram matrix with a
sign convention; determinism is now verified across processes by a test.

**C2 — the binding estimand was missing.** T0A itself established that a
chemistry-free document oracle outscores the whole chemical head, which makes any
*all-pair* gain document-confounded by construction. I applied that control to
the own-head arm and **failed to apply it to the transfer arm**, which is where
the headline claim lived. Recorded now for every selected-head and k-shot arm.

**C3 — "trains nothing" was false.** Closed-form ridge heads are fitted from
labels for every source target and every recipient support set. What is absent is
a gradient-trained model. Corrected in the code, the artifact and this report.

**C4 — support size was misstated.** Full support averages **198.1** labels per
target (median 93), not "~64".

**C5 — artifact hashing was broken.** The runner hashed the JSON, appended an
artifact block containing that hash, then rewrote the file — so the recorded
digest never matched the file. Revision 1 quoted `7f54d1c4…` for a file that
hashed to `f05f1fc0…`. The payload is now hashed once, after the artifact block,
and a test verifies the recorded digest against the file on disk.

**C6 — the information account was overstated.** Retracted as a closure; retained
only as an order-of-magnitude heuristic, with its caveats now shipped inside the
artifact. See §5.

---

# 1. T0A — the headroom is measurement context, and the chemical remainder is unresolved

Probe targets, within-target scaffold-disjoint splits, 52 targets / 50 components,
paired component bootstrap (2 000 draws), deterministic basis.

| Arm | CI gain over frozen base | 95 % interval |
|---|---:|---:|
| own fitted head, all query pairs | +0.0519 | [+0.0285, +0.0739] |
| own fitted head, **cross-document pairs** | +0.0814 | [+0.0374, +0.1370] |
| own fitted head, **same-document pairs** | +0.0290 | [**+0.0046**, +0.0557] |
| **document-mean oracle** (labels, no chemistry) | **+0.0610** | [+0.0386, +0.0824] |
| document-mean oracle, same-document pairs | +0.0000 | [+0.0000, +0.0000] |

**FACT.** A chemistry-free predictor knowing only a compound's ChEMBL document and
that document's mean residual scores **+0.0610**, beating the full 26-dimensional
per-target chemical head (+0.0519). **This arm does not depend on the ligand basis
at all and is fully deterministic**, so it is unaffected by C1 and is the single
robust finding of this gate.

**FACT — T0A fails.** The same-document lower bound is **+0.0046**, below the
registered 0.005 threshold. Under the pre-fix randomized basis across six seeds it
ranged **−0.0029 to +0.0052**, clearing the threshold in no seed. The offset-free
chemical headroom is therefore **unresolved**, not established: its point estimate
is consistently positive (+0.024 to +0.030) but it does not clear its own
admission bar.

**INFERENCE.** Q1's long-standing "information only at Tanimoto ≥ 0.55" is
plausibly the same stratum as "same measurement context", since a document reports
a congeneric series. This remains an inference, not a measurement.

# 2. The split is not leakage-free

Measured, not assumed (`provenance_audit`, now in the artifact):

| Statistic | Value |
|---|---:|
| targets sharing ≥ 1 document between support and query | **52 of 52** |
| query rows reusing a support-side document | **91.1 %** |
| query rows reusing a support-side assay | **88.8 %** |
| targets where *every* query row is from a support-seen document | **21** |

Murcko-scaffold disjointness is not document, assay, or congeneric-series
disjointness. Every prior within-target result in this programme inherits this.

# 3. T0B — transfer does not survive the same-document control

| Arm | CI gain over frozen base | 95 % interval |
|---|---:|---:|
| head selected on the target's own support, **all pairs** | +0.0248 | [+0.0031, +0.0476] |
| the same, **same-document pairs** | **−0.0183** | [−0.0435, +0.0054] |
| the same, minus the pooled head (all pairs) | +0.0320 | [+0.0138, +0.0521] |
| best source head selected on eval labels (oracle) | +0.1386 | [+0.1176, +0.1594] |
| **median** source head | −0.0255 | [−0.0402, −0.0118] |

**FACT.** The all-pair gain is positive and the same-document gain is **negative**,
with the interval containing zero. Across six seeds the all-pair gain ran +0.019
to +0.033 while the same-document gain ran **−0.021 to −0.005**, negative in every
seed.

**INFERENCE.** What transfers between targets here is consistent with
measurement-context structure, not target-specific chemical ranking. **The claim
"discrete transfer is real" is withdrawn.**

The oracle arm (+0.1386) exceeds the target's own head, which is a maximum over
110 noisy estimates and is reported only to bound it away from belief.

# 4. T0C / T0D — unchanged, both fail

Label-free shortlists are null (protein +0.0049 [−0.0051, +0.0158]; chemotype
−0.0097 [−0.0214, +0.0011]). k-shot selection at every budget is null or negative
against the base, and **every same-document arm is null**:

| k | all pairs, vs base | same-document, vs base |
|---:|---:|---:|
| 1 | −0.0502 [−0.0938, −0.0191] | +0.0070 [−0.0124, +0.0284] |
| 3 | −0.0183 [−0.0428, +0.0022] | −0.0057 [−0.0216, +0.0124] |
| 5 | −0.0087 [−0.0285, +0.0103] | −0.0044 [−0.0219, +0.0155] |
| 10 | −0.0024 [−0.0244, +0.0171] | −0.0034 [−0.0188, +0.0147] |
| 20 | +0.0112 [−0.0038, +0.0272] | −0.0130 [−0.0329, +0.0066] |

The k = 20 all-pair interval crosses zero, so revision 1's "break-even at k = 20"
was itself not a resolved effect.

**Positive control passes**: in a planted world, k=5 selection recovers the
generating head at 18.2 % versus 0.9 % chance and returns +0.100 CI; k=1 recovery
is exactly chance. The gate has power; the negatives are measurements.

# 5. The information account, retracted as a closure

Revision 1 claimed the bit budget predicted the measured break-even "with no free
parameter". It does not, and the artifact now carries these caveats:

- `tau` is a held-out ridge projection scale, **not** the separation among
  candidate operators — which is what a selector actually decodes;
- `sigma` is total residual dispersion, not conditional observation noise;
- support contrasts are not exchangeable channel uses, because supports share
  chemical design;
- no empirical mutual information or decoder-error bound was estimated.

**The decisive dependency.** `M_useful` was defined post hoc as heads reaching half
the (noisy) oracle maximum, giving 12.8 of 110 and `log2(110/12.8) = 3.10` bits.
Defining it as heads that simply beat the base gives **44.8 of 110** and
**`1.30` bits** — a factor of 2.4 from an arbitrary choice. The apparent agreement
with k ≈ 20 used the same data on both sides. Retained as an order-of-magnitude
statement only: **the support-label channel is small relative to the hypothesis
space**, without a quantitative closure.

# 6. What survives

1. **The document confound** (§1) — deterministic, basis-independent, and the one
   result this gate establishes. A chemistry-free document oracle beats the full
   per-target chemical head.
2. **The provenance audit** (§2) — scaffold-disjoint evaluation on this substrate
   is 91 % document-overlapping.
3. **A new negative** (§3) — selected-head transfer does not survive the
   same-document control.
4. **Mandatory controls** — every future ranking result must report the
   same-document contrast and the document-mean oracle, alongside the
   magnitude-matched and random-selection controls.

# 7. Decision

- **Gate F1 does not run as written.** It was probe-dependent and rested on a
  transfer premise now withdrawn.
- **Nothing is promoted** to `model/` or `script/`. No breakthrough.
- **T0 is retained as exploratory evidence** that scaffold-disjoint evaluation can
  remain strongly document-confounded.
- The successor must use **deterministic, nested `fit`-only development tasks with
  simultaneous target, scaffold, document and assay separation**, and must score
  adaptation on same-document pairs drawn from documents absent from the support.

Artifacts: `reports/active/a2s_transfer_object_gate_2026-08-02.json` — verify with
`content_sha256` over the payload minus that key (a test enforces this);
`..._records_2026-08-02.parquet`; lock `6bcf6edc…`.

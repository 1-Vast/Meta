# Phase 1 — the registered B5 discriminator: result

Date: 2026-08-10.
Preregistration `PREREG_S7_L2B_UNIFIED.md` (sha `2c333f22…`, commit `ce186f4`).
Phase 0 integrity artifacts committed as `139effd` **before** B5 was scored.

```text
ALL SIX REGISTERED GATES PASS
```

## 1. Phase 0 — every blocker discharged

| Item | Result |
|---|---|
| Atom correspondence | **verified at 375,311 positions**; 14,585 admitted, 4 enumerated and quarantined; quarantine wired into the data contract |
| Tie-aware AP + sealed per-pair predictions | **52,062,975** held-out cells sealed per arm as hashed float16; AP reported lexicographic / optimistic / pessimistic + Monte-Carlo expectation |
| Negative sampler | **contract satisfied** — exactly six per positive, unique within each positive block, **zero** negatives that are actually positives |
| Determinism | **verified bit-identical** on identical data in one process |
| Publication/time closure | **built and frozen** from RCSB; 14,426 of 14,447 entries returned, 14,103 with a document key |
| ESM2-650M | **acquired by range-resume and SHA-256 verified** (`c874668852…`, revision `08e4846e`), then run offline |

Two corrections were forced during Phase 0 and are recorded rather than smoothed
over. First, the earlier determinism comparison was **confounded** — the I-1
quarantine changed train from 9,758 to 9,757 between runs, so the differing
checkpoint hash was expected and was *not* evidence of non-determinism; a
same-data test then showed bit-identical state dicts. Second, the ligand-only
baseline turns out to be **massively tie-dependent** (optimistic 0.199 versus
pessimistic 0.003) because a ligand-only model assigns every residue the same
score; only the tie-aware expectation is a defensible point estimate for it.

## 2. Result — held-out A, 2,409 complexes, 196 protein components

Complete residue × heavy-atom matrix, tie-aware macro-AP over components.

| Arm | macro-AP |
|---|---:|
| `B0` prevalence | 0.00319 |
| `BM5` motif shuffle | 0.00451 |
| `BP5` wrong protein | 0.00464 |
| `BL` ligand-only | 0.00572 |
| `BX5` wrong ligand | 0.01968 |
| `B4` non-PLM residue features | 0.02325 |
| **`B5` frozen ESM2-650M** | **0.06960** |

| Gate | Contrast | Δ | LCB95 | |
|---|---|---:|---:|---|
| G1 | B5 − B0 | +0.06642 | +0.05998 | **PASS** |
| G2 | B5 − BL | +0.06388 | +0.05751 | **PASS** |
| G3 | B5 − BP | +0.06496 | +0.05849 | **PASS** |
| G4 | B5 − BM | +0.06509 | +0.05876 | **PASS** |
| G5 | B5 − BX | +0.04992 | +0.04424 | **PASS** |
| G6 | B5 − B4 | +0.04635 | +0.04039 | **PASS** |

Every lower bound clears the frozen 0.02 threshold by at least a factor of two.
Only the residue features changed; the atom branch, head, rank 32, projected
dimension 128, sampler, optimiser, learning rate, weight decay, epochs, seeds,
split, evaluation mask and tie policy are identical to B4.

**This localises the earlier failure definitively: it was biological
representation, not objective or optimisation.** The trainability control had
already excluded optimisation; changing only the residue representation tripled
pair AP and cleared every Gate.

## 3. The gain is entirely residue-side — and largely generic

Marginal decomposition computed from the **sealed** predictions, nothing retrained.

| Arm | residue-marginal AP | atom-marginal AP |
|---|---:|---:|
| `BL` ligand-only | 0.0313 | 0.7246 |
| `B4` | 0.0879 | 0.6895 |
| **`B5`** | **0.2651** | 0.6796 |
| `BX5` wrong ligand | 0.2453 | 0.5097 |
| `BP5` wrong protein | 0.0434 | 0.6595 |

* residue marginal, B5 − B4 = **+0.1772 [LCB +0.1601]**
* atom marginal, B5 − B4 = **−0.0099 [−0.0218, +0.0044]** — no gain, interval spans zero

Two readings follow, and the second is the important caveat.

**(a) ESM2 buys residue localisation and nothing else.** Residue-marginal AP
triples; atom-marginal AP does not move. Atom propensity is already near its
ceiling for a ligand-only model (`BL` = 0.7246, the highest of any arm), which
matches I-2's finding that the true atom marginal is worth only 0.0085 pair AP.

**(b) The residue localisation is largely ligand-independent.** Swapping in a
wrong ligand leaves residue-marginal AP at **0.2453 against B5's 0.2651** — about
**92.5 %** of it survives. B5 is predicting a *generic pocket*, not a
ligand-conditioned one. The pair-level G5 gap (+0.0499) appears in pair scores,
but may still be explained by additive residue and atom marginals rather than
exact coupling; it does not show that residue localisation is ligand-specific.

This is exactly the "strong `pi` alone means generic pocket localisation"
outcome, and it was predicted in advance by I-2, which measured coupling beyond
degree-preserving margins at a median z of only +0.41.

## 4. What is and is not established

**Established.** A frozen-ESM2 residue representation identifies binding
residues far beyond ligand-only, wrong-protein, motif-shuffle and wrong-ligand
controls, under protein-component closure with ligand-graph disjointness, on
sealed and recomputable predictions. Residue localisation is retained as a
**biological statistic candidate**.

**Not established.** Exact residue–atom coupling. The ligand-conditioned part of
residue localisation is small, and coupling beyond marginals was already measured
as weak. No claim of exact coupling is made.

**Not tested.** Affinity direction, transfer, few-shot section identifiability,
`z` admission. The confirmation cohort was never opened.

**Also recorded.** Publication/time closure shows only **2** additional-PDB
entries were released on or after the frozen 2019-01-01 cutoff, and none would
qualify at 2024-01-01. MONN was assembled in 2020 from PDBbind-v2018-era
structures, so it contains **no time-forward holdout by construction**. A
document-closed confirmation cohort remains constructible (707 documents, only 4
shared with development); a time-forward one does not, and would need a different
source carrying residue-atom interaction labels.

## 5. Boundary classification

| Class | Verdict |
|---|---|
| DATA/LABEL INSUFFICIENCY | for marginals: excluded. For **coupling**: confirmed weak (I-2, median z +0.41) |
| BIOLOGICAL REPRESENTATION FAILURE | **resolved for residue localisation** by the frozen PLM |
| OBJECTIVE/OPTIMIZATION FAILURE | excluded — trainability control 0.759, and determinism verified |
| SUPPORT-SECTION NON-IDENTIFIABILITY | not reached; no adaptation attempted |

## 6. Phase 2 precondition

Phase 2 may be entered only after a **new preregistration is committed**. The
evidence above sets its expectation in advance: the residue-first component
should improve, and the coupling component is predicted to be small. A large
reported coupling gain on this corpus should first be suspected of marginal
leakage into the coupling term.

Frozen surfaces unmodified; no affinity, DAVIS or recipient label read; all code
under `research/`.

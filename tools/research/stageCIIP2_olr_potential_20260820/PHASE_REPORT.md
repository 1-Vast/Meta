# CIIP-2 OLR-Potential phase report (Phases 3-5)

Date: 2026-08-20. Base prereg a7b17e8a...; ADD-1 aa8d06af...; ADD-2
91e2cb3a.... No production code touched. All real-data results below are
interpreted under the instrument-qualification ceiling established first.

## 1. Structural qualification (Phase 3)

12/12 structural tests green (tests/test_structure.py): antisymmetry,
cycle-zero, centering identity, coordinate-free deployed signature,
capacity <= 2.0M, permutation/fold wiring, leakage guards (train-parent-only
gain weights; parent-grouped folds), determinism, erased-state equality
49/49, no-closed-form grep, frozen splits (S2/S3/SPB parent-disjoint).

## 2. Instrument qualification (Phase 4 entry; synthetic, frozen protocol)

Plant: parent-DEVIATION interaction field linear in real ESM mean states,
cross-parent variance scaled to the measured between-parent component
(134.8 %^2), plus pair-specific planted component 44.85, noise 44.85, and a
shared ligand pattern 50.0 - matching the audit-measured variance
decomposition of the real panel.

| estimator | test R2 (full centered target) | Delta vs A0 prior |
|---|---:|---:|
| A0 ligand-pattern prior | 0.5361 | - |
| A1 bilinear (raw target) | 0.1046 | -0.4315 |
| A2 router (raw target) | -0.3272 | -0.8633 |
| A5 composed (residual + gains) | 0.5602-0.5657 | +0.026..0.030 |
| A5 rank 2 | 0.5517 | +0.0156 |
| A5 last-epoch checkpoint | -0.3788 | -0.9149 |

A0 prior recovery of the planted panel: sign acc 0.886. Train R2 of A5
0.67-0.68 (the field is learnable on train); test transfer collapses.

Verdict: INSTRUMENT_UNDERPOWERED, robust across rank {2,8}, selection rule
(val early-stop is optimal; last epoch far worse), and estimator family
(bilinear/router/composed). The 0.25-Delta-R2 recovery standard is not met
by any configuration (max +0.03).

Conclusion (frozen consequence): the 49-pair covered panel cannot adjudicate
deployable protein-conditioned interaction learning at the pre-registered
standard. Real-data Phases 4-5 are reported with interpretation ceiling
UNRESOLVED (power) for the deployment claim (R3).

## 3. Phase 4 smoke (real data, single seed 11, split S1)

GATES (prereg section 8, ADD-2 primary metric = full centered target):
- pipeline end-to-end with finite metrics: PASS (PHASE4_SMOKE.json)
- A0-prior sanity in [0.08, 0.18]: PASS (R2 = 0.1313, matches audit value)
- nonconstant coverage 9/9: PASS
- C-perm destroyed: **FAIL** - C-perm composed R2 0.1818 >= A5 0.1388

Key numbers (test = 9 covered pairs / 6 parents):

| arm | R2 | Spearman | sign | var-rec | note |
|---|---:|---:|---:|---:|---|
| A0-prior | +0.1313 | 0.329 | 0.697 | 0.100 | analytic baseline (matches audit) |
| A1-bilinear | +0.1087 | 0.267 | 0.670 | 0.127 | deployable, no router |
| A2-router | +0.0252 | 0.133 | 0.581 | 0.135 | deployable router, raw target |
| A3-oid | +0.0252 | 0.133 | 0.581 | 0.135 | IDENTICAL to A2 (identity theorem, AM-2) |
| A4-cfoie | +0.1732 | 0.359 | 0.746 | 0.229 | composed; best arm |
| A5-gain | +0.1388 | 0.351 | 0.746 | 0.291 | composed + gain weights |
| C-perm | +0.1818 | 0.380 | 0.744 | 0.108 | permutation control ranks FIRST |

- Model-contrast increment over the prior: A4 +0.042, A5 +0.0075 R2 -
  single-seed, both far below the instrument recovery standard (0.25) and
  below the permutation control. A3 == A2 exactly, confirming the
  orthogonal-decomposition identity theorem in live training.

**Phase-4 verdict: FAIL (gate b). Phase 5 is NOT authorized.**

Interpretation under the instrument ceiling: consistent with
INSTRUMENT_UNDERPOWERED. The deployable increment that would need to be
>= 0.25 R2 (instrument standard) is +0.0075 R2 (real, single seed), and the
permutation control cannot even rank the correct arm above chance.

## 4. Phase 5 SPB parent-disjoint evaluation (real data)

NOT EXECUTED. The frozen chain (prereg section 8) authorizes Phase 5 only
after Phase-4 PASS; Phase 4 failed gate (b). Independently, the instrument
qualification already established that the 49-pair panel cannot recover a
planted transferable interaction field at the pre-registered standard, so
a multi-seed SPB run could only add noise-level measurements. The SPB
split construction and code remain frozen and tested for any successor
programme with adequate data.

## 5. Verdict ladder (programme level)

- R1 representation (mutation sensitivity of deployable residue features):
  SUPPORTED (CIIP-1A context audit: site delta 4.01 vs context 0.057;
  49/49 pairs; X-erasure exact; random window not mutation-null).
- R2 identification (pair-level; mutation information changes
  ligand-conditioned prediction beyond controls): NOT SUPPORTED at this
  power. No deployable arm exceeds the shared-ligand-pattern baseline
  (CIIP-1A correct 0.0075 vs baseline 0.1313; OLR-Potential A5 increment
  +0.0075 with permutation arm higher still).
- R3 deployment (parent-disjoint transfer above ligand prior): UNRESOLVED
  (power). The synthetic instrument planted a transferable,
  sequence-linear parent field of realistic magnitude and no tested
  estimator recovered it beyond +0.03 R2 (standard 0.25); real-data Phase 4
  shows the corresponding increment is +0.0075.
- R4 bridge (zero/few-shot DTA gains): BLOCKED (requires R3 SUPPORTED).
- R5 mechanism/binding-affinity interpretation: NOT CLAIMED (functional
  single-dose endpoint; endpoint separation enforced throughout).

## 6. Stop-rule status

User termination rule 2 ("ligand permutation does not destroy results") is
FORMALLY TRIGGERED at the smoke level: there is no interaction result above
the shared-pattern prior for permutation to destroy. This is not an assay
artifact alarm (the A0 sanity and nonconstant gates pass; the prior itself
behaves exactly as audited); it is the honest absence of detectable
protein-conditioned interaction signal at this panel scale, with the
instrument showing the absence cannot be upgraded to a strong negative
claim either.

## 7. Successor requirements (what would unblock R3)

1. Same-endpoint panels pooled to >= 100 independent mutation conditions
   across >= 30 parents (single-dose % inhibition only, never mixed), or a
   Ki/Kd DeltaDeltaG endpoint corpus (Platinum-scale) with its own prereg.
2. Replicate/noise characterization to separate the mutation-specific
   variance component (currently 40% incl. noise, no replicates).
3. The OLR-Potential estimator, SPB split, controls, and instrument are
   frozen and reusable as-is; first successor action is re-qualification
   on the larger panel before any biological claim.

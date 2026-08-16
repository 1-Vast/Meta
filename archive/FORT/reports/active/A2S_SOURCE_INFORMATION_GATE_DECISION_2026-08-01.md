# A2S Source Information-Gate Decision

Date: 2026-08-01  
Status: source-only diagnostic; no recipient or locked-role labels opened  
Model family: `a2s_information_gate_v1`

## Decision

**STOP before implementing another adaptation model.** The new source-only
probe does not establish the two necessary quantities
`Delta_label > 0` and `Delta_assign > 0` on the new probe role. A synthetic
positive control is recovered at all three budgets, so this is not an obvious
implementation-power failure. The locked role remains untouched and must not
be opened for model selection after this failed probe gate.

This is a bounded information-gate result, not an information-theoretic proof
that no meta-adaptation mechanism can exist.

## FACT: Frozen Source Lock

The metadata-only lock is
`reports/active/a2s_source_information_gate_lock_2026-08-01.json` with content
SHA-256 `ee73601e061d735fcaf9ae2cafc7f40c25300086de666b1174e75297c5457307`.

- The registry was restricted to `dual_cold_split=train`, `endpoint=pKi`.
- The lock requested no affinity column and used target homology plus every
  pipe-delimited document/assay token as a provenance closure.
- There are 559 source targets, 141 closed components, and 185,591 source
  metadata rows.
- Roles are 484 fit targets / 84 components, 41 probe targets / 28 components,
  and 34 locked targets / 29 components.
- Target, homology, document-token, and assay-token overlaps between every pair
  of roles are all zero.
- The largest provenance component contains 380 targets (68.0% of all source
  targets). This is a material power limitation, not a model result.
- A label-free topology check shows that homology-only closure has 517
  components with maximum size 4, but adding exact document cells reduces this
  to 163 components with a 347-target (62.1%) giant component. Assay cells do
  not connect targets in this registry. The imbalance is therefore caused by
  the available document provenance, not by the tokenization implementation
  alone.
- Metadata-only feasibility counts are 37/33/30 probe targets with at least
  6/8/10 rows for the nominal k=1/3/5 plus query requirement. The actually
  constructible probe episodes are smaller because support/query closure is
  stricter.

## FACT: Label and OOF Firewall

The executable gate is `research/a2s_information_gate.py`.

- Only labels for the new `fit` and `probe` roles were requested through a
  target-filtered Parquet read.
- The `locked` target list was not requested; recipient roles were not present
  in the source registry read.
- The base is a target-balanced ligand design plus a fixed 16-dimensional
  protein projection ridge. Every fit component receives a prediction from a
  model trained without that component. Probe rows are predicted by a model
  trained only on fit components.
- The OOF artifact is
  `reports/active/a2s_source_information_gate_oof_2026-08-01.npz`.
  It contains 183,379 finite source predictions (181,575 fit and 1,804 probe)
  and no object arrays or pickled metadata.
- Five component folds were used for fit-role OOF. No in-sample fit residuals
  were used by the probe.
- The OOF folds are severely unbalanced because the 380-target component is
  intact: fold 0 holds 176,193 rows and trains on only 5,382 rows; the largest
  fold is 97.0% of all fit-role held-out rows. This satisfies the exclusion
  contract but makes the real G0/G1 result **non-confirmatory**.

## FACT: Episode and Probe Counts

The completed source diagnostic contains 858 label-bearing episodes:

| role | k=1 | k=3 | k=5 |
|---|---:|---:|---:|
| fit | 327 | 259 | 228 |
| probe | 20 | 14 | 10 |

The probe evaluation has 10, 7, and 4 independent provenance components at
k=1, 3, and 5 respectively. These counts are the statistical units for the
bootstrap intervals below; rows and query compounds are not independent units.

## FACT: Information-Gate Results

The loss is `1 - NDCG@10`; positive `Delta_label` or `Delta_assign` means the
label-aware/correct-assignment arm has lower loss. Intervals are component
bootstraps over the probe role.

| k | Delta_label mean [95% CI] | Delta_assign mean [95% CI] | interpretation |
|---|---:|---:|---|
| 1 | -0.00123 [-0.00277, 0.00035] | undefined | no within-support assignment permutation exists |
| 3 | -0.00358 [-0.01500, 0.00167] | -0.00443 [-0.01031, 0.00076] | neither label information nor assignment specificity is admitted |
| 5 | -0.00076 [-0.00798, 0.00216] | +0.00132 [-0.00614, 0.00684] | point estimates are effectively null and underpowered |

Secondary RMSE deltas do not rescue the gate. At k=3, for example,
`Delta_label_RMSE = +0.0961` (worse for G1), while the ranking interval still
crosses zero. Ranking, rather than calibration, is the primary estimand.

## FACT: Headroom and Positive Control

The high-data target-specific oracle has positive point headroom at k=1/3/5,
but rank-loss intervals cross zero because the probe role has only 10/14/10
oracle targets. The RMSE headroom interval is positive at k=3 and k=5. Thus
the frozen base is not shown to be perfect, but the present evidence is not a
powerful headroom certificate either.

The synthetic support-code control uses the same episode geometry and a
predeclared query-dependent label channel. It is recovered at every budget:

| k | synthetic Delta_label (rank-loss) |
|---|---:|
| 1 | +0.4181 |
| 3 | +0.4219 |
| 5 | +0.2619 |

This control demonstrates that the diagnostic pipeline can detect a known
support-label signal. It is not evidence that real ChEMBL support labels carry
the required signal.

## INFERENCE: Failure Diagnosis

1. **Incremental support-label information is not admitted by this probe.**
   The correct-label channel does not improve query ranking beyond the same
   support-free design and OOF base, and permuting residual assignments does
   not produce a stable penalty.
2. **The OOF residual geometry is severely unbalanced.** The giant provenance
   component occupies 97% of all fit-role OOF holdout rows in one fold. This is
   leakage-safe but not an adequately balanced source-learning regime, so the
   null cannot be treated as a clean confirmatory estimate.
3. **The current episode construction may be too chemically/global for k<=5.**
   The probe uses strict support/query closure and the available source role
   has only 4 independent components at k=5. The null can therefore reflect
   weak identifiability, distribution shift, or limited power rather than a
   universal impossibility.
4. **The provenance closure is conservative but costly.** A 380-target giant
   component is necessary for the current no-overlap contract, yet it leaves a
   small independent probe. A different, biologically justified provenance
   source would be required before interpreting a null as definitive.
5. **The high-data oracle does not imply k-shot recoverability.** It only shows
   that many target labels can sometimes improve a target-specific fit; it does
   not show that one to five labels identify a transferable query correction.

## HYPOTHESIS: What Remains Open

- A same-assay/MMP-connected source stratum could contain assignment-specific
  information that is diluted by the global closure. This must be predeclared
  and tested without opening the locked role.
- A richer full encoder/head could change the bounded-probe result, but such an
  expansion is not justified until a valid information stratum or a new
  provenance-rich dataset is obtained.
- The current result cannot select CMAL, CSRIO, STOP/SWAP, SAR grammar, or any
  other final mechanism. It only blocks another high-capacity attempt under the
  present global episode distribution.

## Required Next Action

Do not read locked labels and do not implement a discrete rank policy. Preserve
the lock and OOF artifact. If research continues, use the English handoff in
`A2S_SOURCE_GATE_FAILURE_HANDOFF_2026-08-01.md` to obtain an independent review
of identifiability, provenance closure, and the choice of source strata. Reopen
only after a new pre-registered information test has either (a) a biologically
defensible, adequately powered support/query stratum or (b) a new dataset with
assay and campaign provenance.

## Verification

- `D:\anaconda\envs\drug\python.exe -m py_compile research\a2s_source_lock.py research\a2s_information_gate.py` passed.
- The CUDA gate completed successfully on the RTX 4060 Laptop GPU.
- No formal recipient training, recipient evaluation, model-folder promotion,
  commit, or GitHub push was performed.

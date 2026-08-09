# Phase 2A preregistration — computational amendment 01

Parent: `research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2A.md`
(SHA-256 `4e01401d0468e3241bd05cde93b2a517919209d770f8f75ab471e42228f1b12e`).

Written 2026-08-10, **after** Phase 0 (contract audit) and Phase 1 (census)
completed, and **before** any Phase 2 teacher-conditionality metric was
computed. The parent file is not edited, so its hash remains valid.

None of the items below is result-dependent. They are numerical and budget
declarations required to execute sections 6–8 of the parent, and each is fixed
here before the corresponding quantity is observed.

## A1 — expected AP under ties is computed in closed form

The parent registers tie-aware AP. For a source ranking with large tied blocks
(a binary residue mask ties every residue into two blocks), the Monte-Carlo
expectation used in Phase 0/1 is noisy and expensive. Phase 2A therefore uses
the exact expectation of average precision under a uniformly random ordering
within each tied block.

For a tied block of `n` items containing `k` positives, preceded by `a` items
containing `b` positives, the expected sum of precisions over the block's
positives is

```text
E = k * (1/n) * sum_{j=1..n} ( b + 1 + (k-1)(j-1)/(n-1) ) / (a + j)
```

and `E[AP]` is the sum over blocks divided by the total positive count. This is
exact, not an approximation. It is validated in-run against the Monte-Carlo
estimator; the maximum absolute discrepancy is reported and must be consistent
with Monte-Carlo error.

## A2 — minimum pair count for a within-construct rank correlation

`T5`/`T6` compute a Spearman correlation across the pairs of one construct. A
rank correlation is not numerically defined below three pairs and is
uninformative just above it. Constructs contributing fewer than **5**
within-construct scaffold-distinct pairs are therefore excluded from `T5`/`T6`
only, their count is reported, and they still contribute to `T1`–`T4` and `T7`.

This is a numerical-validity condition on the estimator, not a selection on
observed values.

## A3 — ligand fingerprint construction

Tanimoto distance uses RDKit Morgan fingerprints, radius 2, 2048 bits (already
frozen in the parent, section 6/T5), computed on the sanitized heavy-atom RDKit
molecule — the same molecule object from which `graph_key` and the Murcko
scaffold were derived, so chemical distance and ligand identity refer to one
object.

## A4 — enumeration is exhaustive

All within-construct pairs are enumerated exhaustively. No subsampling, no cap,
no per-construct truncation. Constructs with many ligands therefore contribute
many pairs, which is precisely why the inference unit is the closure component
and why every headline number is a component-macro.

## A5 — the foreign-ligand corruption control does not exist at teacher level

Residue indices are comparable only within one `seq_key`. A foreign ligand's
residue mask belongs to a different construct and cannot be transferred. The
arbitrary-foreign-ligand control that exists for the model arms (`BX5`) has no
teacher-level counterpart. The corruption floor reported for the teacher is
therefore the prevalence/random-mask baseline, and this is stated in the output
rather than substituted silently.

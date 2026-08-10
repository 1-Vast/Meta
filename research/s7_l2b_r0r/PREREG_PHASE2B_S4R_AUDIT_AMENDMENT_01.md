# Amendment 01 — P1R2B-PHASE2B-S4R-A

Parent: `PREREG_PHASE2B_S4R_LIGAND_REPRESENTATION_AUDIT.md`, SHA-256
`8c3be16973957c6d1e7e735a7c3214d1d7e4f3b5d59f791e3f72271894130138`,
committed `8a643e8`.

Written 2026-08-10, before any audit code was written and before any audit
statistic, A-gate outcome or S4R metric existed.

## 1. Defect in the parent

Parent section 4 specifies

```text
GetMorganGenerator(fpSize = 2**20) ... GetCountFingerprint
```

and asserts that "no folding collision occurs". That assertion is not
guaranteed. `GetCountFingerprint` folds the Weifeiler-Lehman environment hash
into `fpSize` buckets, so with 41,244 distinct radius-2 environments the
expected number of colliding pairs in a `2**20` space is of order `8 x 10^2`,
not zero. A collided vocabulary entry would merge two chemically different
environments into one coordinate and would silently understate collapse.

## 2. Correction

The audit uses the **unfolded** sparse identifiers

```text
GetMorganGenerator(radius = r).GetSparseCountFingerprint(mol).GetNonzeroElements()
```

whose keys live in the full 32-bit environment-hash space. `fpSize` is not used
and no folding is performed. Every other element of parent section 4 —
`radius r in {1, 2}`, `d in {128, 256, 512}`, per-heavy-atom normalization,
train-only top-`d` document-frequency vocabulary with ascending-identifier
tie-break — is unchanged.

## 3. Additional reported evidence

The audit must additionally report, per radius and non-gating:

```text
n_environments_unfolded      distinct identifiers in the 32-bit space
n_environments_folded_2p20   distinct identifiers after folding to 2**20
folding_collisions           the difference between the two
```

This preserves the parent's intent — a vocabulary of exact chemical
environments — and replaces its unverified assertion with a measurement.

## 4. Unchanged

Sections 1, 2, 3, 5, 6, 7, 8, 9, 10 and 11 of the parent stand unmodified. No
threshold, statistic, selection rule, firewall or terminal verdict is altered.
No label, affinity value or S3R metric is opened by this amendment.

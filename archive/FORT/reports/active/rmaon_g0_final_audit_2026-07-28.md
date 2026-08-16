# R-MAON G0 final adversarial audit

**Date:** 2026-07-28

## Accepted evidence

The accepted artifact remains `reports/active/rmaon_g0.json`, SHA-256
`3276a854339cf9ebb7684df2295bb526681004f78847a8dda34dc6a03dacd4cc`.
All 1,000 replicate records were independently recomputed from the JSON. Counts, rates, recovery
medians, and the exact binomial intervals in `history.md` agree. The generated decision is a
byte-identical reconstruction of the accepted JSON.

Dynamic tracing found two affinity-bearing Parquet calls. Both used a reader-level
`dual_cold_split == train` filter, and the coefficient loader still completed when
`DualCold.panel()` was replaced by an exception. No development, Davis, confirmation, or sealed
affinity was read.

## Findings and disposition

1. The TRAIN-only audit fields were self-reported and the regression suite did not dynamically verify
   the Parquet filter. The post-run safety hardening centralizes both reads in
   `read_train_registry`, validates every returned split row, and adds a monkeypatched regression
   test that records the actual `columns` and `filters` arguments and rejects a non-TRAIN return.
   This does not change the accepted numerical artifact.
2. The accepted JSON binds the run-time `task.md` SHA-256
   `d9e20347fdc141adea05f1e5542f642d8ece72911b301d2f421f132595307677`.
   The first result-updated ledger observed before this audit note and safety hardening had SHA-256
   `5c3e40224102be08d82a07e4383af3c91cc8160ad999f1b2948136c0a3d27aa8`.
   The final post-audit `task.md` has SHA-256
   `38add393ca875f9a1c04dd9bef4539cd1b2f8dd3d4057bf0aaa00359d3bc0f99`.
   No byte-exact snapshot of the earlier task ledger was retained, so that one provenance input
   cannot be independently reconstructed from the current workspace. The preregistration, R1
   coordinate, and R1 report hashes do match their current files. This limits provenance
   completeness but does not alter the deterministic replicate records or their interpretation.

## Final boundary

I2 is supported only as a synthetic regular-null estimator module. I1 and real strict dual-cold
prediction remain untested. The final category remains
`SIGNAL_PRESENT_EVIDENCE_INSUFFICIENT`, and G0-A continues to block real-label training.

## Post-hardening verification

Using `D:\anaconda\envs\drug\python.exe`, the R-MAON/K-LBP suite returned `77 passed` and the full
repository suite returned `352 passed`. A real-data reader smoke loaded 12,574 affinity rows and
12,574 accession rows; both scans contained only the `train` split. The complete CPU coefficient
loader also finished with 111 eligible targets and 100 homology components. `git diff --check`
reported no whitespace error (only pre-existing line-ending warnings).

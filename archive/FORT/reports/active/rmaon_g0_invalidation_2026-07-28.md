# R-MAON G0 pre-A1 artifact invalidation

Two concurrent pre-A1 formal writes are invalid and must not be cited.

1. The first observed artifact had SHA-256
   `abd60b7291b5411b56f917c4c46ca2cb96d7e5509724f52f70f018e1f8d9efd0`. Its loader called
   `DualCold.panel()`, which materialized development affinity before using only component metadata.
   The artifact incorrectly declared `panel_development_labels_read=false`.
2. Before it could be archived, a concurrent agent overwrote that path with a repaired-loader attempt.
   The replacement was generated before Amendment A1 explicitly froze the repaired protocol and is
   therefore also not accepted. It is preserved as
   `reports/active/rmaon_g0_invalidated_pre_a1.json`, SHA-256
   `254536f40e55040be5fe3e32fe49e1694a8145f0185fc083579289476f25c932`.

Neither write is interpreted. The scientific thresholds, seed, regimes, and replicate count were not
changed in response. The first accepted result must be generated after A1 using the reader-filtered
TRAIN loader, exact PSD covariance factors, passing amended tests, and bound input hashes.

`sealed_test_consumed=false`; no Davis, confirmation, or sealed label was read.

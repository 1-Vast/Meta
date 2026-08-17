# Final boundary audit — verification report (no training)

Date: 2026-08-18 (night). This audit re-derives the load-bearing numbers of
report/BOUNDARY_20260817_NIGHT.md from the recorded artifacts, with no
training and no meta_test access. Authority: FINAL_BOUNDARY_AUDIT.json.

## Verified numbers

| claim | recomputed | stored authority | match |
|---|---|---|---|
| T2 k=0 MSE / level^2 / centered | 2.5961 / 1.7314 / 0.8648 | boundary doc | exact |
| within-document level transfer R^2 | +0.4515 (210 targets) | D0b_DOC_TRANSFER.json | exact |
| K2 pooled k=0 MSE contrast | -0.1118 [-0.1851, -0.0490] | K2_multiseed_contrast.json | bitwise |
| K2 pooled k=1..5 MSE contrasts | -0.0480 / -0.0273 / -0.0218 / -0.0122, all hi < 0 | same | bitwise |
| meta_test seal | 104 RESULT.json artifacts; 0 evaluated; only 2 legacy R14 artifacts record included=True (pre-cleanup, disclosed) | seal_record in every artifact | consistent |
| preregistration ordering | all 7 training stages have PREREGISTRATION.md alongside their results | per-stage dirs | consistent |

## What this establishes

The final bounded conclusion is reproducible from the raw evaluation rows:
the level/shape decomposition, the assay-history transfer measurement, and
the three-seed pooled contrast of the strongest mechanism (K-REG) all
re-derive exactly. No current-stage artifact opens the meta_test seal
(logical exclusion after parsing; 768 cells withheld), and every trained
stage was preregistered before its results existed.

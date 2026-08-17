# Stage I report — live ESM-2 150M LoRA lane: REJECTED (ranking-only gain)

Development evidence, single seed, meta_val read once after freezing;
meta_test sealed. Authorities: I_meta_val.rows.summary.json,
IFROZEN_meta_val.rows.summary.json, I_vs_IFROZEN.contrast.json,
I_vs_T2.contrast.json, PREREGISTRATION.md, per-arm RESULT.json.

## Verdict

**Rejected by the preregistered gates; nothing promoted.** G2 fails: no MSE
gain at k in {2,3,5} resolves against the frozen T2 baseline (k=2 -0.0459
[-0.123, +0.024], k=3 -0.0383 [-0.093, +0.013], k=5 -0.0182 [-0.054,
+0.017]). Stop rule S1 fires. G1/G3/G4 pass: k=0 is not degraded by a
resolved interval, no ranking metric is degraded by a resolved interval, and
all support/protein controls stay above correct.

Two resolved POSITIVE findings, neither of which the gates require:
- LoRA training vs the matched frozen live path (I - I-FROZEN): k=2 Spearman
  +0.0744 [+0.0118, +0.1420] resolved;
- I vs frozen T2: k=3 Pearson +0.0553 [+0.0007, +0.1205] resolved.

LM conditioning therefore marginally sharpens within-target ranking without
moving MSE, and the trained adapters slightly WORSEN level versus the frozen
live path (k=0 level^2 +0.113, unresolved). The lane is closed for the
level bottleneck; the k=2/k=3 ranking hints are recorded as observations, not
promotion evidence.

## Numbers (frozen meta_val banks, component-weighted, restored pK^2)

| arm | k=0 MSE | k=1 MSE | k=2 MSE | k=3 MSE | k=5 MSE | k=5 Spearman | k=5 CI |
|---|---|---|---|---|---|---|---|
| T2 (frozen bank) | 2.5961 | 1.7712 | 1.3245 | 1.2197 | 0.9859 | 0.3141 | 0.6188 |
| I-FROZEN (live, frozen) | 2.3371 | 1.6595 | 1.2457 | 1.1543 | 0.9441 | 0.2902 | 0.6042 |
| I (live, LoRA trained) | 2.4744 | 1.7018 | 1.2786 | 1.1814 | 0.9677 | 0.3320 | 0.6267 |

Cost: I trainable = T2 1.80M + LoRA 1.23M (1.68x, within the 2.0x budget);
wall 931 s; peak VRAM ~5.6-6.6 GB during training (bounded by the
first-chunk gradient rule); GPU verification recorded in every arm.

## Engineering findings (recorded for any future LM lane)

- Long proteins x chunked LoRA backward exhausted 8 GB VRAM and hard-killed
  the process (silent CUDA OOM); bounding gradients to the first 1022-residue
  chunk fixed it. Any future LM-finetuning lane must carry this bound or
  gradient checkpointing.
- The live frozen encoder beats the cached fp16 bank in this single seed
  (k=0 -0.259, resolved); whether that is precision, checkpoint selection or
  seed noise is unresolved and left as a disclosed open question.

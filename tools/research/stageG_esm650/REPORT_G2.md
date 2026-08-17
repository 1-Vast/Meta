# Stage G2 report — multi-seed confirmation of the ESM-650M lane: NOT CONFIRMED

Development evidence; meta_val read once per seed after each checkpoint was
frozen; meta_test sealed and never opened. Authorities:
G2_multiseed_contrast.json, per-seed row summaries (G_s15/16/17, T2_s15/16/17),
PREREGISTRATION_G2.md.

## Verdict

**The lane is not confirmed; no meta_test is opened.** Per the preregistered
gates: G2-1 fails (k=0 centered MSE is not lower in all three seeds:
G 0.8252/0.8473/0.8583 vs T2 0.8648/0.8636/0.8130; pooled -0.0035
[-0.0244, +0.0168]); G2-2 fails (no pooled resolved MSE improvement at any
k). The single-seed Stage G pattern was a favourable-seed artifact — exactly
the failure the multi-seed protocol exists to catch.

## Per-seed numbers (frozen meta_val banks, component-weighted)

| arm-seed | k=0 MSE | k=0 centered | k=0 Spearman | k=5 MSE | k=5 Spearman |
|---|---|---|---|---|---|
| T2 s15 | 2.5961 | 0.8648 | 0.0790 | 0.9859 | 0.3141 |
| T2 s16 | 2.9811 | 0.8636 | 0.0594 | 1.0072 | 0.3098 |
| T2 s17 | 2.4581 | 0.8130 | 0.0732 | 0.9458 | 0.2855 |
| G s15 | 2.3826 | 0.8252 | 0.1305 | 0.9442 | 0.3361 |
| G s16 | 2.2391 | 0.8473 | 0.0371 | 0.9464 | 0.2870 |
| G s17 | 2.7904 | 0.8583 | 0.0373 | 0.9865 | 0.2972 |

Pooled 3-seed G-minus-T2 (component bootstrap): k=0 MSE -0.2078
[-0.7001, +0.2128]; k=5 MSE -0.0206 [-0.0586, +0.0134]; k=5 Spearman
+0.0036; k=0 CI -0.0030. Every interval crosses zero.

## What is and is not concluded

- ESM-650M residue inputs do NOT produce a resolved improvement over the
  governed 150M bank in any k/metric across three seeds. The lane is closed
  as an external-representation improvement.
- The retraining spread remains the dominant effect: k=0 MSE spans
  2.46-2.98 for T2 and 2.24-2.79 for G across seeds.
- The k=5 band sits at 0.944-1.007 (T2/G, three seeds): at or just below
  1.00 with honest controls, unchanged by the lane.

## Follow-up

No tuning of this lane is authorized. The falsification ledger now includes
both ESM-150M and ESM-650M frozen embeddings as inputs; the remaining
unrun lanes (MSA: no governed UniRef snapshot; structure/pocket: 15/499
exact holo coverage locally; ESM fine-tuning; contrastive coembedding) are
recorded as open questions in the boundary document, not as active stages.

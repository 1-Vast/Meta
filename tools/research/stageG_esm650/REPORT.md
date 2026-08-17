# Stage G report — ESM-650M residue-input trunk: single-seed screen, promising lane

Development evidence, single seed, meta_val read once after freezing;
meta_test sealed. Authorities: G_meta_val.rows.summary.json,
G_vs_T2.contrast.json, PREREGISTRATION.md, RESULT.json under
report/meta_fewshot/stageG_esm650_20260817/G/.

## Result

The ESM-650M residue-input trunk improves the frozen T2 baseline on MSE,
level^2, centered MSE, Spearman, Pearson, CI and cliff sign at EVERY k — the
first arm in this project's record to move MSE and ranking in the same
direction across all five support sizes. One interval resolves: k=0 centered
MSE -0.0396 [-0.0772, -0.0018]. The k=0 MSE gain (-0.2136) and the k=5 MSE
gain (-0.0417) do not resolve; Spearman gains do not resolve.

## Numbers (frozen meta_val banks, component-weighted, restored pK^2)

| k | T2 MSE | G MSE | G-T2 MSE | G-T2 centered | G-T2 Spearman | G-T2 CI |
|---|---|---|---|---|---|---|
| 0 | 2.5961 | 2.3826 | -0.2136 [-0.77, +0.17] | -0.0396 [-0.077, -0.002] R | +0.0515 | +0.0138 |
| 1 | 1.7712 | 1.6878 | -0.0834 | -0.0396 R | +0.0515 | +0.0138 |
| 2 | 1.3245 | 1.2723 | -0.0522 | -0.0298 | +0.0393 | +0.0080 |
| 3 | 1.2197 | 1.1782 | -0.0416 | -0.0287 | +0.0200 | -0.0001 |
| 5 | 0.9859 | 0.9442 | -0.0417 | -0.0265 | +0.0220 | +0.0071 |

Cliff sign: G 0.6191/0.6383/0.6407/0.6874 at k=0/2/3/5 vs T2
0.5299/0.5710/0.5677/0.6086.

Controls (no inversions): matched-wrong 3.59/4.52/5.06/5.82 at k=1/2/3/5;
permuted 1.79/1.57/1.37 at k=2/3/5; wrong-protein above correct at every k.

Cost: 2,044,593 trainable (1.137x T2 — the G5 parameter criterion as written
was not met; wall 172.6 s, peak VRAM 516.7 MB, both below T2's).

## Gates

G1 PASS (k=0 improved), G2 FAIL (no resolved MSE gain in k in {2,3,5}), G3
PASS (no resolved ranking degradation), G4 PASS, G5 partially failed on the
parameter count criterion only.

Per the preregistered stop rules the lane stops at the single-seed screen.
Because this is a new INPUT lane with a resolved centered gain and a
consistent all-k pattern (not an inert mechanism), the continuation is a NEW
preregistered stage — the multi-seed confirmation — rather than a tuning step:
see Stage G2 PREREGISTRATION.

# QPSMP-BPSF v2 research screen

## Decision

Neither candidate is promoted to `model/` or `scripts/`. The relevance-weighted pair readout improved
the endpoint point estimate in one short development seed, but failed the correct-protein control.
The shared-latent candidate did not improve the corrected baseline. In every arm, learned SAR was
effectively zero after 60 optimizer steps.

These runs are consumed-development screens, not inferential Gates.

## Protocol repair completed before the screen

- residue, atom, and adjacency padding is masked before and after local propagation;
- support and query ligand identities are fail-closed disjoint;
- validation uses nested support prefixes and a common query set across `k=1,2,3,5`;
- training caches are refreshed after one pass instead of being replayed indefinitely;
- training diversity is recorded;
- nested evaluation now records query-level CI coverage, RMSE, and Spearman in addition to MSE.

All arms used seed `20260871`, 60 GPU steps, 2 episodes per step, dynamic caches of 8, mixed
`k=1,2,3,5`, and the same consumed development sampling contract. The test subset was limited to at
most two targets per component. Results are therefore diagnostic only.

## Results

### Protocol-faithful two-stage screen

After the initial diagnostic, the runner was corrected to implement the requested separation:
Stage A used only zero-support endpoint MSE; Stage B inherited and froze that checkpoint and used
only centered query error. The second screen used seed `20260872`, 30+30 GPU steps, the same dynamic
episode contract, and two targets per component.
Each stage exposed 60 episodes spanning 53 distinct targets. Peak allocated CUDA memory was about
2.04 GB (baseline), 2.16 GB (relevance), and 1.94 GB (shared latent).

| variant | Stage-A val MSE | Stage-B centered val MSE | k=5 full MSE | k=5 level MSE | k=5 SAR gain | k=5 wrong-state gap | wrong-protein zero-shot gap |
|---|---:|---:|---:|---:|---:|---:|---:|
| corrected baseline | 2.4482 | 0.665872 | 1.1717 | 1.1173 | -0.0000006 | -0.00000006 | +0.3128 |
| relevance weighted | 2.4488 | 0.665755 | 1.1620 | 1.1173 | +0.0000003 | +0.00000002 | +0.8426 |
| shared latent | 2.4481 | 0.665869 | 1.1717 | 1.1173 | -0.0000001 | -0.00000001 | +0.3129 |

The relevance arm improved k=5 endpoint MSE by only `0.0097 pK^2` against the research baseline,
while still losing to level-only by `0.0446 pK^2`. Its positive wrong-protein zero-shot diagnostic
shows that the endpoint does depend on the correct protein, but the support-derived state remains
indistinguishable from wrong state and contributes only about `3e-7 pK^2`. This does not satisfy the
requirement that the meta-section be a material performance source.

### Initial joint-loss diagnostic (superseded for decision making)

| variant | best validation MSE | k=5 full MSE | k=5 level MSE | k=5 SAR gain | k=5 wrong-state gap | wrong-protein zero-shot gap |
|---|---:|---:|---:|---:|---:|---:|
| corrected baseline | 1.3562 | 1.4200 | 1.3338 | 0.0000001 | -0.00000001 | +0.3448 |
| relevance weighted | 1.3238 | 1.3489 | 1.3338 | -0.0000149 | +0.0000021 | -0.0967 |
| shared latent | 1.3589 | 1.4221 | 1.3338 | +0.0000004 | +0.00000004 | +0.2946 |

Positive `SAR gain` means `MSE(SAR-cut)-MSE(full)`. Positive wrong-protein gaps mean the correct
protein or state is better.

## Failure diagnosis

1. The relevance-weighted pair readout produced a small endpoint improvement and positive
   wrong-protein zero-shot diagnostic in the corrected two-stage run, but it remained worse than
   level-only and did not generate a useful or protein-specific support state. The earlier
   joint-loss screen even reversed that control. A sigmoid relevance weight is not mathematically
   sparse and its short-run behavior is not stable enough for promotion.
2. Sharing one latent representation did not improve endpoint or SAR at this budget. It reduced
   parameters but did not create support-identifiable target-specific variation.
3. The section operator collapsed to a near-zero correction in all arms. Consequently wrong-state
   controls are numerically indistinguishable and cannot demonstrate meta-learning specificity.
4. The level baseline remains stronger than the full model at `k=5` in all arms. The remaining
   performance deficit is in population query shape/calibration and protein-specific section
   geometry, not in insufficient adapter complexity.
5. The validation curve was non-monotone and selected early checkpoints, consistent with a small
   independent-component regime. Expanding capacity or adding MQSA, support dropout, ranking loss,
   geometry loss, or ridge would add axes without resolving the observed specificity failure.

## Stop decision

The research branches remain under `research/meta_fewshot/`. MQSA and pair feedback are not run in
this generation because the prerequisite representation/specificity signal did not pass. The next
valid generation needs either a new source-only supervision signal that improves correct-vs-wrong
protein section geometry or fresh dependency-closed confirmation data. Reusing the consumed test
components for more architecture selection would not supply confirmatory evidence.

# Stage K2 preregistration — three-seed confirmation of the K-REG configuration

Frozen before any K2 arm trained. Date: 2026-08-18. Development evidence;
meta_test sealed and never opened in this stage.

## Selection context (preregistered Stage K arms)

The Stage K screen preregistered two training objectives on the coembedding
branch: K (episodic InfoNCE) and K-REG (positive/negative regression
alignment). K fails G2; **K-REG passes the single-seed gates**:
- k=2 MSE -0.0744 [-0.1413, -0.0070] RESOLVED,
- k=3 MSE -0.0739 [-0.1209, -0.0279] RESOLVED,
- k=5 MSE -0.0465 [-0.0759, -0.0160] RESOLVED,
- k=0 centered MSE -0.0601 [-0.1146, -0.0088] RESOLVED,
- k=0 Pearson +0.1025 [+0.0176, +0.1896] RESOLVED,
- no resolved degradation anywhere; k=0 MSE -0.1369 (unresolved, positive).

The promoted configuration therefore is the K-REG arm: framework innovation
= contrastive coembedding branch (128-dim, on the shared trunk); training
innovation = the positive/negative coembedding regression alignment. The
InfoNCE objective is recorded as tested and inferior.

## Design (frozen)

K-REG exactly as trained in Stage K: 1,200 steps, 3 episodes/step, Stage B
loss recipe + 0.5 * coembedding regression (pos -> 1, neg -> 0), seed
20260815/20260816/20260817, leak-free internal checkpoint selection, GPU
verification. Baseline: T2 at the same three seeds (20260815 from Stage D,
20260816/17 from Stage G2 — all already trained and frozen).

## Gates (pooled across seeds, component bootstrap)

K2-1. At least two of k in {2,3,5} MSE improved with pooled resolved
      intervals (hi < 0) across the three seeds.
K2-2. k=0 centered MSE improved with a pooled resolved interval, and k=0
      MSE is not degraded by a pooled resolved interval.
K2-3. No metric at any k degraded with a pooled resolved interval
      (Spearman, Pearson, CI, cliff sign, level^2, centered, MSE).
K2-4. In every seed: permuted/matched-wrong support MSE above correct and
      wrong-protein MSE above correct (no inversion).

Stop: any gate fails -> the configuration is unconfirmed; no meta_test.
All pass -> freeze architecture/hyperparameters/checkpoints, record the
frozen set, then open meta_test exactly once with a written authorization.
Davis/KIBA independent training follows a later stage only after
confirmation.

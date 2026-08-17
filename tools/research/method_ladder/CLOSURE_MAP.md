# Method-ladder closure map (2026-08-18)

The method-ladder cycle listed eight named families. The 2026-08-17/18
research cycle superseded the ladder; this map records, for each family, the
measured successor stage and its verdict, so the ladder is formally closed.

| family | ladder plan | measured successor | verdict |
|---|---|---|---|
| 1. multimodal representation collapse + basis reallocation | M1 input gate | Stage K/K2 (contrastive coembedding, regression alignment; collapse 0.99859 -> 0.9908) | K-REG: first all-k resolved MSE across 3 seeds; centered gain did not survive pooling -> NOT CONFIRMED |
| 2. Gradient Blending / OGM | M0-M5 | Stage E (orthogonal level/shape routing) + Stage Q (decoupled head) | level/ordering conflict on one shared trunk is fundamental; no routing family tested separates them -> closed by measurement |
| 3. Disentangled Gradient Learning | M0-M5 | Stage E routing ablations | same as family 2 -> closed by measurement |
| 4. attention MIL / Set Transformer / adaptive pooling | M0-M5 | Stage E/J panel set-context heads (mean/max pooling; trained panel heads) | panel heads improve k=0 level but degrade k>=1 ranking (E/J/L/Q) -> closed |
| 5. DrugBAN-style bilinear interaction | M0-M5 | Stage F (pairwise learned transport over embed pairs) | inert vs fixed Tanimoto; fifth learned-kernel family to fail -> REJECTED |
| 6. FS-CAP-style episodic scale | M0-M5 | Stage J (paired cross-target level alignment, BatchDTA-inspired) | paired term added nothing measurable; assay-aware head rejected on ranking -> closed |
| 7. AdaMBind-style task valuation / label-noise robustness | M0-M5 | Stage A/B (inner/outer-loop meta-learning; task selector A2) | selector REJECTED (gradient agreement = redundancy); framework NOT PROMISING then REJECTED -> closed |
| 8. MMP-cliff transformation learning | M0-M5 | Stage F pairwise signed-gap supervision (the measured Stage L2 pairwise direction, r +0.270) | pairwise operator inert vs fixed Tanimoto; cliff-sign never improved by a resolved amount -> closed |

The ladder harness (tools/research/method_ladder/_shared/) is retained as
tooling; no family remains pending.

# Q2c-0 harness self-audit report (2026-08-18)

Preregistration SHA-256: 1027ccde8c8946aa8314ebd7642af89a6abbc3366afd965e8ab43f0da5a26a5c.
Artifact: Q2C0_PROJECTION_AUDIT.json. All numbers below are seed-0 recomputes
that reproduce the frozen X0c artifact where comparable (correct_protein raw
dz 0.455 vs artifact 0.4545; oracle 0.674 vs artifact 0.6742).

## 1. Oracle alignment — DONE
ORACLE_ALIGNMENT_TABLE.md: the diag21 P@U probe (dz 0.68-0.76, cited in the
X0c REPORT) differs from the in-artifact oracle_protein arm in truth seeds
(one shared truth vs three), restart/selection rule (final-state vs 8-restart
val-selection) and checkpoint choice. The X0c REPORT clause citing 0.68-0.76
is withdrawn. Admissible oracle evidence: dz 0.607-0.674, Spearman 0.331-0.390
— below the frozen 0.70 threshold in all three seeds.

## 2. ANOVA projection — DEMOTED and REMOVED from diagnostics

Mechanism of the artifact's 0.51-0.54 projected Pearson on dead heads is now
explained (two defects, neither of them a spurious-signal property of the
projection itself):
(a) the artifact projected the FULL fitted output (mu+pm+lm+head) against
    I+noise instead of the interaction-head output;
(b) the no_interaction_head arm's inter_scale parameter kept gradients during
    training and drifted (raw dz of the 'dead' head after training: 0.457), so
    the artifact's dead head was not dead; it was only zeroed at predict time.

Also: interior eval cells = 0/545 (eval parents AND eval ligands are never
seen in train), so on the eval graph the projection operator degenerates to
subtracting the train grand mean — it removes no main effects at all there.
A dead head plus the true planted main effects reproduces projected Pearson
0.093 (not 0.5), confirming the projection is not manufacturing correlation.

Negative controls (head-level projection, frozen bound |pearson|<0.15 and
permutation p>0.05):
- N1 tau*=0: pearson 0.000, p=1.00 — PASS
- N2 no_interaction_head (truly frozen scale): pearson 0.021, p=0.637 — PASS
- N3 ligand_only: pearson 0.114, p=0.015 — FAIL (absorbs noise ligand-marginal)
- N4 random bipartite graph: pearson -0.064, p=0.169 — PASS
- N5 operator train-only: verified (fit on train mask only; eval labels never enter)
- N6 raw-head dead-zone metric without projection: reported for every arm
- N7 pre/post projection double reporting: in artifact per arm

Verdict: because N3 fails, the ANOVA projection is NOT a safe diagnostic and
is removed from all gates. Primary metrics remain the explicit interaction
head (Spearman / dead-zone sign accuracy / gap vs ligand_only).

## 3. Endpoint / scale findings (new, material)

- Endpoint distortion: on determinate eval cells the latent interaction I and
  the interaction implied by the observable scale (logit of % activity minus
  planted main effects) correlate at only Pearson 0.59 / Spearman 0.555.
  The sigmoid + quantization + censoring pipeline destroys roughly 40% of the
  interaction's rank information at the observable level.
- Minimal linear no-censoring Q2 (identity link, linear head, same graph):
  raw dead-zone sign accuracy 0.52, Spearman -0.075 — no recovery. The linear
  learner cannot read the interaction from the one-hot pocket input even
  without endpoint distortion.
- Truths are saved separately (mu/pmain/lmain/I/z/y/determinate) with SHA-256
  fc9bb8052d7d6a2a4ad311ebfa73a83a51d865af2b740b0a6798da010426f3f3.

## 4. What Q2c-0 establishes

1. The X0c anova_projection numbers are explained by an ill-defined metric
   (full-yhat projection) plus a non-frozen negative-control arm; they do NOT
   show that the planted-signal harness itself produces spurious interaction
   recovery.
2. The remaining Q2 failure decomposition now has measured components:
   endpoint distortion (rank corr 0.59) and one-hot pocket unreadability for
   linear learners (dz 0.52); the oracle MLP arm (dz 0.674) is consistent with
   an additional representation/optimization gap of unknown size.
3. Open question carried into Q2c-1: does an oracle MLP on the z-scale
   (no sigmoid, no censoring) reach dz >= 0.70? If yes, the endpoint pipeline
   is a primary blocker; if no, graph power / training protocol is.

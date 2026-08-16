# Stage 2 preregistration: interaction-grammar discriminator

Written before the Stage 2 runs. One changed variable: **the architecture**.
Seed, budget, optimizer, losses, sampler and evaluation banks are identical
across arms.

## Arms

| arm | `--arch` | model |
|---|---|---|
| control | `bpsf` | retained `QPSMPBioModel` (pair trunk + label-locked residual kernel) |
| candidate | `grammar` | `InteractionGrammarModel` (atom-to-residue contact grammar + transferability-gated transport) |

Matched settings: seed 20260812, 60 steps, 4 episodes/step, hidden 192,
task 48, pair/embed 96, latents/contact types 24, heads 8, ligand layers 4,
query 12-20, lr 3e-4, backbone scale 0.25, binding weight 1.0, warmup 0.

## Evaluation

* **Frozen protocol bank**: `evaluation_seed=73101`,
  `eval_targets_per_component=1` (6 episodes per k). Retained comparator.
* **Wide bank**: identical construction with all eligible meta-test targets
  (50 targets, 7 components). Primary Stage 2 discriminator, because the Stage 0
  audit showed the 6-episode bank cannot resolve differences below ~0.05 MSE.

## Gates (all evaluated on the wide bank)

A smoke run can reject but cannot admit. Continue to Stage 3 only if **R1 and
R2 pass and none of R3-R6 fails**.

### Amendment, recorded before any Stage 2 run produced data

R1 was originally "zero-shot k=0 MSE at least 0.10 lower than the control arm".
It was calibrated on the assumption, drawn from the real-corpus audit, that the
retained pair trunk is capacity-limited for protein-conditioned interaction.
`stage0_audit_20260815/TRUNK_CAPACITY*.json` falsified that assumption: on a
noiseless synthetic protein-by-ligand bilinear task, the retained trunk reaches
relative held-out MSE 0.0028-0.0085 at learning rates 3e-4 and 1e-3, and only
collapses to a constant at 3e-3 and above. Both trunks can express the
interaction; the retained trunk is 4.4x more expensive per step and has a
narrower stable learning-rate band.

Since neither arm is converged after 60 steps, an absolute k=0 MSE threshold at
that budget measures noise. R1 is therefore replaced by a
degeneracy test (R1) plus an explicit cost test (R1b), which is what a 60-step
discriminator can actually resolve. R2-R6 are unchanged.

| id | requirement | control value to beat |
|---|---|---|
| R1 | zero-shot k=0 MSE no worse than the control arm, and zero-shot spread across queries inside an episode > 0.20 pK | control arm same budget; retained baseline spread 0.065 pK |
| R1b | seconds per training step lower than the control arm | control arm, same budget |
| R2 | wrong-protein zero-shot MSE gap >= 0.05 | 0.0119 (retained baseline) |
| R3 | k=1 `sar_cut_mse_pk != full_mse_pk` (query-specific channel live) | identically equal in the control |
| R4 | magnitude-matched permuted MSE > full MSE at every k in {1,2,3,5} | — |
| R5 | no zero-gradient trainable tensor at any k | 17-23 dead tensors in the control |
| R6 | k in {2,3,5} full MSE not worse than the control arm by more than 0.02 | control arm, same budget |

## Rejection handling

If R1 or R2 fails, the mechanism is rejected: `--arch grammar` stops being a
candidate, `RESULT.json` and this file are retained, and the source change is
reverted before the next hypothesis is tested.

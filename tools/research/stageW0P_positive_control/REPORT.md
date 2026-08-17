# Stage W0-P report — local point-mutant positive control: FAILED

Preregistration SHA-256:
`ba0b51ec419b0275a129e69e4cb45db1bccbdd138000893ee7daf881e7bacbf1`.
Artifacts: `W0P_PANEL.json`, `W0P_RESULT.json`. No production model was
trained; the only model is the frozen low-capacity bilinear diagnostic.

## Panel

6 near-identical BindingDB sequence pairs (1–5 residue mismatches, >=3 shared
ligands), 32 ligand rows. One pair has 7 ligands with all-positive delta
(mean +1.42 pK); the other five are near zero.

## Test

Low-capacity bilinear `phi(ligand) * psi(protein-change)` trained with AdamW,
300 steps, **leave-one-pair-out** (the earlier leave-one-row-out version was
discarded because it let every arm memorise pair identity; this is recorded,
not hidden). Three seeds.

| arm | Pearson | Spearman | sign accuracy | MSE |
|---|---:|---:|---:|---:|
| correct positions | -0.046 | +0.117 | **0.240** | 0.606 |
| random positions | -0.143 | +0.090 | 0.156 | 0.610 |
| BLOSUM-approximate unrelated positions | -0.187 | -0.034 | 0.125 | 0.609 |
| global ESM pooled difference | -0.103 | +0.193 | **0.760** | 0.609 |
| random protein | +0.055 | +0.226 | 0.625 | 0.597 |
| ligand-only (structural zero) | 0.0 | 0.0 | 0.0 | 0.605 |

## Verdict

**W0-P FAILS.** The correct mutation positions do not beat the corrupted
position arms, and the global-pooled arm has the best sign accuracy. Per the
frozen positive-control rule this is a rejection: the pipeline does not
recover the signed point-mutant change with the required specificity, and no
subsequent cold-protein null may be interpreted as biological absence.

The global-pooled sign accuracy (0.760) is recorded as an unexplained,
under-powered observation, **not** a positive result: with 6 pairs it cannot
be separated from a lucky continuous feature mapping.

## Consequences

- W1 remains **NO-GO**.
- Before retrying W0-P: (a) enlarge the point-mutant panel, (b) use
  leave-one-pair-out and report pair-level bootstrap, (c) test a stronger but
  still gradient-trained representation, (d) pre-register the pass rule from
  the enlarged panel's effective sample size.

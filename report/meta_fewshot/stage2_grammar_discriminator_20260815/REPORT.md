# Stage 2 result: interaction-grammar discriminator

Matched seed 20260812, 60 steps, 4 episodes/step, identical optimizer, losses,
sampler and evaluation banks. One changed variable: `--arch`.

Numerical authority: `arm_bpsf/RESULT.json`, `arm_grammar/RESULT.json`,
`WIDE_arm_bpsf.json`, `WIDE_arm_grammar.json`.

## Wide bank, 42 episodes over all eligible meta-test targets (primary)

| arm | k=0 | k=1 | k=2 | k=3 | k=5 |
|---|---:|---:|---:|---:|---:|
| `bpsf` | 4.085 | 2.631 | 2.146 | 1.866 | 1.611 |
| `grammar` | **3.996** | **2.204** | **1.829** | **1.616** | **1.372** |
| difference | -0.089 | -0.428 | -0.317 | -0.250 | -0.239 |

## Frozen protocol bank, 6 episodes

| arm | k=0 | k=1 | k=2 | k=3 | k=5 | wall s | peak MB |
|---|---:|---:|---:|---:|---:|---:|---:|
| `bpsf` | 2.497 | 1.840 | 1.391 | 1.220 | 1.206 | 405 | 6,053 |
| `grammar` | **2.329** | **1.643** | **1.296** | **1.169** | **1.123** | 350 | **1,698** |

## Gate outcome

| id | requirement | outcome |
|---|---|---|
| R1 primary | k=0 no worse than the control | **pass** (3.996 against 4.085) |
| R1 secondary | zero-shot query spread > 0.20 pK | **fail** (0.107 pK) |
| R1b | seconds per step lower | **pass** (4.70 against 5.57) |
| R2 | wrong-protein zero-shot gap >= 0.05 | **fail** on the wide bank (-0.019); pass on the frozen bank (0.052) |
| R3 | k=1 query-specific channel live | **pass** (`sar_cut` 2.538 against `full` 2.204) |
| R4 | permuted worse than correct at every k | **pass** |
| R5 | no dead trainable branch | **pass** (Stage 1 gate) |
| R6 | k in {2,3,5} not worse than the control | **pass** (0.24-0.32 better) |

## Decision

The **few-shot mechanism is not rejected**: it dominates the control at every
k on the primary bank, and it is the first configuration in this lineage whose
k=1 correction is query-specific rather than a scalar level shift.

The **zero-shot claim is rejected at this budget**: the endpoint is still close
to a constant (0.107 pK spread against a 0.93 pK label spread) and protein
specificity is not established on the wide bank.

R1-secondary and R2 are carried into Stage 3 and Stage 4 as blocking gates for
any zero-shot or protein-specificity statement. Stage 3 targets the zero-shot
degeneracy directly.

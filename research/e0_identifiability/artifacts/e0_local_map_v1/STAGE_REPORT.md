# P1R2B-E0 Local Mechanistic Affinity Potential

Updated: 2026-08-06

Decision: `SYNTHETIC_PRE_GATE_FAIL_CLOSED; E0_S_NOT_RUN; DAVIS_NOT_ACCESSED;
T/P2-P4_FROZEN`.

## Registered Test

E0 freezes the complete P1B frontend and applies a width-32 nonlinear map to
atom-local state, atom chemistry, residue-local state, residue chemistry,
contact probability and five distance probabilities before contact-weighted
pooling. Global protein/ligand embeddings and identity fields are prohibited
MAP inputs. The output is a Mechanistic Affinity Potential, not physical
binding free energy.

The synthetic control was frozen at seed 17, eight tasks per existing fold,
20 ligand states per task, folds 0-3 training, fold 4 holdout, 60 epochs,
AdamW learning rate 1e-3, weight decay 1e-4 and batch size 4. PASS required
correct CI >=0.80, correct-minus-ligand >=0.10, correct-minus-deranged >=0.10
and atom-permutation error <=1e-6.

## Label-Blind Inputs

The input audit discarded affinity fields while parsing. It found 154,165
governed rows; 152,934 satisfy the frozen 128-atom and standard-residue input
contract. Exclusions were 1,111 oversize-ligand rows and 120 unsupported-residue
rows. The retained floor is 3,783 tasks, 681 proteins, 245 closure components
and 573 tasks in the smallest fold, which passes the D0/D1 floor.

The label-free model manifest contains 152,737 rows, 3,783 tasks, 93,761 full
InChIKey ligand states and 681 proteins. Full InChIKey is required for cached
state identity because connectivity keys merge stereoisomers.

## Frozen Cache

The exact P1B checkpoint SHA is
`90b0010b81fa2758a2dbdd1a8dbe06adae2e05acbbc267ccb62ceee6ff6c4f37`.
The exact ESM model revision is
`a695f6045e2e32885fa60af20c13cb35398ce30c`. Six protein and 46 ligand shards
occupy 637,642,084 bytes. Pooled ligand states are marked ligand-prior-only;
MAP access is prohibited.

## Synthetic Result

| Arm | Holdout CI |
|---|---:|
| Ligand | 0.48454 |
| Correct MAP | 0.68553 |
| Deranged MAP | 0.64934 |

Correct-minus-ligand was `+0.20099` and passed. Correct-minus-deranged was
`+0.03618` and failed the registered `+0.10` threshold. Correct CI failed the
registered `0.80` threshold. Atom-permutation error was `0.0` and passed.

The model learned synthetic incremental signal and is permutation invariant,
but partner-specific held-out recovery is insufficient under the frozen
control. This is not a real-affinity result.

## Stop Decision

No real ChEMBL affinity value was read. E0-S closure-OOF training did not run.
DAVIS metaval and recipient labels were not accessed. The synthetic run is not
eligible for a hyperparameter retry under its registration. Any redesign of
the synthetic control, or authorization of typed interaction T, requires a new
decision. P2-P4 remain frozen.

Key SHA-256 values:

- input manifest: `2c45a3684ffb70c59dc4979e5121c95af6ded03c767a6c5b799d82c076877467`
- ESM manifest: `a08744bce119e2db0e64b1a8a229d51630542e5564d2d7b299b6aeedeee1e0a4`
- local-state manifest: `4a924c67ec4ed58f31b7478b2773cfecc1d637f569d0d2297ca9843979c75d50`
- input audit: `cf87188cd362e976a4753fe526742fdebea68edbc01285c2e0d95340a8c12e6f`
- synthetic gate: `2688285ffdab6b0d9f6a71e17d6d5d42de0cc69fe6cefb59259589474af297b9`
- MAP implementation: `471e0a33d9ad4488c73982c0f325a799155ab12870ec28b11a565a80abdff18f`
- synthetic runner: `b98fe14313b11512219e113f814c6901f708f14a6279c13c5ad47c9203ac668c`
- E0 tensor contract: `751c08eaa95ca59993ef30917bf51270aa30ff26eca4c1adc82ef6f37ea5cc55`

Verification: the complete `drug` regression suite passed
`205 passed in 55.43s`.

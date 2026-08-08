# T-DIR-P0 Result

Date: 2026-08-07

## Verdict

```text
PILOT_LEARNABILITY_SIGNAL_NOT_OBSERVED
```

This was a `24/8/8` structure-only feasibility pilot. It used 40 distinct
homology groups and scaffold-disjoint selected splits, PLIP 3.0.1 weak labels,
frozen P1B features and fixed logistic probes. ChEMBL affinity, DAVIS and
recipient-label reads were all zero.

## Primary Result

All `40/40` selected structures completed the executed annotation pipeline.
The hydrophobic primary channel was evaluable with `53/26/24` positive pairs
from `15/7/6` train/validation/test complexes among `17,556` oracle-near
candidate pairs.

| Arm | Train AP | Validation AP | Test AP |
|---|---:|---:|---:|
| D0 contact + distance | 0.00903 | 0.01906 | 0.00987 |
| D1 + atom/residue chemistry | 0.03954 | 0.04445 | 0.03735 |
| D2 + frozen local states | 0.25496 | 0.01699 | 0.01996 |

Test prevalence was `0.00620`. D2 exceeded prevalence and D0 on test, but it
failed the registered `AP - prevalence >= 0.10` condition and was below D0 on
validation. The train-to-held-out collapse is consistent with high-dimensional
small-complex overfit. Because arm capacity is not matched, this is a
description rather than a causal attribution to local-state biology.

Protein-donor H-bond was also evaluable. Its D2 test AP was `0.01344` at
prevalence `0.00528`; this is a small descriptive signal, not a formal PASS.
Ligand-donor H-bond, halogen, salt, pi-stack and cation-pi did not satisfy the
fixed prevalence/support criteria or had invalid group mapping.

## Post-Run Triage

The executed primary hydrophobic mapping covered `144/147 = 97.96%` PLIP
events. Directional H-bond mapping covered more than 98% of events. A recursive
namedtuple bug prevented all group-based salt/pi/cation-pi ligand atom sets from
mapping; code and a regression test were corrected after the run, but the same
test panel was not rerun.

Two additional preregistration deviations were found: shuffle controls and a
full Open Babel ligand-chemistry round-trip audit were not implemented. These
omissions cannot rescue the already negative primary result, but they prevent a
claim that the complete annotation contract passed. The original report remains
immutable; `postrun_audit.json` binds the executed and corrected source hashes.

## Scientific Meaning

The run shows that the open-coordinate -> PLIP -> canonical atom/residue ->
frozen-P1B feature pipeline is operational for direct hydrophobic and H-bond
events. It does not show replicated typed-interaction learnability: the primary
D2 signal did not satisfy the frozen feasibility criteria and overfit strongly.

This result does not identify affinity energetics, interaction energy, novel
target/scaffold generalization or a biological coordinate admissible to `z`.
Nothing moves to `model/` or normal `scripts/`. A future attempt requires a new
sealed panel, complete mapping/chemistry controls, lower-dimensional named
features or materially more independent complexes, and a separately registered
formal T-DIR Gate.

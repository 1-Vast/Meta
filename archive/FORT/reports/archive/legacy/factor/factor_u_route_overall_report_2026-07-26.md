# FACTOR-U strict unlabeled expansion: overall report before reopening exploration

Date: 2026-07-26

## Why this route was run

FACTOR-C failed at representation identifiability, not at proven real chemical support. The user's
next honest option was to expand unlabeled chemistry under a global evaluation firewall. This route
did not consume an agent-proposed candidate slot.

## U0: PLINDER-only corpus

After deleting every KIRHub, Reinecke and Papyrus-Christmann connectivity/scaffold, PLINDER retained
29,382 molecules and 17,003 nonempty scaffolds with complete element/role and size-domain support.
It failed only because its largest scaffold contributed 1.266%, above the frozen 0.5% maximum.
Verdict: `FACTOR_U0_CORPUS_INELIGIBLE_STOP`.

## U0-B: mixed public corpus

Without relaxing that threshold, adding predicate-filtered ChEMBL-37 `train` structures produced
138,805 molecules and 60,021 scaffolds. The maximum scaffold fraction fell to 0.367%; 26,682
PLINDER-exclusive and 109,423 ChEMBL-exclusive structures remained. All eligibility and firewall
gates passed. Verdict: `FACTOR_U0B_PASS_AUTHORIZE_U1_PREREGISTRATION`.

## U1-L: 50,000-scaffold local mechanism proof

The exact F0-C1 masked GINE was trained once on 45,000 globally disjoint scaffolds and selected by
5,000 unlabeled validation scaffolds. It reconstructed local chemistry strongly and encoded
pharmacophore roles, but the same failure mode remained in every evaluation fold:

| Diagnostic | F0-C1 | U1-L | Gate |
|---|---:|---:|---:|
| participation rank | 8.15--10.05 | 8.77--8.91 | >=16 |
| rank ratio | 0.064--0.078 | 0.0685--0.0696 | >=0.10 |
| atom false-decoy coverage | 0.0660--0.0748 | 0.0670--0.0742 | <=0.05 |
| inner median coverage | 0.7735--0.8203 | 0.7345--0.7947 | >=0.85 |
| external median/q10 | 0.6710/0.4789 | 0.6231/0.4078 | 0.90/0.70 |

Verdict: `FACTOR_U1L_REPRESENTATION_UNIDENTIFIED_STOP`.

## Root cause and closed options

Corpus availability is not the current bottleneck: a large, licensed, globally firewalled and
scaffold-diverse corpus exists. The ligand-only masked reconstruction objective compresses
chemistry into approximately nine effective directions that are excellent for common structural
and pharmacophore labels but unsuitable as a calibrated functional-substitution atlas. More
unlabeled structures strengthen reconstruction without supplying the missing
protein-conditioned SAR geometry.

Therefore:

- do not run full 138,805-molecule or uncapped-pair replication on the larger machine;
- do not tune width, rank penalties, epochs, seed, bandwidth or gates;
- do not authorize F1-C or an affinity model;
- do not interpret low coverage as proven lack of real chemical support;
- keep the existing ChEMBL confirmation partition quarantined.

The overall category remains **② `SIGNAL_PRESENT_EVIDENCE_INSUFFICIENT`**: true-minus-decoy and
structural/role signals are strong, but they are not the transferable protein-conditioned
interaction signal needed by the task.

Only a substantively different supervision source is now justified: repeated ligand
transformations measured across proteins/documents, explicit residue--atom interaction labels with
target-necessity controls, or the already designed prospective factorial panel. Another
ligand-only carrier representation is not justified.

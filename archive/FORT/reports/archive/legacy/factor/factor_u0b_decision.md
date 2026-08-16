# FACTOR-U U0-B decision

Date: 2026-07-26  
Verdict: `FACTOR_U0B_PASS_AUTHORIZE_U1_PREREGISTRATION`

The mixed PLINDER plus ChEMBL-37-train unlabeled corpus passed every preregistered eligibility gate
after the global connectivity and Murcko firewall against KIRHub, Reinecke and
Papyrus-Christmann.

- 138,805 unique retained molecules;
- 60,021 distinct nonempty scaffolds;
- largest scaffold 510 molecules, 0.3674%, below 0.5%;
- 26,682 PLINDER-exclusive, 109,423 ChEMBL-train-exclusive and 2,700 shared molecules;
- all evaluation element classes and pharmacophore roles supported;
- retained heavy-atom q01/q99 10/66;
- zero connectivity and scaffold overlap with every evaluation source.

U0-B only authorizes a separately preregistered U1 representation pilot. It does not authorize F1-C
or any affinity model.

PLINDER projected only its seven frozen ligand structure/quality fields. ChEMBL used Parquet
predicate pushdown for `dual_cold_split == "train"` and materialized only `conn/scaffold`;
confirmation rows materialized = 0 and no activity/affinity/protein field was projected. The
historical ChEMBL confirmation partition remains permanently quarantined.

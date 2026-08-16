# FACTOR-U U0 decision

Date: 2026-07-26  
Verdict: `FACTOR_U0_CORPUS_INELIGIBLE_STOP`

The PLINDER-only strict unlabeled expansion passed every scale, firewall and domain gate except the
frozen scaffold-concentration gate.

- 29,382 unique molecules retained after deleting every KIRHub, Reinecke and
  Papyrus-Christmann connectivity or Murcko scaffold;
- 17,003 distinct nonempty scaffolds;
- zero connectivity and scaffold overlap with every evaluation source;
- all evaluation element classes and pharmacophore roles supported;
- retained heavy-atom q01/q99 7/58, covering evaluation q05/q95;
- largest scaffold 372 molecules, or 1.266%, versus the frozen maximum 0.5%.

The size result is encouraging, but U1 was explicitly contingent on all U0 gates. The PLINDER-only
corpus therefore cannot authorize training, and no after-result scaffold downsampling or gate
relaxation is allowed.

This failure is a corpus-composition result, not a memory limit and not an encoder result. A broader
public-unlabeled source mixture is scientifically distinct and remains allowed under the user's
post-F0-C data route, provided it receives a new preregistration and preserves the same global
evaluation firewall.

U0 projected only seven ligand structure/quality fields; activity, affinity, protein, pocket and
PLINDER split fields were unread. The existing ChEMBL confirmation partition remains permanently
quarantined (`project_historical_confirmation_labels_read=true`);
`current_run_confirmation_labels_read=false`; `sealed_test_consumed=false`.

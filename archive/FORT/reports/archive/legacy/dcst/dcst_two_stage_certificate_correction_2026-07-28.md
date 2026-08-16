# DCST source-certificate correction

Date: 2026-07-28  
Applies before: any 4,000-step accepted DCST run  
Supersedes only: the per-rank-one certificate unit in
`dcst_two_stage_preregistration_2026-07-28.md`

## Trigger

The 500-step engineering pilot produced one source spectral component with
held-source residual Spearman 0.2498, but its wrong-target Spearman was also
0.2498. This equality is structural, not a threshold or optimization issue.

For a rank-one term

```text
g_j(t,d) = s_j (z_t^T u_j) (z_d^T v_j),
```

fixing `t` leaves a constant scalar times one global ligand projection.
Replacing the target changes only that scalar. Within-target Spearman is
unchanged whenever the sign does not flip. Consequently a single rank-one
term cannot support the preregistered target-destruction certificate.

No 4,000-step ChEMBL result existed when this defect was found. The correction
is based on the source-only pilot and the algebra above. The downstream MDE,
losses, data, firewall, encoder width, rank, seeds, arms and pass gates remain
unchanged.

## Corrected certificate unit

The rank-eight SVD is partitioned before scoring into four non-overlapping
ordered spectral bands:

```text
(1,2), (3,4), (5,6), (7,8).
```

Each band's score is the sum of its two rank-one terms. A target changes the
relative mixture of two ligand projections, so wrong-target destruction can
now alter a within-target ordering. True, wrong-target and ligand-deranged
utilities, margins and confidence mapping are otherwise identical to the
original preregistration.

One Stage-2 gate is shared by both components in a band. A band is inactive
unless its source-only destruction margin is positive. The complete
`FullTransferResidual` arm assigns confidence one to all four bands, and
`DCST-CertShuffle` permutes the four band confidences.

This is an identifiability correction, not a relaxation: a two-direction band
must still beat both destruction controls, while the invalid rank-one test
could pass or fail primarily through a target-dependent sign flip.


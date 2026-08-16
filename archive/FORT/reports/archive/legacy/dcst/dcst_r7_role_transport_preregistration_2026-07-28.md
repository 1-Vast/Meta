# DCST-R7 content-addressed role transport preregistration

Date: 2026-07-28  
Status: frozen before implementation and training

## Hypothesis and innovation

R7 adds one constrained operator to R6: **Content-Addressed Role Transport
(CART)**. The SIFTS-aligned structural head still predicts
`P_phi(segment, interaction_type)` on 32 sequence segments. A shared router
assigns each segment to one of eight latent protein roles using only that
segment's frozen ESM content embedding:

```text
A_phi(role | segment_content)
P_phi(role, type) =
    sum_segment P_phi(segment, type) A_phi(role | segment_content)
g(target, ligand) = <Theta_role, P_phi(role, type)>.
```

The router never receives the absolute segment index or the learned positional
embedding. `P_phi(role,type)` remains normalized and `Theta_role` is an
`8 × 8` exact-null interaction-energy matrix.

This makes transfer equivariant to the source protein's arbitrary normalized
segment numbering while retaining segment-local structural supervision. The
eight learned roles are shared across proteins and are addressed by ESM
content, so a certified direction denotes a content-defined protein role ×
interaction-type energy rather than "relative sequence segment 19".

The nearest recent cold-start methods align entity or link representations,
and fine-grained binding pretraining predicts atom/residue interactions. The
registered novelty claim is narrower: exact-entity structure-privileged
pretraining, content-addressed transport of the predicted interaction measure,
and held-source destructive spectral certification as the only cross-stage
affinity interface.

## Frozen architecture and budget

- eight content roles; no role-count search;
- same exact UniProt/SIFTS source rows, 32 segments, source firewall, ESM and
  Morgan inputs, seed 1729, optimizer, and 4,000-step budgets as R6;
- unchanged segment-level absolute, centered, counterfactual, and retrieval
  privileged losses;
- matched CART-NoPriv with identical router and capacity but no privileged
  structural labels;
- unchanged two-direction SVD bands, certificate rule, mechanism thresholds,
  source base, and split policy.

The existing R6 absolute-position SMB is retained as a historical matched
architectural control. No R6 development value is used to select the role
count, losses, learning rate, or gate threshold.

## Source gate

Before another downstream training run:

1. the segment-level privileged mechanism probe passes both registered
   destruction margins;
2. CART certifies at least one held-source role×type spectral band;
3. CART certifies strictly more bands than CART-NoPriv.

Failure stops R7 before another ChEMBL affinity-label load.

## Conditional Stage 2

If the source gate passes, reuse the already frozen and hash-verified ChEMBL
32-segment cache. The source segment head, content router, certified spectral
directions, and their destruction identities are frozen in the transfer
branch. The existing ChEMBL train/development arms, paired bootstrap, MDE,
negative-transfer audit, RMSE gate, target/ligand destruction, confirmation,
and sealed policies remain unchanged.

R7 must pass the complete pre-existing Stage-2 gate; confirmation remains
unauthorized otherwise.


# DCST-R6 shared mechanism bottleneck preregistration

Date: 2026-07-28  
Status: frozen before implementation and training

## Hypothesis and innovation

R6 tests a constrained cross-stage interface, the **Shared Mechanism
Bottleneck (SMB)**. Given an exact UniProt target and ligand, the Stage-1
network predicts

```text
P_phi(segment, interaction_type | target, ligand)
```

over the SIFTS-aligned `32 × 8` structural grid. The affinity interaction
residual is restricted to

```text
g(target, ligand) = <Theta, P_phi(target, ligand)>.
```

`Theta` is a `32 × 8` interaction-energy matrix with the exact regular null
`Theta=0`. Thus Stage-1 privileged supervision and affinity supervision act
through the same pair-specific object; affinity cannot bypass the structural
map through separate pooled vectors.

This is materially different from a generic pretrained encoder or adapter:
the transferred directions have an explicit segment-by-interaction-type
meaning, are available for unseen targets and ligands from ESM segments and
Morgan environments alone, and can be destroyed and certified before
downstream training.

## Stage 1

Use the unchanged exact-UniProt/SIFTS PLINDER firewall, seed 1729, 32 target
segments, source base, 4,000 base steps, and 4,000 interaction steps. Train
`P_phi` jointly with the source affinity readout and the already frozen
absolute, centered, counterfactual, and bidirectional privileged losses.

Matched `SMB-NoPriv` uses the same architecture, initialization, affinity
labels, optimizer, and steps but receives no privileged structural labels.
No decoupled teacher or random-frozen control is used: R5 already rejected
that interface.

Decompose `Theta` by SVD. Each pre-existing two-direction spectral band is
scored additively on held-source exact targets. For a rank-one direction
`k`, its per-pair contribution is

```text
s_k * sum_(segment,type)
    P_phi(segment,type) * u_k(segment) * v_k(type).
```

Target destruction recomputes `P_phi` using the registered exact-target
derangement; ligand destruction applies the registered within-target
derangement. Certificate scale, positive-utility rule, and all mechanism
thresholds are unchanged.

## Frozen source gate

Before any new ChEMBL affinity-label load:

1. the SMB privileged joint-map probe passes the existing absolute and
   centered target- and ligand-destruction gate;
2. at least one held-source SMB spectral band is active;
3. SMB certifies strictly more bands than `SMB-NoPriv`.

Failure stops R6 at Stage 1.

## Stage 2, conditional on source pass

The certified Stage-1 map generator and spectral directions are frozen. Only
their registered group gates plus a downstream SMB residual may train on
ChEMBL-37 strict dual-cold train data. Existing B0, scratch, naive fine-tune,
frozen-map, full-transfer, certificate-shuffle, and no-privileged controls;
paired target-level bootstrap; MDE; RMSE; destruction; confirmation; and
sealed-test rules remain unchanged.

This preregistration changes only the cross-stage affinity interface and the
mathematically corresponding spectral component calculation. It does not
change data access, outcome labels, thresholds, seed, training budget, or
split policy.

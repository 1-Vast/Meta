# Preregistration XP5 — Fixed Named Typed Interaction Basis

Registered: 2026-08-08, before any XP5 arm was scored. All prior verdicts stand.

## 0. Which failed assumption this addresses

XP2 closed the *free latent section* route. Its failure mode was specific and
measured: the query ligand's loading `u(L)` had to be **predicted** from
chemistry, and although that prediction is genuinely informative
(`R2 = +0.199`), the noise it injects into `<uhat, v>` destroys the
target-specific component faster than the shared one. Under double closure the
derangement control was indistinguishable from the correct protein and random
ligand features reproduced the whole gain.

XP4 then closed the *many-independent-panels* route for a different reason: in
the BindingDB curated-articles design the within-panel interaction sits **below
the measurement noise floor** (per-report `sigma = 0.777`, `gamma` sd `0.406`,
within-panel chemistry-neighbour ceiling `R2 = -0.539`).

XP5 tests the one remaining rung that changes an actual assumption rather than
adding capacity: on the panel that **does** contain reproducible interaction
(Metz), replace the free latent basis with a **fixed, biologically named basis
that is computed analytically from both partners** and therefore injects no
prediction noise and carries no gauge ambiguity.

```
XP1/XP2:  gamma_hat = < uhat(L) , v(S) >      uhat is PREDICTED, latent, gauge-free only up to GL(d)
XP5:      gamma_hat = sum_c  w_c * x_c(P,L)   x_c is COMPUTED, named, fixed before fitting
```

## 1. Hypothesis

> A small fixed set of named physicochemical complementarity channels
> `x_c(P,L)`, each computed analytically from the aligned binding-site residues
> of `P` and the structure of `L`, predicts the protein-by-ligand interaction
> residual on held-out kinase groups **and** held-out ligand scaffolds, beating
> capacity-matched random channels and every pairing control.

## 2. Data

`BLK-METZ-XP2` exactly as frozen in `PREREG_XP2.md` §4: 928 compounds x 147
kinases, 32,849 measured cells, 258 scaffold components, 8 KLIFS groups, panel
SHA-256 `beb97a0e125ccabf...`. Metz is a **development panel** at this point; a
pass here would require external replication before any claim, and is registered
as such.

## 3. The basis (fixed and named BEFORE any fitting)

For each kinase, the KLIFS 85-residue aligned pocket gives a residue multiset.
For each ligand, RDKit gives structural counts. Each channel is a product of one
protein-side and one ligand-side quantity, both standardised on training cells:

| `c` | channel | protein side | ligand side | mechanism |
|---|---|---|---|---|
| 1 | H-bond donor complementarity | pocket acceptor-capable residues (D,E,N,Q,S,T,Y,H) | ligand HBD count | donor-acceptor pairing |
| 2 | H-bond acceptor complementarity | pocket donor-capable residues (K,R,N,Q,S,T,Y,W,H) | ligand HBA count | acceptor-donor pairing |
| 3 | electrostatic, anionic ligand | pocket net positive charge (K,R minus D,E) | ligand negative formal charge | ionic |
| 4 | electrostatic, cationic ligand | pocket net negative charge | ligand positive formal charge | ionic |
| 5 | hydrophobic burial | pocket hydrophobic fraction (A,V,L,I,M,F,W,C) | ligand cLogP | hydrophobic effect |
| 6 | aromatic / pi stacking | pocket aromatic fraction (F,W,Y,H) | ligand aromatic ring count | pi-pi |
| 7 | steric gatekeeper fit | gatekeeper (KLIFS position 45) side-chain volume | ligand heavy-atom count | back-pocket access |
| 8 | steric hinge fit | mean side-chain volume at hinge positions 46-48 | ligand rotatable bonds | hinge accommodation |
| 9 | polar surface complementarity | pocket polar fraction | ligand TPSA | desolvation proxy |
| 10 | size complementarity | pocket total side-chain volume | ligand molecular weight | overall fit |

KLIFS position 45 = gatekeeper and 46-48 = hinge were verified against
CDK2/ABL1/EGFR/SRC in XP1 (DFG at 81-83; CDK2 Phe80 maps to position 45). Ten
channels, all deployment-observable from sequence plus SMILES, all named before
fitting, none latent.

## 4. Arms

| arm | definition |
|---|---|
| `Z0` | `gamma_hat = 0` |
| `TYPED` | ridge on the 10 named channels, fitted on training cells |
| `RAND-C` | 10 capacity-matched random channels, same product form |
| `PERM-P` | protein side of every channel permuted across proteins in training |
| `FOREIGN-P` | at test, channels computed against a deranged held-out protein |
| `SHUF-PAIR` | protein-ligand correspondence shuffled at test |
| `LATENT-REF` | the XP2 `SEC` arm at `k=5`, for reference only |

Predictions are double-centred within the held-out block before scoring, so no
arm can win by reintroducing a main effect.

## 5. Split, metric, inference

Double held-out: KLIFS **group** closure x scaffold-component closure, the same
folds as XP2-D. Metric `R2_gamma = 1 - SSE(arm)/SSE(Z0)` on the interaction
residual computed from training cells only. Intervals are the **wider** of a
protein-component and a scaffold-component cluster bootstrap, 2,000 resamples.

## 6. Gate (frozen, same floors as XP2)

1. `R2_gamma >= 0.05`, CI lower bound `> 0.02`;
2. `Delta` vs `RAND-C` CI lower bound `> 0`;
3. `Delta` vs `FOREIGN-P` CI lower bound `> 0`;
4. `Delta` vs `PERM-P` CI lower bound `> 0`.

One registered configuration; ridge strength by nested group-CV inside training.
No sweep. Thresholds fixed. A positive but negligible effect
(`R2_gamma < 0.05`) is a failure.

## 7. What each outcome means

| outcome | verdict |
|---|---|
| Gate passes | `TYPED_BASIS_INTERACTION_IDENTIFIED`; proceed to `k <= 5` section and external replication |
| `TYPED` ~ `RAND-C` ~ 0 | the deployment-observable pose-free typed basis does not carry the interaction; combined with XP2 and XP4 this closes the sequence-plus-2D rung, and only structure-dependent rungs `B2`/`B3-with-pose` remain |

No XP5 outcome authorises DAVIS access, `model/` promotion, or modification of
the frozen theory, CSMO, Band, `K` or the mesh.

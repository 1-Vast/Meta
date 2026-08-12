# Preregistration: quotient training corpus and linear T-BASIS witness

Stage: `E-AFF-CQ-R2_BINDINGDB_TBASIS_LINEAR_WITNESS`

Status: frozen after CQ-R1 development PASS and before split metrics, T-BASIS
feature generation or model fitting.

## Scope

Primary training source is exact uncensored Ki from BindingDB-curated Articles
202608. Kd, patents, PDSP, ChEMBL, profiling panels, DAVIS and recipient data are
not training inputs in this stage.

## Frozen biological instrument

Restore the previously validated T-BASIS statistic:

```text
8 ligand chemistry classes x 6 residue chemistry classes x 6 radial RBFs
= 288 coordinates
```

The P1B checkpoint, ESM2-t30 provider, ligand graph schema, residue-slot policy,
RBF centers, radial calibration and normalization are frozen. First-stage
prediction is only:

```text
q(P,L) = w.T @ phi_288(P,L)
```

with 288 trainable coefficients. No new PLM, GNN, attention, adapter or named
four-channel mechanism is added. The 288D statistic remains research-only.

## Identity and closure

Build one immutable cell table after replicate averaging. Define dependencies by:

- source document;
- protein sequence clusters at global identity >= 0.40 using the frozen
  Needleman-Wunsch/BLOSUM62/gap-open-10/gap-extend-1 implementation;
- stereo-aware ligand connectivity and Bemis-Murcko scaffold.

Panels sharing any dependency belong to one conflict component. Components are
assigned whole to train or development by a deterministic SHA-256 ordering;
no row or panel may cross a split. If strict union produces a giant component,
report it rather than pretending quotient coordinates are independent. A later
DataSAIL-style discard split may be separately registered but cannot replace
this primary split after results are known.

Development training readiness requires non-empty train and development sets,
at least five development conflict components and non-zero retained quotient
rank in each split. This is an engineering Gate, not a population inference
Gate.

## Quotient implementation

Each batch contains complete panels. For panel `p`, compute the canonical
float64 additive projection and loss:

```text
L_p = ||P_perp (y - q)||^2 / retained_rank_p
```

then average panel losses. Tests must establish additive-null annihilation,
disconnected-graph rank, numerical orthogonality, permutation invariance,
non-zero autograd and CPU/GPU agreement. A row minibatch that splits a panel is
invalid.

## Linear witness

Fit one ridge-regularized linear response. The ridge value is selected using
train components only from the fixed grid `1e-6, 1e-4, 1e-2, 1`; development
labels are opened once after selection. Report:

- quotient RMSE and explained fraction versus zero interaction;
- correct pair versus ligand-only, protein-only, foreign-ligand and
  deranged-protein controls on identical masks;
- coefficient norm, effective rank, gradient participation and replay hashes;
- GPU preprocessing throughput separately from the linear solve.

No threshold may be changed to produce a PASS. A positive development witness
only authorizes an independently registered confirmation source; it does not
admit the statistic to `model/` or biological `z`.

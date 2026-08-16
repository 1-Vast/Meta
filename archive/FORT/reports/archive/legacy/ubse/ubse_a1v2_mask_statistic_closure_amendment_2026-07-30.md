# UBSE-A1-v2 mask and statistic closure amendment

Date: 2026-07-30

Status: binding static amendment; no coordinate, CCD, parse, event, student,
affinity, confirmation, or sealed-data authorization

Parent contract:
`manifests/ubse_a1v2_stage1_extractor_contract.v1.json`

Parent SHA-256:
`e2e66f1143c93ce886af86a976abe4dcb28f5677ab332451998b9d301e638c92`

Decision:
`FREEZE_A1V2_MASK_AND_CONDITIONAL_TETRAD_CLOSURE_KEEP_DATA_LOCKED`

## 1. Reason for amendment

The parent contract defines an event-blind structural mask using a deposited
ligand within 6.5 A. That is valid for measuring cross-deposition event
reliability, but it is unavailable for a strict dual-cold student. Passing
that structural mask to Stage 3 would reveal the true binding neighborhood.

The parent also leaves the primary statistic as "log-score or deviance" and
uses a lexical independent cycle basis for weighting. Those choices are not
unique enough for a preregistered test. The basis is suitable for reporting
the dimension of the legal cycle space, but an arbitrary basis is not a
biological sampling unit and should not define a nonlinear likelihood.

This amendment closes those ambiguities before any event value is read.

## 2. Separate masks

Three masks have different roles.

### 2.1 Reliability mask

`M_rel(c,i,g,k)` is the parent contract's deposited-complex mask:

- the residue is resolved in model 1 of the exact receptor auth chain;
- a residue heavy atom is within 6.5 A of the selected ligand for the general
  channels, or within the diagnostic 8.5 A water envelope;
- the residue, FG orbit, ligand role, and event channel are chemically and
  atomically evaluable.

`M_rel` is allowed only for A1-R extractor/repeatability and manual-chemistry
measurement. It is never a Stage-3 input, normalization support, pocket
proposal, or evaluation mask.

### 2.2 Deployment mask

The no-P0A primary uses:

```text
M_deploy(t,l,i,g,k) = 1
```

if and only if:

- `i` is a real, unpadded position in the frozen target sequence;
- `g` is a frozen FG automorphism orbit in the ligand 2D graph;
- the FG ligand role is compatible with event channel `k`.

`M_deploy` is constructed without a PDB identifier, deposited coordinate,
resolved-residue flag, ligand pose, true pocket, template, P0A score, event
value, affinity, or confirmation outcome. The primary mask includes all
sequence positions; optional monomer or pose features may change scores but
may not remove cells or change normalization support.

### 2.3 Supervision observation mask

For A1-S and the later one-shot A1-C evaluation:

```text
M_obs(c,i,g,k) =
  M_deploy(t,l,i,g,k)
  and exact sequence-to-auth residue mapping exists
  and required deposited atoms are resolved
  and FG mapping is unique
  and the extractor completed.
```

There is no 6.5 A prefilter in `M_obs`. A value is zero only when
`M_obs=1` and no PLIP event maps to the exact cell. Missing residues, missing
required atoms, ambiguous mappings, or extractor failures remain null.
Sparse positive output alone still cannot define negatives: the complete
`M_obs` grid must be enumerated first.

Student predictions and `Pi0/Pi1` are computed on `M_deploy`. Supervised
cell or tetrad loss is evaluated only where `M_obs` is true. `M_rel` must not
be loaded by a Stage-3 data loader.

## 3. Unique conditional tetrad statistic

For an elementary legal 2-by-2 tetrad `q=(i,j,g,h,k)`, define:

```text
Omega_q =
  s[i,g,k] + s[j,h,k]
  - s[i,h,k] - s[j,g,k].
```

After extraction, an informative fixed-margin outcome has direction:

```text
y_q = +1 for [[1,0],[0,1]]
y_q = -1 for [[0,1],[1,0]].
```

Conditioning on one event in each row and column gives:

```text
P(y_q | margins, s) = sigmoid(y_q * Omega_q).
```

The pair-conditioned rank-one null is exactly `Omega_q=0`, hence probability
`0.5`. The per-tetrad primary contribution is the held-out deviance gain:

```text
d_q = 2 * (log sigmoid(y_q * Omega_q) - log 0.5).
```

The reported primary statistic is `d_q`; "log-score or deviance" is no
longer an analyst choice. Directional accuracy against `0.55` remains an
auxiliary materiality diagnostic and is not the conditional null.

The score is generated in five dependency-component-disjoint folds. It is
averaged equally over all informative elementary tetrads within each
`target x channel`, then equally across frozen dependency components.
Overlapping tetrads are not counted as independent observations.
Uncertainty and power use the frozen component resampling scheme.

## 4. Rank is not likelihood weighting

The exact rational rank and lexical greedy basis remain the unique topology
report:

```text
candidate legal cycle space -> exact rank -> lexical basis certificate
```

Training and the primary conditional statistic use all informative elementary
2-by-2 tetrads, with the target/channel normalization above. They do not use
only a lexical basis. General even circuits may be registered later as a
secondary sensitivity analysis, but cannot replace a missing 2-by-2 primary
or select members, channels, thresholds, or ranks.

## 5. Destruction invariants

Every wrong-target, wrong-ligand, FG, residue, position, or event destruction
must preserve, bit for bit:

- the outcome and observation mask being evaluated;
- `M_deploy`;
- the event-channel and FG-role multiset;
- real-real row and column margins;
- dustbin mass, score, mask, and real-real normalization;
- target/channel/component evaluation weights;
- the set of elementary tetrads entering the statistic.

Only the preregistered representation or assignment is changed. A
wrong-entity control is eligible only when a deterministic type/shape-matched
bijection preserves these objects. An ineligible match is null and cannot be
made legal by changing the mask, padding weights, or margins.

## 6. Stage consequences

- A1-R reliability may continue to use `M_rel` after a separately authorized
  A1-R acquisition and parse/extraction sequence.
- A1-S coordinate acquisition, event read, Stage-2 fitting, and Stage-3
  training remain blocked until an affinity-blind `M_deploy` ledger and its
  role-local hash certificate are frozen.
- Stage-3 code must prove that deleting `M_rel` from its inputs leaves outputs
  byte-identical.
- The no-P0A arm must not load P0A features or weights.
- A later pose arm may add features on `M_deploy`, but it cannot change the
  mask or use pose confidence derived from affinity.
- A1-C cannot select the mask, statistic, tetrad set, control matching,
  material threshold, or model setting.

## 7. Binding next action

The hash-first authorization packet may cover A1-R coordinate bodies and the
CCD snapshot only. It must bind this amendment but does not authorize its
network execution.

Before any A1-S body acquisition is proposed, freeze and test:

1. the deterministic `M_deploy` builder;
2. the sequence/auth-position mapping certificate;
3. complete `M_obs` enumeration with `0/1/null` semantics;
4. exact conditional tetrad deviance;
5. destruction-invariant checks.

All current execution authorizations remain false.


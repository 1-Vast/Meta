# Active experiment protocol

## Status

```text
S7_L2B_PHASE1_B5 ................. COMPLETE, DEVELOPMENT PASS 6/6
S7_L2B_PHASE2A_AUDIT ............. COMPLETE, AUDIT-ONLY, NOTHING TRAINED
S7_L2B_PHASE2B_PREREGISTRATION ... WRITTEN AND HASHED, NOT AUTHORIZED
NEW MODEL TRAINING ............... NOT AUTHORIZED
DOCUMENT-CLOSED CONFIRMATION ..... SEALED
TIME-FORWARD MONN CONFIRMATION ... INFEASIBLE BY DATA CONSTRUCTION
SOURCE AFFINITY / DAVIS / KIBA / z  FROZEN
GIT COMMIT ....................... NOT AUTHORIZED, NOTHING COMMITTED
```

## Phase 2A protocol as executed

Registered by `research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2A.md`
(SHA-256 `4e01401d…`) with amendments 01–03, each frozen before the phase it
governs. Every output artifact embeds those hashes.

| phase | executed | outcome |
|---|---|---|
| 0 contract and artifact audit | 26 artifacts hashed, 7 fail-closed checks | `PHASE2A_CONTRACT_PASS` |
| 1 data identifiability census | full corpus, label-blind power | `DATA_IDENTIFIABLE` |
| 2 teacher ligand-conditionality | replicate floor vs alternative ligand | teacher **is** ligand-conditioned |
| 3 marginal/coupling decomposition | weighted ALS, orthogonality `1.2e-9` | additive dominates |
| 4 matched attribution battery | 9 arms, degree-preserving rewiring | `BC = false` |
| 5 label semantics | type census, dense-distance comparator | not ambiguous |

The load-bearing Phase 0 check was `C3`: Phase 1's marginal decomposition
indexed the B5-family memmaps with the B4-family offset table. Both tables were
rebuilt independently and proved identical key-for-key, so the Phase 1 B5/BX5/BP5
marginal numbers were correctly aligned.

## Terminal verdict and mandated action

```text
LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING
  -> preregister one ligand-conditioned residue residual head
```

Precedence was applied in the registered order; rules 1–4 did not fire, and
`BC = false` with `TC = false` selects rule 7. The verdict was not chosen for
interest.

## Phase 2B registration boundary

`research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2B_RESIDUE_RESIDUAL.md`
(SHA-256 `ae6d1a01…`) freezes, before any Phase 2B number exists:

- the frozen prior `b_r(P)` and the single trainable residual `delta_r(P, L)`;
- rank `K <= 8` low-rank bilinear form over existing frozen states only;
- the mandatory projection away from constant, pocket-prior and ligand-only
  directions, with tolerance `1e-8`;
- the differential objective on symmetric-difference residues;
- the split, inference unit and differential-AUPRC metric;
- the replicate-oracle ceiling, so a modest absolute number is read correctly;
- Gates `D1`–`D5` with margins and lower bounds;
- a fail-closed module-participation audit.

It authorizes nothing by itself. No Phase 2B code exists.

## Decision rule for Phase 2B

```text
D1 fails
  -> the sequence-plus-2D residue-differential route is closed; do not
     substitute a larger model

D1 passes but D2 or D4 fails
  -> the head is fitting construct or annotation structure, not ligand
     chemistry; reject it

all gates pass and the module-participation audit passes
  -> a structural ligand-conditioned residue statistic is established, and
     nothing more
```

No branch authorizes affinity, DAVIS, KIBA, support adaptation, production `z`,
or any modification of the frozen probability-law operator.

## Historical record

Superseded and failed protocols are recorded in `history.md` and
`report/EXPERIMENTAL_EVIDENCE_LEDGER.md`. They are not part of the active
protocol.

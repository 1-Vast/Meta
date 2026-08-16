# UBSE-G1 centered residue-contact student decision

Date: 2026-07-29  
Decision: `STOP_UBSE_G1_NO_DEPLOYABLE_INTERACTION_RESIDUAL`

## Outcome

UBSE-G1 passed the frozen substrate and execution/firewall gates, but failed
every semantic and matched-destruction gate. A small frozen-ESM2 plus
two-dimensional-ligand student did not predict a held-domain,
ligand-conditioned residue-contact residual that was distinguishable from the
exact additive null or from either destruction control.

This is a Stage-1 semantic failure. It is not evidence that protein-ligand
interactions are absent, and it does not authorize ChEMBL affinity loading,
UBSE-G2, Stage-2 fitting, confirmation access, or sealed-test access.

## Frozen execution

- Source SHA-256:
  `7905e4edf88073f564baa4b2d4fb50d496432bc4e15e97cccbfa0766b1b0638d`.
- G0PB manifest SHA-256:
  `4fea01e332eb3c60e41d76d5062d33cc95b13bc2e96b01df226532f78fe1b371`.
- Only `target_key`, `sequence`, `pubmed`, `scaffold`, `conn`, and
  `binding_residues_reindexed` were loaded from the source.
- No affinity field or value, development/confirmation feature or label,
  structure archive, or sealed outcome was loaded.
- ESM2 inference and all six student fits ran on
  `NVIDIA GeForce RTX 4060 Laptop GPU` under torch `2.6.0+cu124`.
- Three seeds (`1729`, `1730`, `1731`), 30 fixed epochs, batch size eight,
  and identical within-seed initialization and panel order were used.
- Wall time was 263.172 seconds. A training-phase sample showed approximately
  7-8% SM utilization, 12 W, and 2.84 GB framebuffer use. The run was valid
  but underutilized the GPU; future preregistrations must use larger
  length-bucketed batches or a more vectorized interaction kernel and must
  record utilization and power telemetry prospectively.

## Substrate

| Role | Valid panels | Panels with contact contrast |
| --- | ---: | ---: |
| fit | 1,259 | 1,138 |
| validation | 64 | 57 |
| untouched audit | 88 | 81 |

One panel was rejected because two ligand records became duplicates under the
frozen non-isomeric canonicalization. All 88 audit identities remained
present. Audit homology components, scaffolds, and PubMed identifiers were
each unique, and all three resources had zero overlap with fit or validation.

The final audit ledger contains 81 contrast panels times three seeds times
four controls, or 972 complete rows. Every seed-panel has exactly one null,
cross, ligand-deranged, and protein-free-position evaluation.

## Held-audit results

Median across the three frozen seeds:

| Arm or control | Directional accuracy | Centered cosine | Assignment Recall@1 | Centered MSE |
| --- | ---: | ---: | ---: | ---: |
| cross | 0.5107 | 0.0689 | 0.5201 | 0.2452 |
| exact additive null | 0.4804 | -0.0409 | 0.4388 | 0.2499 |
| ligand deranged | 0.4868 | -0.0697 | 0.4169 | 0.2681 |
| protein-free-position | 0.5467 | 0.0871 | 0.5355 | 0.2438 |

The panel-over-seed directional deltas and frozen 2,000-replicate
panel-bootstrap intervals were:

| Contrast | Estimate | 95% interval |
| --- | ---: | ---: |
| cross minus exact additive null | 0.0219 | [-0.0513, 0.0952] |
| cross minus ligand derangement | 0.0263 | [-0.0705, 0.1224] |
| cross minus protein-free-position | -0.0209 | [-0.0838, 0.0425] |

The cross arm's seed-specific directional scores were 0.4907, 0.5366, and
0.5107. Its apparent training-fit advantage therefore did not transport:
the final training loss was approximately 0.47-0.48 for cross versus
1.06-1.08 for the null, while the independent audit metrics remained close
to chance and unstable across seeds.

## Frozen gates

| Gate | Result | Evidence |
| --- | --- | --- |
| S1 substrate | pass | 1,138/57/81 contrast panels and 88/88 audit identities |
| S2 nontrivial cross signal | fail | directional 0.5107 < 0.60; cosine 0.0689 < 0.10 |
| S3 exact-null increment | fail | delta 0.0219 < 0.05; lower bound -0.0513 |
| S4 ligand destruction | fail | delta 0.0263 < 0.05; lower bound -0.0705 |
| S5 protein-position destruction | fail | delta -0.0209; lower bound -0.0838 |
| S6 pair assignment | fail | 0.5201 < 0.60 and below strongest control 0.5355 |
| S7 execution/firewall | pass | finite, parameter-matched, CUDA, no forbidden access |

## Independent implementation review and correction

An independent agent verified:

- zero homology/scaffold/PubMed leakage among fit, validation, and audit;
- complete paired audit rows with no duplicate `(seed, panel, control)`;
- identical null/cross parameters, initialization, panel order, and update
  count;
- matched within-panel ligand cycling and mask-aware residue-position
  destruction; and
- the preregistered directional, centered-cosine, assignment, and tie rules.

The review found one non-outcome-changing discrepancy. The original S3-S5
point checks subtracted control-specific cross-seed medians, whereas the
preregistration requires first averaging each panel over seeds. The bootstrap
code already implemented the preregistered estimand. The gate code was
corrected to use the bootstrap estimate, and the existing ledger was
recomputed without retraining. All three gates still fail by wide margins, as
shown above.

The code also now binds reusable ESM caches to explicit CUDA-creation
provenance instead of allowing a future hard-coded CUDA pass. The accepted
run did not reuse a cache, reported `esm.cuda=true`, and is unaffected.

## Interpretation and next boundary

G0R showed that observed BioLiP contact labels contain ligand-specific
information. G1 now shows that this information was not recoverable on the
closed held domains from the chosen deployable inputs and student. In
particular, a position-free protein control outperformed the full cross arm.
The observed score can therefore be explained by panel/ligand-level
regularities without evidence that the model used residue position.

The stopped route must not be rescued by post-hoc threshold, seed, epoch,
width, loss-weight, or attention changes. A successor needs a genuinely new
deployable information object or training supervision and must again beat:

1. an exact additive target-plus-ligand null;
2. ligand destruction;
3. protein-position or structure destruction; and
4. strict homology, scaffold, PubMed/provenance, and outcome firewalls.

Until such a successor passes an affinity-blind semantic gate, the only
authorized work is further source/representation audit, prospective
measurement design, and preregistration. Affinity loading and Stage-2 training
remain locked.

## Artifacts

- Preregistration:
  `reports/active/ubse_g1_centered_contact_student_preregistration_2026-07-29.md`
- Result:
  `reports/active/ubse_g1_seed1729_1731.json`
- Audit ledger:
  `reports/active/ubse_g1_audit_ledger.parquet`
- Implementation:
  `research/ubse_g1.py`
- Tests:
  `tests/test_ubse_g1.py`


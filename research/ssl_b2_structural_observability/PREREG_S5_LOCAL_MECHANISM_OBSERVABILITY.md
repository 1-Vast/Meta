# S5 preregistration — local mechanism observability

Status: **AUTHORIZED FOR RESEARCH DESIGN AND LABEL-FREE STRUCTURAL TESTING ONLY**.
Real affinity values, DAVIS/recipient labels, production `z`, CSMO/Band changes,
and P2–P4 remain frozen.

## 1. Why S5 is necessary

S4 tested a linear map from a mean-pooled ESM representation and a 1024-bit
ECFP to six complex-level aggregate pseudo-teacher values.  That experiment is
a valid negative result for that particular aggregate probe.  It is not an
upper bound on all sequence+2D models.  The production P1B path contains
atom-local GINE states, residue-local ESM states, and atom-by-residue predicted
contact/distance tensors.  Those pair-local observables were absent from S4.

The historical S4 result is therefore retained, but its strongest admissible
verdict is:

```text
AGGREGATE_ESM_ECFP_PROBE_NOT_PROTEIN_SPECIFIC
PAIR_LOCAL_P1B_OBSERVABILITY_NOT_TESTED
```

`POSE_FREE_DEPLOYMENT_INPUTS_INSUFFICIENT` is not an identified terminal
conclusion until S5 evaluates the actual pair-local input contract.

## 2. Scientific question

S5 asks exactly one label-free question:

> Do frozen, deployment-observable atom-local ligand states, residue-local
> protein states, and P1B contact/distance predictions retain enough
> correct-protein information to reconstruct at least one prespecified local
> structure-mechanism channel under protein-homology and ligand-scaffold
> closure?

It does not ask whether a channel predicts affinity, and a structural PASS does
not admit any coordinate into `z`.

## 3. Data roles and sealing

1. The 14,906 historical holo records are development/training material.
2. The 1,118 S4 complexes are now development-exposed because their teacher
   outputs and test metrics were opened.  They may be used for diagnostics and
   power estimates, never as a new untouched confirmation set.
3. A new S5 confirmation block must be selected score-blind from RCSB entries
   absent from both exposure registries.  Selection, exclusions, PDB IDs,
   protein clusters, scaffolds, file hashes, and the split hash are frozen
   before any S5 prediction is scored.
4. Every split closes protein homology, exact protein sequence, PDB ID, exact
   ligand identity, Bemis–Murcko scaffold, and teacher/development exposure.
5. The inference and bootstrap units are closure components, never atom–residue
   pairs.

No affinity value is read in S5.

## 4. S5-A — deployment-contract and mapping audit

Before fitting any probe, construct deterministic mappings:

```text
CCD atom serial <-> canonical ligand-graph atom index
mmCIF polymer residue <-> canonical target sequence index
canonical sequence index -> P1B 128-slot index
```

The primary stratum contains one deployment target chain and one non-covalent
drug-like ligand.  Multi-chain interfaces, covalent ligands, ambiguous residue
mapping, modified residues without a declared parent, and ligand atom mappings
below 95% are excluded by rules frozen before outcomes are inspected.  The
homo-oligomer stratum is reported separately because a single target sequence
does not specify oligomer geometry.

Required audit outputs include atom mapping coverage, residue mapping coverage,
residues per occupied slot, mixed-residue-class slots, chain multiplicity,
discard reasons, and all hashes.  Failure is
`S5_DATA_OR_MAPPING_CONTRACT_FAIL_CLOSED`.

## 5. S5-B — pseudo-teacher fidelity and slot ceiling

The existing S2 values are deterministic structural pseudo-labels, not physical
free-energy components.  S5-B upgrades their chemistry inputs using the CCD
bond graph, formal charge, donor/acceptor flags, and aromatic ring membership.
Protein donor/acceptor, charge, and aromatic definitions are versioned.  Waters,
metals, protonation uncertainty, solvation, and entropy remain explicitly
unobserved.

The upgraded teacher must pass rotation, translation, atom permutation,
residue permutation, and repeatability tolerances of `1e-9`.  Exact-residue
pair contributions are then deterministically summed into P1B slots.  The
following ceilings are reported:

- exact teacher aggregate reconstructed from exact pair contributions;
- slot teacher aggregate reconstructed from slot-summed contributions;
- oracle slot-chemistry reconstruction using explicit slot composition;
- loss caused by replacing explicit slot composition with the averaged ESM
  slot state.

This stage separates an input/mapping failure from a model realization failure.

## 6. S5-C — zero-training observability ladder

All arms use identical rows, splits, task weights, and metrics.

| Arm | Inputs | Purpose |
|---|---|---|
| `C0` | train mean | null |
| `C1` | ligand-only local chemistry | ligand shortcut |
| `C2` | pooled ESM + ECFP | historical S4 continuity |
| `C3` | atom chemistry x explicit slot chemistry | pair chemistry without P1B geometry |
| `C4` | `C3` + frozen P1B contact and 5-bin distance probabilities | geometry retention |
| `C5` | `C4` + frozen atom-local and residue-local states | full deployment-observable P1B input |
| `CO` | exact holo pair coordinates and teacher chemistry | oracle ceiling, not deployable |
| `CR` | capacity-matched random features | realization/control |
| `CD` | `C5` with a score-blind `<40%` wrong protein from the same split | nuisance derangement control |
| `CG` | within-complex pair-geometry shuffle | geometry necessity control |

`CD` is never called a biological non-binder.  The map is frozen before scores,
uses different closure components, limits reuse, and records identities and
hashes.

Simple deterministic linear or low-rank probes are used here; no architecture
or hyperparameter search is allowed.  If `CO` fails, the teacher/mapping or
evaluation pipeline is defective.  If `CO` passes but `C5` fails, the tested
pose-free P1B observables are insufficient for that channel.  If `C5` passes,
lightweight GPU distillation is authorized.

## 7. S5-D — synthetic trainability control

Before real structural pseudo-label training, generate a known pair-local
teacher from the frozen `C5` inputs.  The proposed head must recover it on held
out closure components and must beat ligand-only, deranged-protein, and
pair-shuffle controls.  This detects broken optimization, normalization, or
aggregation without reading affinity.

Failure is `S5_HEAD_OR_OPTIMIZATION_NOT_IDENTIFIED`; it is not evidence that
biology is absent.

## 8. S5-E — lightweight local mechanism distillation

All upstream modules are frozen: ESM, ProteinEncoder, ligand GINE, P1B bridge,
contact head, and distance head.  The only trainable object is a small shared
pair map:

```text
u_ij = [atom-local state, residue-slot-local state,
        atom chemistry, explicit slot composition,
        p(contact_ij), p(distance_ij)]
e_ij = g_theta(u_ij) in R^6
E(P,L) = invariant masked aggregation of e_ij
```

`g_theta` is a 1–5M parameter MLP or an equally small separable low-rank map.
It is evaluated in chunks so a dense `[B,N,128,H]` tensor is not retained on
the RTX 4060.  Batch size is 2–4 with AMP and gradient accumulation.  Frontend
states are cached; no ESM or GINE update is allowed.

The primary loss is train-only robust-scaled Huber on pair-local channel
contributions, averaged within complex and then across closure components.  A
secondary aggregate Huber term may be used with its coefficient frozen before
validation.  No affinity loss, target ID, ligand ID, assay ID, wrong-protein
training loss, or dataset-specific constant is allowed.

## 9. Structural feasibility criterion

For a prespecified core channel, all of the following must hold on the sealed
S5 confirmation block using component bootstrap:

```text
R2(correct vs train mean) >= 0.02, 95% LCB > 0
R2(correct) - R2(ligand-only) >= 0.02, 95% LCB > 0
R2(correct) - R2(deranged) >= 0.02, 95% LCB > 0
R2(correct) - R2(pair-shuffle) >= 0.02, 95% LCB > 0
```

Pair-local reconstruction and complex-level aggregate reconstruction must both
have the same favorable direction.  Channels failing prevalence or mapping
preconditions are `NOT_EVALUABLE`, not zero-valued failures.

Permitted terminal verdicts are mutually exclusive:

```text
S5_DATA_OR_MAPPING_CONTRACT_FAIL_CLOSED
S5_TEACHER_OR_SLOT_CONTRACT_DEFECT
S5_HEAD_OR_OPTIMIZATION_NOT_IDENTIFIED
P1B_PAIR_LOCAL_MECHANISM_NOT_OBSERVABLE
P1B_PAIR_LOCAL_STRUCTURAL_MECHANISM_OBSERVED
```

Only the last verdict authorizes a separately registered real source-affinity
stage.  It does not authorize DAVIS, production integration, or biological `z`.

## 10. Downstream decision tree

```text
S5 structural PASS
  -> freeze named structural channels
  -> source-only ChEMBL/BindingDB closure-OOF affinity calibration
  -> require correct > ligand and correct > deranged by >=0.03 with LCB>0
  -> only then study a k<=5 rank-aware support section
  -> only then define compact, observable, bounded biological z coordinates
  -> only then test the unchanged K(B(z)F(z)) operator

S5 C5 FAIL but oracle CO PASS
  -> pose-free P1B inputs insufficient for the named channel
  -> separately register pose-aware S7

S5 oracle CO FAIL
  -> repair teacher/mapping/evaluation; do not add model capacity
```

The final support adapter must obey the frozen identifiability theory: the
adapted coefficient lies in the support-design row space, the identifiable rank
is at most `k`, query coverage and conditioning are reported, and the model
abstains when the query direction is not covered.  Structural reconstruction
and affinity ranking are not claimed to be covered by the frozen point-valued
regression theorem.

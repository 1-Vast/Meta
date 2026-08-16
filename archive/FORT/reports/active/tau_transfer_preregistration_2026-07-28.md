# Tau transfer: label-free A0/T0 preregistration

**Date:** 2026-07-28

**Scope:** conceptual and metadata-only feasibility. No affinity, development,
confirmation, or sealed-test labels may be read.

## Question

The transferable abstraction from *tau: Learning Touch-Augmented
Vision-Language-Action Models from Future Visual Supervision* is:

```text
current interaction state + a verifiable intervention
    -> change in interaction state
```

The audit asks whether the current DTA program has both:

1. a label-free bilateral intervention graph with canonical protein and ligand
   edits; and
2. an admissible frozen teacher that carries pair-specific interaction
   information rather than target-only, ligand-only, family, source, pose
   confidence, or training-membership shortcuts.

It does not treat time windows, future offsets, ligand pairs, or affinity
quartets as independent samples.

## Frozen stages

### `TAU-A0`: conceptual isomorphism

Pass only the narrow design principle: training-only, intervention-conditioned
delta supervision may improve representation efficiency. This pass provides no
biological, affinity, power, or strict dual-cold performance credit.

### `TAU-G0`: bilateral intervention topology

Read only TRAIN metadata columns:

```text
target, conn, endpoint, scaffold, assays, docs, accession, hcluster,
dual_cold_split
```

Protein transitions require explicit `base_protein`, `construct_sequence`, and
`directed_substitution` fields. Ligand transitions require an explicit,
canonical `directed_ligand_edit`. Dense same-target panels, shared scaffolds, or
target identifiers are not substitutes. If either side is unavailable,
bilateral factorial contrasts are not identifiable and the route stops.

### `TAU-T0`: teacher admissibility

Each teacher is tested separately. Admission requires all of:

- frozen version;
- available weights or deterministic engine;
- resolved training lineage or a verified no-training construction;
- frozen rights;
- inputs constructible without privileged test entities;
- fold-local construction;
- pretraining-overlap exclusion;
- pair-specific output;
- executable destruction controls.

The required negative controls are a separable target-plus-ligand embedding,
matched random pair teacher, and constant teacher. They cannot qualify as
semantic teachers.

## Candidate order

1. PBCNet2.0;
2. a frozen physical pose/contact teacher;
3. LEXOR-MC qualitative constraints;
4. generic separable embeddings;
5. random and constant controls.

## Downstream unlock

`TAU-S0` synthetic calibration is allowed only if both `TAU-G0` and `TAU-T0`
pass. It must reuse the existing R-MAON direct regular-null carrier; synthetic
success is engineering evidence only.

Real teacher association (`TAU-I0`) additionally requires OpenMut `I0` to pass.
Predictive integration (`TAU-M0`) requires `TAU-G0`, `TAU-T0`, `TAU-S0`,
`TAU-I0`, and coordinate gate `OMUT-C0`. The inference model, basis, capacity,
split, weights, and nested null must remain identical with and without the
training-only auxiliary loss.

This order applies only to a non-structural, independently measured teacher.
Any structure-, pose-, docking-, contact-, or complex-conditioned teacher
remains deferred until `OMUT-M1` passes and may enter only through `OMUT-R0`.
Tau terminology cannot reopen the earlier physical-structure route.

Failure is terminal for the current teacher. It may not be rescued by changing
the number of offsets or teacher views, loss weight, rank, backbone, capacity,
epochs, seeds, or coordinate concatenation.

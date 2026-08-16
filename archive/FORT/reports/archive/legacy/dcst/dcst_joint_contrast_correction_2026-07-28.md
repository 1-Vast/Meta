# DCST within-target joint-contrast correction

Date: 2026-07-28  
Status: frozen before joint-contrast training  
Supersedes only: the absolute-only joint privileged loss and mechanism gate in
`dcst_joint_privileged_correction_2026-07-28.md`

## Trigger

The 500-step absolute-joint pilot produced two certified source affinity
bands while the matched no-privileged model produced none. The intermediate
joint mechanism probe, however, showed:

```text
true joint CE          3.235
wrong-target joint CE  3.304  (+0.069)
wrong-ligand joint CE  3.250  (+0.015)
```

The joint target was still dominated by a target's common pocket-interaction
distribution. Absolute cross-entropy did not sufficiently identify which
ligand changes which segment–interaction cells.

## Corrected loss

Within each target-balanced Stage-1 episode, covered rows are converted to:

```text
Delta P_i = P_i - mean_ligands(P).
```

The predicted softmax distributions are centered identically. Their flattened
cosine alignment is maximized. The frozen privileged objective becomes:

```text
0.25 * absolute joint cross-entropy
+ 1.00 * (1 - centered joint cosine alignment).
```

The weaker absolute term retains a calibrated distribution; the main term
removes the target-common pocket prior and supervises ligand-specific
structural changes. Episodes with fewer than two covered ligands contribute
only the absolute term.

## Corrected source mechanism gate

On firewalled source-development targets with at least two covered ligands:

- true centered alignment must be positive;
- true alignment must exceed wrong-target alignment by more than 0.05;
- true alignment must exceed within-target wrong-ligand alignment by more
  than 0.05;
- true absolute joint CE must remain below uniform `log(64)`.

No downstream threshold, arm, rank, data row, seed, step count or confirmation
rule changes.


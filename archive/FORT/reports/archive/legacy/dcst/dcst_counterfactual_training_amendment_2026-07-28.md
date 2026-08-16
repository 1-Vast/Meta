# DCST training-time counterfactual amendment

Date: 2026-07-28  
Status: frozen before the counterfactual pilot  
Supersedes only: the positive-only centered privileged objective in
`dcst_joint_contrast_correction_2026-07-28.md`

## Trigger

The 500-step centered-joint pilot improved ligand destruction but did not pass
the source mechanism gate:

```text
true centered alignment          +0.0407
wrong-target centered alignment  +0.0657
wrong-ligand centered alignment  -0.0663
```

Positive-only alignment could therefore learn a globally useful ligand
mapping without requiring the correct target. This is a Stage-1 source
identifiability failure, not a ChEMBL development-score adjustment.

## Frozen correction

For every Stage-1 target-balanced episode with privileged rows, define:

```text
L_true = centered joint-distribution alignment loss for the correct pair
L_T-   = the same loss after replacing the target representation
L_L-   = the same loss after cyclically deranging ligands within the episode
```

Add the margin loss:

```text
L_cf = relu(0.10 + L_true - L_T-) + relu(0.10 + L_true - L_L-)
```

The complete privileged term is frozen as:

```text
0.25 * absolute joint cross-entropy
+ 1.00 * L_true
+ 0.50 * L_cf
```

The wrong-target representation is produced by the existing deterministic
target derangement. The wrong-ligand control changes only the ligand input;
the target and privileged target distribution remain fixed. No ChEMBL label
enters this loss.

## Interpretation and unchanged controls

The amendment makes correct target–ligand pairing a training requirement,
rather than merely a post-training probe. It does not change the source
firewall, spectral bands, source certificate threshold, Stage-2 arms, seeds,
step counts, downstream success thresholds, or confirmation policy.

`DCST-NoPriv` remains the load-bearing-information control. The source
mechanism probe must still pass both destructive margins before the privileged
route can be interpreted as having learned joint information.

The runner's `--stage1-only` execution mode is the required engineering gate
for longer source optimization. It may read label-blind downstream entity
metadata to rebuild the cross-source firewall, but it must not load or score
downstream affinity labels.

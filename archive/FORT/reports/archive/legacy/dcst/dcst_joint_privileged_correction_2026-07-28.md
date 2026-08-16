# DCST joint privileged-target correction

Date: 2026-07-28  
Status: frozen before joint-target training  
Supersedes only: the two separate marginal privileged losses in
`dcst_privileged_contact_amendment_2026-07-28.md`

## Trigger

The 500-step marginal-target mechanism probe showed:

| Auxiliary target | True CE | Wrong-target CE | Wrong-ligand CE | Uniform CE |
| --- | ---: | ---: | ---: | ---: |
| contact segment | 1.808 | 1.987 | 1.808 | 2.079 |
| interaction type | 1.366 | 1.366 | 1.386 | 2.079 |

The student learned both tasks, but separably: contact location was mostly a
target pocket prior and interaction type was mostly ligand chemistry.
Successful marginal prediction therefore did not establish pair-specific
information and produced zero certified affinity bands.

## Corrected privileged target

Each PLINDER interaction record binds a residue identifier to one interaction
type. The source projector maps that residue's sequence index into one of
eight ordered protein segments and constructs one normalized `8 x 8`
distribution:

```text
P(segment, interaction_type | target, ligand).
```

The student applies an explicit ligand-conditioned head to every protein
segment and emits 64 joint logits. The two marginal losses are removed. On
covered source-train rows the only privileged loss is:

```text
1.00 * cross_entropy(joint_8x8_distribution).
```

The affinity residual loss and all Stage-2 settings are unchanged.

## Frozen mechanism gate

On covered, firewalled source-development targets with at least two ligands:

1. true-pair joint cross-entropy must be below uniform `log(64)`;
2. wrong-target joint cross-entropy must exceed true by more than 0.05;
3. within-target wrong-ligand cross-entropy must exceed true by more than
   0.05.

This gate is evaluated before the affinity-transfer gate. `DCST-NoPriv`
remains the matched architecture control.


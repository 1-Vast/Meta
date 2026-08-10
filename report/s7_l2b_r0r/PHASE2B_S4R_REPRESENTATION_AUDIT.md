# Phase 2B S4R-A ligand representation audit

Terminal verdict: `GRAPH_LIGAND_REPRESENTATION_AVAILABLE_AND_INFORMATIVE`

Label-blind. Residue label reads: 0. Affinity value reads: 0. Trainable parameters introduced: 0.

## Baseline mean-pooled 41-D reference

| statistic | value |
|---|---:|
| effective rank of the embedding | 5.336 |
| numerical rank of the embedding | 33 |
| effective rank of pair differences | 6.183 |
| distinct graphs sharing an identical vector | 687 in 307 groups |
| scaffold-distinct pairs with cosine > 0.999 | 0.00241 |
| R2 of difference norm on heavy-atom-count delta | 0.8517 |

A1 therefore requires a candidate difference effective rank of at least `18.549`.

## Candidate grid

| radius | d | eff.rank(dg) | INC | RET | coverage | INC_perp | W params | admissible |
|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| 1 | 128 | 20.93 | 0.3552 | 0.0780 | 0.99997 | 0.3544 | 163840 | True |
| 1 | 256 | 28.44 | 0.3910 | 0.0514 | 1.00000 | 0.3901 | 327680 | True |
| 1 | 512 | 35.08 | 0.4143 | 0.0262 | 1.00000 | 0.4134 | 655360 | True |
| 2 | 128 | 21.78 | 0.3799 | 0.0802 | 0.99997 | 0.3791 | 163840 | True |
| 2 | 256 | 30.85 | 0.4236 | 0.0510 | 1.00000 | 0.4228 | 327680 | True |
| 2 | 512 | 41.09 | 0.4588 | 0.0283 | 1.00000 | 0.4579 | 655360 | True |

`INC` is the fraction of candidate pair-difference energy that no
linear function of the baseline pair difference can express. `RET` is
the converse loss. `INC_perp` additionally removes heavy-atom-count and
log-count differences. All are heldout-A, label-blind.

## Selection

Capacity-parsimony rule selects radius `1`, vocabulary `128`, `W` parameters `163840` against the baseline's `52480`.

Frozen vocabulary SHA-256 `a200a4b986af1850fdb1d244f2e002c9b5ae707a114d8a3635053edb215ed877`.

## Boundary

This audit measures a ligand representation only. It admits no
statistic to `z`, opens no residue label, and does not modify
`A(F,z) = K(B(z)F(z))`.

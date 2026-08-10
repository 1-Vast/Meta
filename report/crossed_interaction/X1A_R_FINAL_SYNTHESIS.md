# X1A-R direct-DD dependence: final synthesis

## Terminal verdict

```text
X1A_R_DEPENDENCE_PRECONDITION_FAILED
```

The repaired audit used the same statistic planned for the interaction test,
`Z=(DD/2)^2-v_D,U`, on exact-assay-aligned, label-blind selected rectangles.
It fit no target or ligand nuisance model. The inference unit was the frozen
dependency component.

| endpoint | rectangles | components | rho UCB95 | registered maximum | effective n | required |
|---|---:|---:|---:|---:|---:|---:|
| Ki | 827 | 36 | 0.120406 | 0.0915 | 200.43 | 245 |
| Kd | 590 | 12 | 0.101078 | 0.0164 | 61.05 | 245 |

Both endpoints fail both the dependence threshold and the conservative
effective-unit requirement. Cluster dominance itself passed (`0.0387` Ki,
`0.2119` Kd). The result is therefore an information/dependence stop for this
public ChEMBL crossed design, not a model-training or GPU failure.

## Noise and access

The exact-assay pooled measurement-noise upper SD is `0.6928` p-units for Ki
and `0.6967` for Kd. Only 13 Ki rectangles and one Kd rectangle have all four
cells directly replicate-supported; the remainder use the conservative pooled
noise contract.

The run opened 5,986 preselected ChEMBL37 pChEMBL rows. It read no BindingDB,
DAVIS, KIBA, PDBbind or recipient value, trained no model, and used no GPU.

## Consequence

The conditionally preregistered X1B interaction-existence test was not run.
X2, privileged 3D supervision, support sectioning, biological `z`, and the
frozen probability-law operator remain untouched. A future route requires a
new, independently governed crossed source with substantially more independent
components; thresholds may not be relaxed to rescue this corpus.

Repository regression after consolidation: **203 passed** in `drug`.

Machine evidence: `report/crossed_interaction/x1ar_direct_dd/gate.json` and
`dd_rows.jsonl`.

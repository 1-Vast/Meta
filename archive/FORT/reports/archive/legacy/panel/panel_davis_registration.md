# Registration of the independent confirmation panel (Davis 2011 pKd)

Registered 2026-07-25, before any label of this panel was read. This closes the single item that
blocked `ORRC_PD0_READY` in `reports/active/orrc_eb_blueprint_v2.md`: section F requires an
evaluation source the ORRC route has never seen, because ORRC-EB was designed after the Gate PB and
PC development outcomes on the Metz panel were observed.

## Why this panel

Candidate scan over every ChEMBL document in the frozen local extract, restricted to
sequence-resolved targets and to targets outside the main registry's sealed confirmation components,
using the same `(20 ligands, 5 targets)` dense-core rule as the development panel:

| document | endpoint | targets | components | ligands | fill | median ligands/target |
|---|---|---:|---:|---:|---:|---:|
| `CHEMBL1908390` Davis | pKd | 116 | **102** | 68 | 43.5% | 28 |
| `CHEMBL3991601` | pKd | 64 | 57 | 127 | 24.8% | 26 |
| `CHEMBL1909046` | pKi | 24 | 24 | 116 | 48.7% | 57 |
| `CHEMBL5442175` | pKi | 17 | 17 | 53 | 54.3% | 31 |
| `CHEMBL1150977` | pKd | 14 | 14 | 29 | 74.9% | 22 |

`CHEMBL1908390` is the only candidate whose independent-component count (102) is comparable to the
development panel's (101); the statistical unit is the target homology component, so component count
is what determines resolvable effect size. `CHEMBL3991601` is retained as a secondary candidate but
carries a data-quality flag: it contains a `pKd` of `-4.603`, i.e. a reported Kd of ~40 mM, which
needs investigation before use.

## Data-quality audit, performed before the split was chosen

| document | min | max | top-3 value share | median within-target duplicate-value fraction |
|---|---:|---:|---:|---:|
| `CHEMBL1908390` Davis pKd | 5.004 | 10.796 | 5.4% | **0.053** |
| `CHEMBL1201862` Metz pKi (development panel) | 4.000 | 11.100 | 15.0% | **0.817** |

Davis shows no censoring floor: the minimum is `5.004` rather than a spike at exactly `5.0`, and no
single value exceeds 2% of records. This matters because a Kd panel truncated at 10 uM would place a
large tied block at the floor and corrupt a within-target ranking metric.

The same audit records a material property of the *development* panel that was not previously
measured: Metz values are rounded to 0.1 pK, so 81.7% of a target's values duplicate another value
of the same target. Gate PA/PB/PC remain valid — Spearman uses mid-ranks and every arm faced the
same ties — but the development panel's ranking resolution per measurement is much coarser than
Davis's, and that is now on the record.

## The registered substrate

Built by `tools/panel_registry.py`, the generalized builder, which
`tests/test_panel_registry_equivalence.py` proves reproduces the frozen Metz registry cell-for-cell.
The confirmation substrate is therefore produced by the same audited code path as the development
substrate.

`dataset/public/chembl_37/processed/panel_davis/`, registry sha256
`f15daa5478f63a648a07d52d76aee588e4dc6d7275444fc50d204a774a3499fe`:

* endpoint pKd only, document `CHEMBL1908390` only. **Never pooled with pKi.** Panel-specific
  intercept, ligand basis, noise scale and endpoint scale, per blueprint section F;
* train 2,069 cells / 116 targets / 102 homology components / 42 anchor ligands / 42 scaffolds /
  42.5% fill / median 17 anchor ligands per target;
* development 1,360 cells / 116 targets / 102 components / 26 query ligands / median **12** query
  ligands per target; maximum anchor Tanimoto `0.6923`, so the query chemistry is far more distant
  from the anchors than the development panel's `0.9091`;
* ligand axis exactly disjoint: 0 shared scaffolds and 0 shared parent connectivities between train
  and development;
* 21 targets whose homology component appears in the main registry's sealed confirmation split were
  dropped, so this panel cannot contaminate that sealed set;
* no target-axis sealed slice. This panel is **single use**: the whole panel is the sealed asset.

## Independence audit

* different publication, different endpoint (pKd vs pKi), different assay platform;
* **0** Davis cells are main-registry confirmation cells; **0** are main-registry development cells;
* no gate of the ORRC route — PA, PB, PC, the power records, or the void PD0 run — has ever read a
  Davis label;
* overlap with the development panel is reported, not hidden: 74 shared targets, 31 shared parent
  connectivities, 33 shared scaffolds, and 7 Davis *development* ligands that appear among Metz
  *train* ligands. This does not compromise the design, because the confirmation gate is
  **self-contained**: ORRC is fitted on Davis train cells and evaluated on Davis development cells,
  so dual-coldness is enforced within Davis (held-out homology component, scaffold- and
  connectivity-disjoint query ligands). A secondary cross-panel transfer diagnostic restricted to
  the 42 Davis targets that are absent from the Metz panel may be reported, clearly labelled as a
  diagnostic and never as the gate.

## Single-use policy

`consumed=false` in the manifest. It is flipped by exactly one authorized run: the first scoring of
a **target-conditioned arm** on Davis development cells. Until then:

* an **arm-blind** power audit (ligand-only base at four seeds under the leave-component-out
  protocol) is permitted and does **not** consume the panel, because it scores no target-conditioned
  arm and therefore reveals nothing about ORRC. This is the same "power audit first" doctrine that
  task.md already fixes for every gate in this program;
* the go/no-go decision after that audit is a pre-registered function of one number only: the paired
  MDE80. The base's absolute macro Spearman is reported but is explicitly **not** a decision input,
  so knowing it cannot influence whether the confirmation gate is attempted.

## Pre-registered threshold rule for the confirmation gate

Gate PD-C's minimum detectable effect must budget for **arm heterogeneity**, not same-arm retraining
noise. The development panel measured both on the identical protocol: retraining spread `0.0948` per
component and observed `I0 - B0` arm spread `0.2201` per component, an inflation factor of `2.32`;
the `CFRI - B0` and `R0 - B0` contrasts give `1.82` and `1.66`. The registered rule, fixed here
before the Davis audit runs, uses the largest of the three:

```text
MDE80_PD-C = max( 0.03 , 2.32 x MDE80_retraining(Davis) )
```

If `MDE80_PD-C` exceeds `0.0614` — the arm-heterogeneity MDE80 already recorded for the development
panel in `reports/active/panel_power_pd1.json` — then Davis is underpowered relative to the effect
sizes this program cares about, that fact is recorded, and no confirmation gate is run on it. A
failed power check is a result, not a reason to relax the rule.

## What this registration authorizes

The arm-blind Davis power audit, and nothing else. It does not authorize a PD-C run, a PD-M or PD-X
stage, multi-seed runs, Hierarchical MoT, long training, a signed prior mean, or any access to the
main registry's confirmation split.

---

## Arm-blind power audit outcome, 2026-07-25 — Davis is UNDERPOWERED, no confirmation gate runs

`research/panel_power.py --panel davis`, four seeds of the ligand-only base under the identical
leave-component-out protocol, no target-conditioned arm scored, panel `consumed` still `false`.
Report: `reports/active/panel_power_davis.json`.

| quantity | development panel (Metz) | confirmation candidate (Davis) |
|---|---:|---:|
| independent homology components | 101 | 102 |
| median query ligands per target | 43 | 12 |
| B0 macro Spearman by seed | 0.2786 / 0.3022 / 0.2938 / 0.2913 | 0.1438 / 0.1620 / 0.1151 / 0.1846 |
| per-component retraining spread | 0.0948 | **0.2314** |
| retraining MDE80 | 0.0181 | **0.0688** |

Applying the rule frozen above, `MDE80_PD-C = max(0.03, 2.32 x 0.0688)`:

```text
MDE80_PD-C = 0.1596   >>   0.0614   (the development panel's arm-heterogeneity reference)
```

Davis therefore cannot resolve the effect sizes at stake, and by the rule registered before the
audit ran, **no confirmation gate is run on it**. The panel remains sealed and unconsumed; only its
base-arm retraining spread has been read.

The cause is not the component count, which is essentially identical to the development panel (102
versus 101). It is the query depth: 12 query ligands per target instead of 43 inflates the
per-component Spearman noise by roughly `sqrt(43/12) = 1.9`, which is the whole of the observed
`0.2314 / 0.0948 = 2.44` ratio. Component count and query depth are separate resources and only the
second is scarce here.

### What an admissible confirmation source must satisfy

Quantified so the next attempt is not a guess. Combining the four independent candidates by inverse
variance — Davis 102 components at 12 query ligands, `CHEMBL3991601` 57 at 26, `CHEMBL1909046` 24 at
57, `CHEMBL5442175` 17 at 31 — gives about 384 Davis-equivalent components and
`MDE80_PD-C ~ 0.082`, still above `0.0614`, so even pooling every remaining panel in the local
extract does not reach the required resolution. An admissible source needs **both**

* at least ~100 independent homology components, **and**
* at least ~40 query ligands per target after the scaffold-disjoint split,

which is a second Metz-scale panel. The frozen local ChEMBL-37 extract does not contain one: the
document scan found exactly one panel of that density, and it is already the development substrate.

This is a measurement-design conclusion, not a model result. It does not weaken the ORRC-EB
mathematics, and it does not license running PD-C on an underpowered source or relaxing the rule.

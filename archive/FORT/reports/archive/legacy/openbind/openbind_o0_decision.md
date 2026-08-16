# OpenBind O0 local-SAR information audit

Date: 2026-07-26  
Decision: **k=4 signal passes after a strict chemical firewall, but oracle/target evidence is
insufficient; stop before model training**

## Scope

OpenBind EV-A71/CVA16 2A is a dense **single-target** campaign. It can test whether a few
same-series labels contain useful local SAR information, but it cannot test protein conditioning,
protein shuffling, or unseen-protein transfer.

The official release reports 601 compounds with affinity measurements. After applying the
official benchmark filters, requiring a structure-linked pKD, canonicalizing SMILES, and
collapsing duplicate structures, the usable audit set contains 488 unique compounds.

## Protocol

- Morgan fingerprints, radius 2 and 2,048 bits.
- Butina chemical-series clustering at center similarity 0.50.
- Independent statistical unit: chemical series, not molecule or random draw.
- 64 paired support/query draws per eligible series at k = 4, 8, and 16.
- B0: ExtraTrees ligand-only baseline trained outside the held-out series.
- B0 firewall: canonical identity disjoint, Bemis–Murcko scaffold disjoint, and maximum Morgan
  Tanimoto to every held-series compound below 0.50.
- k-shot: B0 plus Tanimoto-weighted residuals from exactly k labeled support compounds.
- Permutation control: the same support compounds with their labels permuted.
- Dense oracle: every other compound in the held-out series is labeled; neighbor count and
  correction strength are selected by nested local leave-one-out Spearman. This is intentionally
  more information than a deployable k-shot method receives.

## Result

| Shot | Independent series | k-shot − B0 Spearman | 95% LCB | Dense oracle − B0 | Oracle LCB |
|---:|---:|---:|---:|---:|---:|
| 4 | 13 | +0.0984 | +0.0138 | +0.1089 | −0.0704 |
| 8 | 8 | +0.0970 | −0.0123 | +0.1377 | −0.0362 |
| 16 | 4 | +0.1338 | −0.0686 | +0.1370 | +0.0153 |

At k=4, the exact-support arm passes its registered performance and mechanism checks: gain exceeds
+0.03, its component LCB is positive, and support beats label permutation by +0.1110
[+0.0285,+0.2038]. RMSE improves by 0.1827 pKD [+0.0142,+0.3973]. k=8 has similar point estimates
but its ranking and permutation lower bounds cross zero. k=16 has only four independent series and
is not adequately powered.

The chemical firewall materially weakens B0, as it should: it removes 8–168 additional
same-scaffold/high-similarity neighbors per held series, leaves 161–471 external training
compounds, produces exactly zero scaffold overlaps, and caps observed train-to-held Tanimoto at
0.4935. The earlier non-firewalled O0 result is invalid for strict chemical OOD and is superseded.

The full-local-label oracle has positive point headroom at every k, but its LCB crosses zero at the
only adequately populated settings (k=4/8). The k=16 oracle interval is positive but rests on four
series and fails the minimum-eight-series gate. This heterogeneity prevents treating one
single-target campaign as evidence for a general support encoder.

## Decision

Verdict: `OPENBIND_K4_SIGNAL_ORACLE_UNCERTAIN_STOP`. The k=4 result is genuine evidence that
chemically proximal same-target support can repair a chemically hard ligand baseline. It is not
evidence for target-conditioned dual-cold transfer: there is only one protein, protein
shuffle/random controls are impossible, and oracle variability across 13 series remains large.
Do not build a support encoder, posterior module or protein model from this target alone. Reopening
requires the same protocol on independent targets or families; changing the current model after
reading this development result is not allowed.

## Reproduction

```powershell
D:\anaconda\envs\drug\python.exe research\openbind_o0.py
D:\anaconda\envs\drug\python.exe -m pytest tests\test_openbind_o0.py -q
```

Machine-readable result: `reports/active/openbind_o0.json`.

Sources:

- OpenBind Zenodo release: https://zenodo.org/records/20026661
- Official benchmark repository:
  https://github.com/OpenBind-Consortium/EV-A71_2A_benchmark

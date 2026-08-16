# DCST privileged-contact Stage-1 amendment

Date: 2026-07-28  
Status: frozen before any privileged-contact Stage-1 or downstream result  
Parent records:
`dcst_two_stage_preregistration_2026-07-28.md` and
`dcst_two_stage_certificate_correction_2026-07-28.md`

## Trigger and non-rescue boundary

The corrected 500-step source-only pilot found zero positive spectral-band
utilities (`-0.02`, `-0.11`, `-0.12`, `-0.15`). Full unfiltered transfer
reduced the short ChEMBL pilot from B0 Spearman 0.1108 to 0.0886. The shared
pooled/segment ESM plus Morgan coordinate therefore did not learn a
target-specific source reordering signal.

This amendment does not relax a gate, tune a certificate threshold, or inspect
a 4,000-step downstream result. It changes the scientific information supplied
to Stage 1: structure-derived PLINDER contacts are used as training-only
privileged supervision and must prove load-bearing against an otherwise
identical no-privileged arm.

## Label-blind privileged projection

The local PLINDER annotation table contains per-ligand:

- interacting protein residues with sequence indices;
- hydrogen bond, hydrophobic, water bridge, salt bridge, pi-stack,
  pi-cation, halogen-bond and metal-complex interaction records;
- pocket and crystallographic validation metadata.

No coordinate cache and no raw affinity column are needed. After the frozen
cross-source firewall, exact `(system_id, ligand parent connectivity)` joining
provides privileged labels for:

- 2,344 / 5,577 allowed Stage-1 train-plus-development rows (42.03%);
- 2,124 / 5,130 Stage-1 train rows;
- 220 / 447 Stage-1 development rows;
- 1,700 covered rows whose source entry passes the PLINDER validation flag.

The projection reads no affinity field. Eight equal-length sequence bins map
the terminal sequence index in each interacting-residue record to the existing
eight ESM segment tokens. Interaction types form a separate normalized
eight-class distribution.

## Corrected Stage-1 student

The ligand embedding now queries the eight ordered protein segment tokens.
Its attention logits predict the observed contact-segment distribution, and a
pair head predicts the observed interaction-type distribution. The
ligand-conditioned protein state and ligand state still enter the same direct
bilinear affinity residual with exact null `Theta=0`.

On covered Stage-1 train rows only:

```text
L_stage1
= L_registered_affinity_residual
+ 0.50 * cross_entropy(contact_distribution)
+ 0.50 * cross_entropy(interaction_type_distribution).
```

Rows without privileged labels retain the registered affinity-residual loss.
No privileged value is used on ChEMBL. At Stage 2 and inference the model
requires only the already available ESM features and ligand features.

## Required new control

`DCST-NoPriv` uses the identical ligand-conditioned architecture,
initialization, data, source base, optimizer steps, certificate rule and
Stage-2 residual, but omits both privileged losses. DCST must have grouped
bootstrap LCB95 greater than zero against `DCST-NoPriv` in addition to the
existing scratch and naive-fine-tune comparisons.

All other Stage-2 arms, rank-eight four-band certificate, MDE 0.0586, source
firewall, seeds, steps, downstream destruction controls, RMSE guard and
confirmation isolation remain unchanged.


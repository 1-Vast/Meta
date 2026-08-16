# SAFSA-G0 decision: family-selectivity anchor is not identifiable

Date: 2026-07-27  
Preregistration SHA-256:
`6ad649657e80973c03d7fc63e63df856a2cfdb37517097b79aa343bab443956e`  
Result SHA-256:
`59d31a613b8e43f6c1820371a072249b1ec630a7d85a41818e70dfa2879a6ee9`

## Verdict

**`SAFSA_G0_FAMILY_SELECTIVITY_NOT_IDENTIFIABLE_STOP`.**

Do not recover structures or sequences for this route and do not train the proposed
Selection-Aware Family-Selectivity Anchor. The point estimate is suggestive in only part of the
taxonomy, but the intended general cross-family anchor is not supported.

## Frozen-gate result

The public-evidence-clean matrix contains 67 evaluable genes from all six target classes and 67,999
drug--gene cells (92.2% tested negatives), so sample eligibility and the inactive-retaining
selection controls pass.

| frozen contrast | family-macro AUPRC difference | hierarchical 95% CI | gate |
|---|---:|---:|---|
| own family - family cold | +0.0466 | [-0.0630, +0.1496] | FAIL |
| own family - wrong family | +0.0423 | [-0.0395, +0.1234] | FAIL |
| own family - selection coverage | +0.2055 | [+0.1228, +0.2972] | PASS |
| own family - global promiscuity (sentinel) | +0.0047 | [-0.0889, +0.1070] | non-gating warning |

Removing all DrugCentral mechanism-of-action cells does not rescue the claim:
own-family minus family-cold is +0.0441 [-0.0645, +0.1467], failing the frozen +0.02/LCB>0
sensitivity gate.

## Failure anatomy

The effect is not a weak but coherent six-family signal:

- GPCR: +0.1519 own-minus-cold;
- nuclear receptor: +0.1568;
- kinase: +0.2069, but only one evaluable gene;
- enzyme: -0.0163;
- ion channel: -0.0833;
- transporter: -0.1362.

Thus equal-family uncertainty correctly crosses zero. A gene-weighted result would be dominated by
the 42 GPCR genes and would answer a different, leakage-prone question.

The selection sentinel passes strongly, so the positive point estimate is not explained only by
which same-family assays were run. However, own-family AUPRC (0.3560) is virtually identical to
global ligand-promiscuity AUPRC (0.3513). The available signal is therefore mostly broad ligand
activity plus taxonomy-specific pockets of coherence, not a generally transferable family
compatibility anchor.

## Scientific consequence

The SPD substrate remains valuable as a systematic, inactive-retaining pretraining source, but it
does not authorize this particular representation target. Training SAFSA now would optimize a
family label whose incremental identifying information is unresolved and whose sign reverses in
half of the broad target classes. Model capacity or a larger GPU cannot repair that statistical
failure.

Local memory is explicitly not treated as a scientific route constraint: once a mechanism passes
its identifying gate, full-scale training may move to another machine. Here the stop is caused by
non-identifiability, not compute.

Firewall status: SPD was used only as a pretraining-anchor audit; current-run FORT development and
confirmation labels were not read; historical confirmation remains quarantined;
`sealed_test_consumed=false`.


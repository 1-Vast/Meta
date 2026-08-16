# FACTOR F0 decision

Date: 2026-07-26  
Protocol verdict: `FACTOR_F0_FAIL_STOP_BEFORE_F1`  
Failure type: `F0_DISCRETE_LIGAND_CARRIER_COVERAGE_FAIL`

The corrected authoritative result is
`reports/active/factor_f0.json` (SHA-256
`A7F7C5239972B5D2E5F918FCB90DC5FB708488FD1A6A5DD21FBD70CA95E1EA45`).
No numerical activity column, quarantined ChEMBL confirmation row or sealed test was read.

## Gate result

| Quantity | Result | Frozen gate | Decision |
|---|---:|---:|---|
| valid 85-position pocket mapping | 98.73% | >=95% | pass |
| ligand carrier mapping | 100% | >=99% | pass |
| anchor weighted median / q10 | 1.000 / 0.956 | >=0.95 / >=0.80 | pass |
| carrier weighted median / q10 | 0.801 / 0.604 | >=0.90 / >=0.70 | **fail** |
| non-novel-primitive observations | 79.89% | >=70% | pass |
| document graph largest component | 30/30 environments, all 3 sources | >=80%, all sources | pass |
| equal-source maximum weight | 33.33% | <=40% | pass |
| grouped design MDE80 | 0.0226 | <=0.03 | pass |

The source-balanced graph is connected and the active-site anchor side is well covered, but the
*exact discrete-token carrier ontology* is not. Of 54,323 eligible observations, 10,922 (20.11%) are
classified as `novel_primitive` under exact count-Morgan/BRICS identity. This result does not yet
distinguish genuine chemical extrapolation from semantic aliasing in an over-discrete tokenization.
The cross-source pair-overlap median is high (0.922) only among the deterministic sampled local pairs;
it cannot compensate for the large low-coverage tail of whole-ligand carrier sets.

## Adapter correction

The first implementation accidentally used the structure-only KLIFS registry and mapped 69.58% of
observations. That non-authoritative output is preserved as
`factor_f0_uncorrected_structural_adapter.json` (SHA-256
`6A0267C1B35F4F1A1767AC89EB10747A3AB099E1D33E9A93FC7F60A7BEC8D32F`).
The authoritative rerun changed only the identity adapter to the already-audited, label-blind
KLIFS kinase-information registry. No source, token, threshold or endpoint rule changed. Mapping rose
to 98.73%, while the pre-existing carrier failure remained.

## Consequence

The original exact-token F1 is not authorized. The next permissible experiment is a separately
preregistered `F0-C0` continuous-carrier audit that keeps the 0.90/0.70 coverage gates, calibrates
distance scales on train-only pseudo-OOD scaffolds, and includes chemistry-broken decoys,
reconstruction proxies and effective-rank checks. Lowering the gate, tuning Morgan/BRICS until they
pass, using Papyrus-ChEMBL31 after its frozen exclusion, or training FACTOR anyway is forbidden.

# UBSE-G0R binding-residue reliability decision

Date: 2026-07-29  
Decision:
`REQUEST_UBSE_G0C_SMALL_COORDINATE_RELIABILITY_PREREGISTRATION`

## Outcome

The existing, firewalled BioLiP2 binding-residue substrate passes the frozen
cross-publication reliability and ligand-specificity gate:

| Metric | Result |
| --- | ---: |
| Eligible closed observations | 66,660 |
| Repeated exact-sequence / ligand-connectivity units | 1,028 |
| Exact targets / ligands | 800 / 648 |
| PubMed IDs / PDB entries | 1,679 / 1,936 |
| Units / targets with legal same-target wrong-ligand controls | 753 / 525 |
| Median correct cross-PubMed contact Jaccard | 0.7500 |
| Median size/scaffold-matched wrong-ligand Jaccard | 0.5000 |
| Median unit ligand-specific margin | 0.10417 |
| Target-median margin, bootstrap 95% interval | [0.1250, 0.2000] |
| Optimistic retrieval queries / targets | 1,076 / 220 |
| Median ligand candidates | 11 |
| Contact-label Recall@1 / random expectation | 0.5604 / 0.0846 |
| Contact-label MRR / random expectation | 0.6698 / 0.2391 |

All six frozen gates pass. The run parsed only the seven allowed structural
metadata columns. It decoded no affinity field or value, downloaded no
structure, and touched no development, confirmation, or sealed outcome.

## What this resolves

The result falsifies the strongest negative explanation for this particular
substrate: the reindexed binding-residue label is not merely a fixed
target/pocket signature. Across independent PubMed records and different PDB
entries, the same non-isomeric ligand connectivity is measurably more
reproducible than a hard wrong connectivity on the same sequence-exact
target. Stereo, charge, tautomer, construct, and biological-assembly identity
remain unresolved and are mandatory in the next gate.

This is the first current-source result that simultaneously has:

- a real ligand-conditioned observation;
- independent cross-publication repeats;
- hundreds of exact targets;
- a direct same-target wrong-ligand control; and
- no affinity dependence.

It therefore survives the failures that stopped the PLINDER absolute moment,
DTIOD tangent, exact-pair RDIB delta, recurring-edit RDIB, and PD-MVR bridge.
Those routes required replicated ligand-pair differences or affinity bridges;
UBSE-G0R instead tests repeatability of an absolute but ligand-specific
interaction state.

## What this does not resolve

The retrieval diagnostic consumes the observed contact labels. It is an
optimistic information ceiling, not a sequence-plus-2D prediction result.
The gate does not show that:

- fine-grained atom/residue interaction events are reliable;
- a sequence/SMILES model can infer the held structure state;
- the signal survives homology, scaffold, chemical-neighbour, and provenance
  closure simultaneously;
- a learned state adds information beyond protein-only, ligand-only, family,
  scaffold, or memorized-complex baselines; or
- the state predicts affinity in strict dual-cold evaluation.

The margin also clears the frozen floor narrowly (`0.10417` versus `0.10`).
No threshold may be relaxed, and the next experiment must preserve
target-cluster uncertainty and hard wrong-ligand controls.

## Binding continuation

The only authorized continuation is to preregister a small coordinate-level
reliability pilot on a deterministic subset of these repeated units. It must:

1. obtain at least two different PDB structures from different PubMed records
   per exact target-ligand unit;
2. define mmCIF entity, assembly, chain, ligand-instance, altloc, metal,
   cofactor, water, missing-residue, mutation, and covalent-ligand rules before
   download;
3. extract fixed atom/residue interaction-event tokens without affinity;
4. compare correct repeats with same-target hard wrong ligands and
   structure/identity-only controls;
5. keep a PubMed-, homology-, and scaffold-disjoint audit partition; and
6. stop before teacher/student training if event reliability or independent
   packing fails.

No 20k-50k pretraining run, affinity readout, or predictive claim is
authorized by G0R.

## Artifacts

- Preregistration:
  `reports/active/ubse_g0r_binding_residue_reliability_preregistration_2026-07-29.md`
- Result:
  `reports/active/ubse_g0r_seed1729.json`
- Unit ledger:
  `dataset/public/biolip2/processed/ubse_g0r_units.parquet`
- Implementation:
  `research/ubse_g0r.py`
- Tests:
  `tests/test_ubse_g0r.py`

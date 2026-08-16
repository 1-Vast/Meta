# FIRE-DTA Literature and Data Audit

## Decision

`FIRE_DTA_GATE0_CONDITIONAL_PASS`

The proposed information flow is mathematically coherent and can be tested without changing the
validated B1 posterior. It is not yet evidence that bound-free channels transfer to the ChEMBL-37
cold-target estimand. IPBind is already a close published standard, and BERP retains standard Bayesian
linear-update mathematics. Both candidates remain claims to falsify, not validated innovations.

## Primary literature

| Work | Year | Relevant mechanism | Consequence for FIRE-DTA |
| --- | ---: | --- | --- |
| [IPBind](https://doi.org/10.1109/OJEMB.2026.3667030) | 2026 | Shared geometric encoding of complex, protein and ligand followed by bound-minus-unbound atomic-potential subtraction | Mandatory strongest BFEO control; bound-free subtraction alone is not novel |
| [PLINDER](https://doi.org/10.1101/2024.07.17.603955) | 2024 | More than 400k annotated protein-ligand systems, similarity-aware splits, and linked apo/predicted structures | Supports a leakage-safe structural subset and genuine apo/holo checks; not every system has an apo structure or affinity |
| [MISATO](https://doi.org/10.1038/s43588-024-00627-2) | 2024 | Roughly 20k curated complexes with QM ligand properties and explicit-water MD trajectories | Supports state stability and strain probes; trajectories are bound-complex dynamics, not binding free-energy labels |
| [BigBind](https://doi.org/10.1021/acs.jcim.3c01211) | 2024 | ChEMBL activities mapped to representative 3D pockets | Feasible structural warm-up, but ChEMBL33 is an ancestor of the current labels and must be quarantined as additional-data evidence |
| [PDBbind CleanSplit](https://doi.org/10.1038/s42256-025-01124-5) | 2025 | Removes structural and ligand similarities that inflate affinity benchmarks | Requires protein, pocket, ligand, pose and label-neighborhood separation, not only target IDs |
| [Cofolding memorisation audit](https://doi.org/10.1101/2025.02.03.636309) | 2025 | Tests whether protein-ligand cofolding has moved beyond pose memorisation | Boltz/other generated states require template and similarity controls |
| [Boltz-2](https://doi.org/10.1101/2025.06.14.659707) | 2025 | Joint structure and affinity model | May generate candidate states only; pretrained affinity outputs and label-guided pose selection are forbidden |
| [UMA](https://arxiv.org/abs/2506.23971) | 2025 | Universal interatomic potential trained across large atomic datasets | Optional ligand/local-energy teacher only; it is not a validated protein-ligand binding free-energy oracle |

The supplied citation `[3]`, arXiv `2504.16261`, resolves to the IPBind preprint. It is not a separate
"LFM" source. The proposal therefore has one, not two, independent literature precedents for its
bound-free claim.

## Data reachability and scale

- PLINDER `2024-06/v2` is publicly reachable. The official GCS inventory contains a 0.948 GiB index,
  0.557 GiB of split assets and 128.818 GiB of system archives. The official tutorial warns that the
  complete resource approaches 1 TB and supports lazy loading.
- MISATO Zenodo record `7711953` is CC BY 4.0. `MD.hdf5` is 132,841,014,019 bytes (123.72 GiB;
  MD5 `9bc6446922cd80e0f2f3f69349bf88ed`); `QM.hdf5` is 343,064,967 bytes.
- BigBind V1.5 is reachable as a 19,254,374,400-byte archive (17.93 GiB). It uses ChEMBL33 labels,
  so those labels cannot be treated as independent validation of a ChEMBL37 model.
- The workstation has 264.2 GiB free on `D:`. A selected PLINDER subset is feasible. PLINDER plus
  full MISATO plus working caches is not a responsible local plan.

## Identifiability boundary

`complex - protein - ligand` is an exact representation contrast, but it is not automatically a
thermodynamic free energy. If the "free" protein and ligand are merely extracted from the bound pose,
the model has no reorganization or solvation observation. FIRE-DTA may use the terms `state channel`
and `bound-minus-free contrast`; it may claim free-energy semantics only after all of the following:

1. linked experimental apo/holo systems show a gain over bound-coordinate extraction;
2. native/decoy and branch-swap controls destroy the signal;
3. an IPBind-style single-state model does not match the multi-state gain;
4. the representation transfers to target-disjoint affinity without ancestor-label leakage;
5. pocket, ligand and pose similarity deletion does not erase the result.

## Novelty boundary

BFEO's candidate novelty is the ensemble operator, confidence-aware statistical marginalization and
anti-bypass data flow. Bound-minus-free subtraction and invariant atomistic encoding are standards.

BERP exactly preserves B1's Cholesky posterior with `phi=[1,z_phys]`. That is scientifically desirable
but not sufficient for an algorithmic novelty claim. BERP becomes load-bearing only if the same BFEO
architecture with ligand-latent B1 loses, physical-channel permutation destroys the gain, and the
effect improves harmful adaptation rather than only the target intercept.


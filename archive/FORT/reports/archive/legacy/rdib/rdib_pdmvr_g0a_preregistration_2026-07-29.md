# RDIB / PD-MVR G0-A preregistration

Date frozen: 2026-07-29  
Status: frozen before candidate enumeration  
Purpose: optimistic, affinity-blind screening of whether BioLiP2 contains
cross-publication replicated same-target ligand-pair differences.

## Fixed source and firewall

- Source:
  `dataset/public/biolip2/processed/closed_registry.parquet`
- Frozen SHA-256:
  `7905E4EDF88073F564BAA4B2D4FB50D496432BC4E15E97CCCBFA0766B1B0638D`
- Permitted BioLiP2 columns:
  `target_key`, `sequence`, `accession`, `pdb_id`, `pubmed`, `conn`,
  and `scaffold`.
- Permitted ChEMBL TRAIN metadata:
  `target`, `accession`, `conn`, `scaffold`, `hcluster`, and
  `dual_cold_split`.
- Prohibited fields:
  BioLiP2 `affinity_presence`, every numeric affinity or relation field,
  confirmation/development outcomes, and sealed-test outcomes.
- This audit is CPU-only set enumeration and graph bookkeeping. No
  affinity model or CUDA operator is authorized.

## Primary unit

For a sequence-exact target `t`, unordered exact ligand parents `(a,b)`,
and a non-empty PubMed identifier `p`, a source-specific difference is
eligible only when the closed registry contains both `(t,a,p)` and
`(t,b,p)`. A replicated candidate block `(t,a,b)` requires at least two
distinct PubMed identifiers, each independently containing both ligands.

Different PubMed identifiers are only an optimistic article-level
separation. They do not establish independent institution, author,
construct, cofactor, resolution, or experimental lineage. Passing this
screen therefore cannot authorize model training.

## Frozen metrics and gates

The audit will report:

1. exact target-ligand observations repeated across at least two PubMeds;
2. exact target-ligand-pair difference blocks repeated across at least
   two PubMeds;
3. exact-target, ligand, scaffold, and PubMed breadth;
4. largest target and PubMed evidence shares;
5. an optimistic conflict-free packing ceiling;
6. exact accession and ligand overlap with ChEMBL TRAIN metadata.

The optimistic packing ceiling is

`min(B, 3*T, floor(L/2), S, floor(P/2))`,

where `B` is replicated block count, `T` exact-target breadth, `L`
ligand breadth, `S` scaffold-token breadth, and `P` PubMed breadth.
Two ligands within one difference may share a scaffold, so scaffold
breadth contributes `S`, not `floor(S/2)`; scaffold tokens may not be
reused across packed blocks.
For this optimistic ceiling only, a ligand with an empty Murcko scaffold
receives its own ligand-specific scaffold token; missing/acyclic
scaffolds therefore cannot create a false early stop.
The target capacity of three is deliberately optimistic and homology,
family, chemical-neighbour, author, institution, construct, and assay
closure can only reduce it.

G0-A passes this screening only if all are true:

- at least 200 replicated exact ligand-pair difference blocks;
- at least 80 sequence-exact targets;
- optimistic conflict-free packing ceiling at least 200;
- no single PubMed supports more than 5% of candidate blocks;
- the firewall records zero protected values loaded.

If any condition fails, stop RDIB structural extraction and PD-MVR
structure-block packing without downloading coordinates. If all pass,
the next required audit is provenance/homology closure plus exact
conflict-constrained packing; PLIP extraction and contact ICC remain
locked until that audit passes.

The recurring-directed-edit route is not silently substituted for the
primary unit. It requires a separately frozen stereochemistry,
tautomer, charge, atom-mapping, and edit-equivalence contract.

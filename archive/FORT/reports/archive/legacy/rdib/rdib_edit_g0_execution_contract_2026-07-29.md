# RDIB-Edit-G0 execution contract

Date frozen: 2026-07-29  
Status: frozen reproduction contract after an explicitly non-binding scout  
Role: affinity-blind optimistic upper bound for the recurring-directed-edit
secondary route.

The multi-agent scout preceded this formal artifact. Its first global edit
count used the wrong independence unit and was withdrawn. This contract
freezes the corrected biological unit and reproduces it in reviewed code.
All thresholds are inherited from the supplied proposals and the existing
MMP-X contract; none were selected from the scout result.

## Source and firewall

- Input:
  `dataset/public/biolip2/processed/closed_registry.parquet`
- SHA-256:
  `7905e4edf88073f564baa4b2d4fb50d496432bc4e15e97cccbfa0766b1b0638d`
- Permitted columns:
  `target_key`, `pubmed`, `pdb_id`, and `conn`.
- Prohibited:
  `affinity_presence`, every affinity value/relation, all
  development/confirmation features and outcomes, and sealed outcomes.
- CPU is required because the work is RDKit fragmentation and exact set
  enumeration, not tensor algebra.

## Frozen optimistic edit definition

Reuse the existing `MMP-X F0` implementation without changing its edit
vocabulary:

- RDKit `rdMMPA.FragmentMol`, `maxCuts=1`;
- common core at least 10 heavy atoms;
- exchanged substituent at most 5 heavy atoms, excluding the dummy atom;
- lexical orientation of canonical dummy-labelled substituents;
- a source-specific molecule pair is eligible only when all shared cores
  imply exactly one edit token;
- zero-token and multi-token molecule pairs are excluded.

The input connectivity was generated with `isomericSmiles=False`.
Consequently this is only a non-isomeric optimistic upper bound. It cannot
pass a final identity gate and cannot justify coordinates or training. If
the upper bound survives, a later gate must reconstruct exact isomeric
identity from `ligand_ccd`, original SMILES, InChI and full InChIKey.

## Correct biological replication unit

For each `(target_key, PubMed)` panel, enumerate eligible ligand pairs and
their unique edit. Collapse to source blocks
`(target_key, edit, PubMed)`.

A replicated biological unit is exactly:

`(target_key, edit)` supported by at least two distinct PubMed identifiers.

An edit occurring once on many different targets is not replicated
interaction supervision. Different targets cannot substitute for a repeated
measurement of the same target-conditioned contact change.

Each replicated unit is classified as:

- exact-pair replicated, if one exact ligand pair occurs under at least two
  supporting PubMeds; or
- edit-only salvage, if replication is created only by different molecule
  pairs sharing the edit token.

## Inherited gates

Report the 88-unit optimistic mechanism floor, the 200-block/80-target
structural floor, and the approximately 423-unit predictive scale
separately.

Stop immediately, without strict stereo reconstruction or coordinates, if
any of these monotone upper bounds fails:

- fewer than 88 replicated `(target_key, edit)` units;
- fewer than 80 sequence-exact targets;
- optimistic resource ceiling below 88.

The resource ceiling is

`min(U, 3*T, floor(P/2), floor(L/2))`,

where `U` is replicated units, `T` exact targets, `P` distinct PubMeds and
`L` distinct ligands. It deliberately ignores homology, families,
scaffolds, chemical neighbours, constructs, author/institution lineages,
stereochemistry and tautomer ambiguity. Every strict closure can only lower
the attainable count.

Passing this reproduction would only request strict identity and conflict
closure. It would not authorize affinity access, coordinates, PLIP, contact
ICC, operator ranks, or model training.


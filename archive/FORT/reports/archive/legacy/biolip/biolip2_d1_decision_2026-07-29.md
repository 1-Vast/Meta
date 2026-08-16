# BioLiP2-D1 decision

Date: 2026-07-29  
Decision: `STOP_BIOLIP2_D1_EXACT_COMPLEX_TOPOLOGY_INADEQUATE`

## Positive evidence

After the protected development/confirmation firewall, the safe registry has:

- 562,794 exact-complex rows;
- 103,969 sequence-exact targets, 25,830 ligands, and 86,675 PDB entries;
- 37,063 PubMed IDs;
- 12,113 targets with different ligands in distinct PDB entries, yielding
  426,410 same-target ligand pairs;
- exact accession support for 420/559 ChEMBL-TRAIN targets (`75.1342%`);
- Morgan support >=0.40 for `49.815%` of the fixed 20,000-ligand TRAIN sample.

This is the first audited public structural layer in the program that passes
both target and ligand support by a wide margin.

## Failure

The target-ligand 2-core has 62,917 edges, 23,652 targets, 7,342 ligands, and
164 connected components, but its largest component contains 61,327 edges
(`97.4729%`). The ordinary ligand graph is therefore joined by generic or
highly promiscuous ligands and fails the frozen 50% limit.

All retained rows have PubMed provenance, but the largest PubMed ID
contributes `2.00535%`, just above the frozen 2% cap. D1 is stopped; neither
threshold is relaxed.

## Firewall

Affinity columns 14-17 were never decoded, parsed, copied, hashed, logged, or
written. Only their nonempty presence bits were retained. No structure
archive, development/confirmation feature or label, or sealed-test value was
loaded.

## Consequence

BioLiP2 remains admissible as an exact-complex/contact identity layer, but D1
does not unlock the directed-edit audit. A separately preregistered closure
audit may remove generic ligands and collapse repeated
`(PubMed, exact target, ligand)` structures before remeasuring independent
topology.

Authoritative machine result:
`reports/active/biolip2_d1_seed1729.json`.


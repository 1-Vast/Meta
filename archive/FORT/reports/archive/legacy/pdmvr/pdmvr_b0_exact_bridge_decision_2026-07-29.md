# PD-MVR-B0 exact bridge decision

Date: 2026-07-29  
Decision: `STOP_PDMVR_B0_EXACT_BRIDGE_UPPER_BOUND_INADEQUATE`

## Result

R18 had reported 348 ChEMBL TRAIN identity links. B0 tested the stronger
claim needed by PD-MVR: an allowed BindingDB row and a BioLiP2 complex had
to agree on exact PDB ID, canonical ligand connectivity, receptor sequence,
and accession.

The 1,267 firewalled BindingDB rows expanded to 24,777 listed PDB links, but
the four-way join retained only:

- 23 exact bridge candidates;
- 5 sequence-exact targets;
- 6 ligands and 6 PDB entries;
- 22 BindingDB DOIs, 16 non-empty institutions, and 5 BioLiP PubMed IDs;
- 2 known ChEMBL TRAIN homology components among the 19 mapped candidates.

The optimistic one-per-target/ligand/PDB/DOI/institution/PubMed ceiling was
only 5, versus the frozen minimum of 40. Broad-family mapping and bridge
common-rank computation were unnecessary because the upper bound had
already failed.

The local ChEMBL document map covered none of the five BioLiP PubMed IDs, so
zero same-article DOI matches were locally verified. This missing coverage
was treated as unverified rather than as a DOI mismatch; it is not the
reason for stopping. Even granting every unverified article and every
unmapped family optimistically cannot raise the target/provenance ceiling
above 5.

## Consequence

The 348 R18 links are identity overlap, not 348 exact dual-modal anchors.
PD-MVR bridge anchoring, coordinate download, contact-latent extraction,
stable common-rank calculation, affinity decoding, and model training are
stopped.

The source parser decoded only safe identity/provenance fields and Ki/Kd
presence booleans. Numeric affinity bytes remained undecoded; no
development/confirmation feature or outcome and no sealed-test outcome was
loaded.

Authoritative machine result:
`reports/active/pdmvr_b0_exact_bridge.json`.


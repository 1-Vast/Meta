# Q0-B alias ledger (first-hand sources only)

Every alias below is per-record evidence. None of it is generalized to
other proteins; a coordinate transform is applied only when its own
evidence chain is recorded here.

## BRAF V599E -> V600E (historical numbering)
- Reported: Duong-Ly Table S1/S2 construct BRAF(V599E); Davies et al.
  2002 Nature 417:949-954 (doi:10.1038/nature00766) reported V599E.
- Canonical: UniProt P15056 SV=4 residue 600 is V (variant annotation
  VAR_018629 p.Val600Glu); residue 599 is T.
- Sequence evidence (computed this stage, BRAF_HISTORICAL_EVIDENCE.json):
  the 1992 reference M95712.1 CDS is 2298 nt (765 aa); the current
  NM_004333.4 CDS is 2301 nt (766 aa) and its translation equals P15056.
  The historical CDS lacks exactly 3 nt (one codon) in the 5-prime
  region (new-sequence insertions at nt 88/90/97), shifting every
  downstream coordinate by +1: historical V599 == canonical V600.
- Applied transform: reported_position 599 -> canonical_position 600,
  kind=historical_alias, old residue V verified at 600 before use.
- NOT generalized: no other protein is given a +1 shift by analogy.

## PDGFRalpha (Duong-Ly)
- S1 lists GenBank NP_006197 (human PDGFRA) and Protein Accession
  Q9DE49. Q9DE49 is Danio rerio pdgfra (UniProt header PGFRA_DANRE).
- Correct human canonical: UniProt P16234 (PGFRA_HUMAN, 1089 aa), also
  resolved by KLIFS kinase lookup PDGFRa/Human -> P16234.
- D842V: S1 clone "Cytoplasmic (668-1210)" exceeds the canonical length
  (1089) -> QUARANTINED (excluded_construct_unresolved); the old residue
  D842 is verified on P16234 but the construct range cannot be mapped
  without silent repair.
- T674I: S1 Mutation column typo "T6741I"; construct name column reads
  T674I; T674 verified on P16234 -> admitted (notation fix recorded).
- V561D: V561 verified on P16234, construct 550-1089 -> admitted.

## KLIFS pocket numbering
- 85 aligned pocket positions; gatekeeper = index 45 verified on the
  known gatekeeper mutants (see Q0B_KLIFS_CENSUS.json).
- Sources: https://klifs.net (numbering help page, accessed 2026-08-18);
  Kooistra et al. 2016, Nucleic Acids Res 44:D365-D371,
  doi:10.1093/nar/gkv1058.

## Duong-Ly notation quirks
- S2 row labels use P38alpha/MAPK14 and TIE2/TEK (parent aliases).
- S1 Mutation column empty for EGFR deletion constructs; the construct
  name is authoritative there (d746-750, d747-749/A750P,
  d747-752/P753S, d752-759, d746-750/T790M).
- FLT3(ITD): "internal tandem duplication aa591-601" -> insertion class.


# Phase 2A preregistration — computational amendment 03

Parent: `research/s7_l2b_r0r/PREREG_S7_L2B_PHASE2A.md`
(SHA-256 `4e01401d0468e3241bd05cde93b2a517919209d770f8f75ab471e42228f1b12e`).
Prior amendments: 01, 02.

Written 2026-08-10, **after** Phases 3–4 and **before** any Phase 5
label-semantics metric. No parent or prior amendment is edited.

## C1 — the dense-distance comparator exists and must therefore be built

Section 9 of the parent says the label-semantics comparator is to be recorded as
`UNRESOLVED` only *if it does not exist*. It does exist: 2,068 of the 14,447
MONN PDB entries already have local mmCIF coordinates under
`dataset/raw/open_structures/` and `dataset/raw/ssl_b2_independent/`, acquired
under earlier governed stages. Declaring the question unresolved when a
comparator is available would be the fabrication the parent forbids, in the
opposite direction.

Phase 5 therefore builds a local dense-distance residue teacher and compares it
with the PLIP-derived binary residue mask. The construction is fixed here.

## C2 — residue index mapping, fail-closed

MONN indexes residues into the UniProt sequence; mmCIF indexes them by
`label_seq_id`. No SIFTS file is present, so the correspondence is established
by an exhaustive integer-offset scan:

- for each polymer chain, collect `(label_seq_id, one_letter_code)` over
  tabulated amino-acid residues;
- for every integer offset `o` in the admissible range, compute the identity
  rate of `seq[label_seq_id - 1 + o]` against the observed code over positions
  that fall inside the sequence;
- require at least **50** matched observed residues and an identity rate of at
  least **0.95**;
- among qualifying chains, select the one with the largest number of matched
  residues; ties are broken by chain name ascending.

A complex with no qualifying chain is **excluded and counted**. The exclusion
count and the identity-rate distribution are reported so the mapping yield is
visible rather than assumed.

## C3 — ligand copy selection is label-blind

A structure may contain several copies of the ligand CCD. The copy used is the
one with the greatest number of protein heavy atoms of the mapped chain within
4.5 A of any of its heavy atoms. This rule never consults the PLIP mask, so it
cannot be tuned toward agreement. Ties are broken by the copy's first atom
serial ascending.

## C4 — the geometric teacher and the reported contrasts

The geometric residue mask at threshold `t` is

```text
R_geom(t) = { r : min over ligand heavy atoms of dist(protein heavy atom of r,
                  ligand heavy atom) <= t }
```

for `t` in **{4.0, 4.5, 5.0} A**, all three reported. Hydrogens are excluded on
both sides. Waters, ions and other heteroatoms are excluded from the protein
side.

Reported against the PLIP mask `R_plip`: Jaccard, precision
`|R_plip ∩ R_geom| / |R_plip|`, recall `|R_plip ∩ R_geom| / |R_geom|`, and the
size ratio — each summarised by component-macro, with a component bootstrap.

The threshold sweep is the registered "dependence on distance threshold" test.
It is descriptive. It cannot by itself demonstrate label ambiguity, because a
distance teacher and an interaction-type teacher are different biological
objects: PLIP requires geometric *and* chemical criteria, so `R_plip` being a
strict subset of `R_geom` is the expected, correct behaviour and is **not**
evidence of missing positives.

## C5 — what would count as demonstrated ambiguity

Unchanged from the parent, section 9: at least 20% of positive edges indirect
**and** the `T1` conclusion reversing when they are removed; or a second frozen
teacher disagreeing on at least 20% of edges. For C4 the disagreement direction
that would matter is `R_plip` containing residues that are **not** within 5.0 A
of the ligand at all — a physically impossible direct contact — since that
would indicate a mapping or label defect rather than a criterion difference.
That quantity is reported explicitly as `plip_positives_beyond_5A`.

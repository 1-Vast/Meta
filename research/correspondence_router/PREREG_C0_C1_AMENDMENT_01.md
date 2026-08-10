# Amendment 01 — C0/C1

Parent: `PREREG_C0_C1_CORRESPONDENCE_INFORMATION_AUDIT.md`, SHA-256
`007f8439609078649cf7751b588716492f59c93bc27e2dad997b11afd7172c1e`,
committed `f844679`.

Written 2026-08-10, after the parent's own fail-closed mapping-equivalence check
rejected rule `M4`, and **before any C0 census or C1 statistic has been read
under a corrected mapping**. The C1 execution that had been started under the
rejected mapping was stopped, its partial output discarded and never inspected.

## 1. What the registered check caught

Parent section 4.2 asserted, as a matter of provenance, that the P1B corpus was
built with `record["sequence"]` equal to the mmCIF `_entity_poly` canonical
sequence, and therefore that

```text
M4  sequence_index(r) = int(label_seq_id) - 1
```

is an exact equivalence. The parent required that claim to be verified against
the parasail mapping on already-exposed structures and to fail closed otherwise.

`C0_MAPPING_EQUIVALENCE.json` reports `23 / 40` agreement. The premise is
**false**. Inspection of the exposed corpus shows the record sequence is
systematically shorter than the entity canonical sequence — `162` versus `164`
for `182l`, `357` versus `358` for `1a05`, `332` versus `348` for `1a0i` — and
is not even a prefix for `1a0q`.

## 2. Actual P1B provenance

`scripts/structure_sources/biolip.py` parses `BioLiP.txt` and takes
`sequence = columns[20]`, the BioLiP **receptor sequence**.
`scripts/build_holo_complex_index.py` then passes that string to
`_protein_sequence_mapping`, which aligns the resolved structure residues to it
with parasail. So the P1B slot denominator is the BioLiP receptor sequence
length, not the entity length, and the index is an alignment coordinate, not
`label_seq_id - 1`.

This is a factual correction about where a frozen input came from. It is not a
relaxation of any scientific Gate.

## 3. Correction

`M3` and `M4` of the parent are replaced by:

```text
M3'  the protein sequence is the BioLiP2 receptor sequence, column 20 of
     BioLiP.txt, for the (pdb_id, receptor_auth_asym_id) of the BioLiP row
M4'  sequence_index(r) is the parasail alignment coordinate produced by
     nw_trace_striped_16(structure_sequence, sequence, 10, 1, blosum62),
     where structure_sequence is built from the resolved residues in
     ascending label_seq_id using gemmi one-letter codes, exactly as
     scripts/build_holo_complex_index._protein_sequence_mapping does
```

Consequently the **system key** also changes, because P1B enumerates BioLiP
rows rather than scanning asym pairs:

```text
S0'  a candidate system is one BioLiP row, giving
     (pdb_id, receptor_auth_asym_id, ligand_comp_id, ligand_auth_asym_id,
      ligand_auth_seq_id)
S1'  protein atoms are those with auth_asym_id == receptor_auth_asym_id; the
     row is excluded unless they resolve to exactly one label_asym_id
S2'  ligand atoms are those with auth_asym_id == ligand_auth_asym_id, matching
     comp id and auth_seq_id
```

Using column 20 remains `ANNOTATION_ONLY`: it is a sequence string, not a
numeric annotation and not an affinity field. Affinity reads stay `0`.

`M5` is unchanged in form — `slot(r) = min(127, sequence_index(r) * 128 // L)` —
but `L` is now the BioLiP receptor sequence length, matching P1B exactly.

## 4. Re-verification requirement

The parent's fail-closed check is retained and re-pointed: the corrected path
must reproduce `scripts/build_holo_complex_index._protein_sequence_mapping`
**slot for slot** on a sample of already-exposed structures. Anything short of
full agreement is `SLOT_ROUTING_ESTIMAND_INVALID`, and this amendment may not be
amended again to rescue it.

## 5. Unchanged

Everything else in the parent stands: the exposure policy, the untouched-corpus
restriction, the additive blacklist, the `6.0 A` contact contract, the selection
rules `S2`-`S8`, the exclusion and no-replacement policy, the census
definitions, the closure rules `E1`-`E4`, the union closure with its
DataSAIL-style fallback, the C0 Gates `G0a`/`G0b`/`G0c` and their power
calculation, the C1 Gates `C1a`/`C1b`/`C1c`, every threshold and margin, the
seeds, the terminal verdict set and the stopping rules.

No Gate, threshold, margin, seed or verdict is altered by this amendment. The
prior C0 census computed under the rejected mapping is superseded and its
artifacts are rebuilt from scratch.

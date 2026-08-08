# S2 — Frozen 3D Mechanistic Teacher Contract

Registered before any teacher value was used in an evaluation.
No affinity label, no target ID, no ligand ID, no PDB ID, no document ID and no
dataset identity enters the teacher. It is a deterministic function of raw
experimental holo coordinates and standard atom chemistry only.

## 0. Inputs and observability

| quantity | source | observable at deployment? |
|---|---|---|
| protein heavy-atom coordinates + element + residue name + atom name | mmCIF `atom_site` | only with a protein structure |
| ligand heavy-atom coordinates + element | mmCIF `atom_site` (non-polymer) | only with a pose |
| formal aromaticity / ring membership of ligand atoms | RDKit on the CCD-derived SMILES, mapped by element+geometry | yes, from 2D |
| waters, metals, cofactors | `HOH`, metal elements, other non-polymer entities | **explicitly excluded**; recorded as unobserved |
| protonation states | not observable in X-ray at these resolutions | **not modelled**; H atoms ignored |

Waters, metals and cofactors are excluded rather than modelled. This is a stated
limitation: channels 1 and 2 therefore describe **direct** polar contacts only and
cannot represent water-mediated bridges.

## 1. Channel definitions (fixed before evaluation)

For a protein heavy atom `p` and ligand heavy atom `l` at distance `d`, write
`f_cut(d; a, b)` for the smooth cutoff that is `1` below `a`, `0` above `b`, and a
cosine ramp between. All channels use `f_cut` so the teacher is continuous in the
coordinates.

| # | channel | chemistry | distance | orientation | range |
|---|---|---|---|---|---|
| 1 | `hbond_directional` | `p` in {N,O} with a bonded heavy neighbour, `l` in {N,O} | `f_cut(d; 3.2, 4.0)` | `cos^2(theta)` where `theta` is the angle at the donor/acceptor between its bonded-neighbour centroid and the partner; `0` if no neighbour | `[0, inf)` |
| 2 | `electrostatic_signed` | formal charge proxy: protein `Arg NH*/NE`, `Lys NZ`, `His ND1/NE2` = `+1`; `Asp OD*`, `Glu OE*`, C-terminal `OXT` = `-1`; ligand N with no bonded O = `+1` if in a guanidinium/amine-like local pattern, ligand carboxylate/phosphate/sulfonate O = `-1` | `f_cut(d; 4.0, 6.0)` | none | signed real |
| 3 | `hydrophobic_burial` | both atoms carbon or sulphur, neither bonded to N/O | `f_cut(d; 4.0, 5.0)` | none | `[0, inf)` |
| 4 | `aromatic_orientation` | `l` in a ligand aromatic ring, `p` in a Phe/Tyr/Trp/His ring | `f_cut(d_ring; 4.5, 6.5)` between ring centroids | `cos^2(alpha)` between ring normals, reported separately for face-to-face (`alpha < 30 deg`) and edge-to-face | `[0, inf)` |
| 5 | `steric_overlap` | any heavy pair | `max(0, r_p + r_l - 0.4 - d)` with van der Waals radii | none | `[0, inf)` |
| 6 | `pocket_burial` | ligand atom vs all protein heavy atoms | count of protein heavy atoms within `6.0 A`, normalised by ligand heavy-atom count | none | `[0, inf)` |

Van der Waals radii (A): C 1.70, N 1.55, O 1.52, S 1.80, P 1.80, F 1.47,
Cl 1.75, Br 1.85, I 1.98, else 1.70.

## 2. Invariances, asserted and tested

| property | why it holds | test |
|---|---|---|
| rotation invariant | every term depends only on distances and on angles between internal vectors | random `SO(3)` rotation reproduces channels to `< 1e-9` |
| translation invariant | same | random translation reproduces to `< 1e-9` |
| atom-permutation invariant | channels are sums over unordered pairs | random atom reordering reproduces to `< 1e-9` |
| deterministic | no randomness, no learned parameters | repeat run reproduces bitwise |
| label-free | no affinity, no identity input | the function signature accepts coordinates and elements only |
| reflection | **NOT invariant by design for channel 4**; ring-normal `cos^2` is reflection invariant, so in practice all six channels are reflection invariant. Recorded, not claimed as a virtue. |

## 3. Outputs

Two artifacts per complex:

- `pair_local`: for every ligand-atom / protein-residue pair within `8 A`, the six
  channel contributions. This is what a pair-local student must reproduce.
- `aggregate`: the six channel means over the complex, each an invariant scalar.

Missing state (no ligand, no protein atoms in range, unparsable CCD) yields an
explicit `null` and the complex is dropped with a recorded reason, never
imputed as zero.

## 4. What this teacher is not

It is **not** an energy, not a free energy, and not an affinity. It is a set of
six named, reproducible geometric-chemical statistics of an experimental
complex. Any later use in an affinity stage must call the resulting quantities
structural statistics, never binding energetics.

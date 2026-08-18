# Q1 interpretation addendum (per independent review, 2026-08-18)

The Q1_SELECTIVITY.json numbers are unchanged; this addendum fixes their
interpretation. The three passing probes are NOT three independent protein
representations:

- `pair_centered_local_esm` — a genuine protein representation
  (ESM-2 local-window embedding centred at the mutation pair).
- `mutation_position_only` — a mutation-position descriptor, not a protein
  representation.
- `substitution_type_only` — a substitution/edit-type descriptor, not a
  protein representation.

Therefore Q1 is re-labelled into three sub-audits:

| Sub-audit | Content | Status |
|---|---|---|
| Q1-A protein representation capability | pair_centered_local_esm (+0.189 [0.033,0.363]); klifs_pocket (-0.086, not significant); local_onehot_window (+0.027, ns); residue_identity_context (+0.054, ns); global_esm / composition / parent_id / family_id (0.000) | PARTIAL PASS: exactly one genuine protein representation passes |
| Q1-B shortcut / edit-descriptor audit | mutation_position_only (+0.110 [0.021,0.230]); substitution_type_only (+0.209 [0.007,0.420]); edit_descriptor (control) | both descriptors pass; they must be treated as shortcuts, never as biological protein signal |
| Q1-C random / probe capacity audit | random representation 0.000; MLP-8 capacity does not inflate the ESM result (0.919 vs 0.954 linear); random-label curve 0.451/0.446/0.644 | PASS with caveat: substitution_type_only random-label probe reaches 0.644 (small-sample overfit) |

Consequence: the admissible Q1 claim is reduced to: pair-centered local ESM
reads pocket-membership-relevant information under leave-one-parent-out.
`substitution_type_only` positive selectivity is NOT interpreted as protein
signal (it plausibly reflects substitution-type priors or task label
structure), and Q1 is never used as a substitute for Q2.

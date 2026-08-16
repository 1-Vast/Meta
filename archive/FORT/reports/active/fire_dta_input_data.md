# FIRE-DTA Input Data Handoff

## Verdict

`FIRE_DTA_INPUT_DATA_READY`

The frozen label-blind PLINDER registry contains 9,264 training systems and 358 validation systems.
Official link provenance certifies 1,828 experimental-apo-linked training systems; predicted structures
remain a separate state class and were not counted toward the apo gate.

Selective extraction and full parsing produced 12,018 validated caches: 9,622 bound systems, 1,628
experimental apo structures and 768 predicted structures. They contain 41,820,542 receptor heavy atoms,
260,933 ligand heavy atoms, ligand bonds and exact chain/residue/atom identities. The cache inventory is
`10f88a116b38d545aeb32a381e5efa5b3bd3253cd96094248adeea1a3126b11f`.

The data boundary is intentionally model-independent. Pocket crops, graph edges, alignments and learned
node features have not been created, so BFEO/BERP can be adjusted without changing the frozen split or
redownloading structures. No model forward pass or training was run. `sealed_test_consumed=false`.

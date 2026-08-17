# Stage W0-P preregistration — local point-mutant positive-control panel

Frozen **2026-08-17, before any W0-P label statistic was read.** Stage W0b
found no named mutation/ortholog panel locally and issued NO-GO for W1
biological interpretation. This stage constructs the closest admissible local
positive control: near-identical BindingDB-Ki target-sequence pairs that
differ by one or a few residues and share measured ligands, then tests whether
a low-capacity gradient-trained model can recover the signed affinity change
across the protein pair.

## 0. Candidate definition (label-blind)

From the governed BindingDB-Ki `meta_train` + development-validation protein
sequences (`proteins.jsonl`, 499 sequences):
* CD-HIT at 98% identity (`-c 0.98 -n 5`) finds near-identical sequence
  clusters;
* a candidate pair has equal length, 1–5 residue mismatches, and >= 3 shared
  ligand identities with governed Ki cells in the development splits;
* the pair is labelled `A` (lexicographically smaller target id) and `B`;
* mutation positions are the zero-based mismatch positions.

This is a **point-mutant / isoform positive-control candidate**, not a
labelled gatekeeper panel. It is disclosed as such; it does not license any
claim about the functional role of the mutated positions.

## 1. Frozen controls and estimand

For each shared ligand in a candidate pair, the label is the signed change
`delta_y = y_B(ligand) - y_A(ligand)`.

Controls for every pair:
* correct mutation positions;
* random positions of equal count;
* BLOSUM-distance-matched but spatially unrelated positions;
* global pooled protein embedding;
* capacity-matched random protein representation.

A positive-control test passes only if a gradient-trained, deliberately
low-capacity bilinear model (`phi(ligand) * psi(mutation-position features)`)
predicts `delta_y` on leave-one-ligand-out within the panel with positive
Pearson/Spearman/sign accuracy, and the correct-position arm beats all
corrupted-position arms with a cluster-bootstrap lower bound above zero.

## 2. Execution order

1. Build `W0P_PANEL.json`: pairs, mutation positions, shared-ligand counts,
   censoring status of the underlying BindingDB cells, and label-sign
   consistency (descriptive).
2. If the panel has fewer than 3 pairs or fewer than 20 ligand rows in total:
   W0-P is recorded `INSUFFICIENT` and W1 remains NO-GO.
3. Otherwise run the gradient-trained positive-control model (ordinary AdamW;
   no ridge, no closed form) with the controls above. Single seed can reject;
   three seeds required for a pass.

## 3. Claim boundary

A pass shows the pipeline can recover protein-conditioned signed changes on
near-identical local sequences with shared ligands. It does **not** prove the
mutated positions are binding-site determinants, and it does not by itself
authorize W1; W1 also requires the W0b censoring and support re-census gates.

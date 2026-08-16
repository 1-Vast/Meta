# RDIB / PD-MVR multi-agent deep review

Date: 2026-07-29  
Status: empirical gates executed; both proposed primary routes stopped

## Combined verdict

| Route | Empirical result | Mathematical result | Current status |
| --- | --- | --- | --- |
| RDIB exact replicated differences | 145 blocks, 59 targets, optimistic packing ceiling 56, only 5 same-ChEMBL-target dual-ligand blocks | The proposed 448-dimensional pair-specific covariance has rank at most one with two lineages; exact construct and source independence are not established | Stop exact-pair route; audit recurring directed edits only under a new frozen contract |
| PD-MVR structural rectangles | The same BioLiP2 result misses 200 blocks/80 targets and has 17.24% maximum PubMed share | Packing cannot create independent information absent from the candidates | Stop structural-block route |
| PD-MVR exact bridges | 23 exact PDB-sequence-accession-ligand candidates, 5 targets, optimistic ceiling 5, 2 known homology components | The proposed bridge loss contains no affinity `Y`; permuting bridge `Y` leaves the objective and gradients unchanged | Stop bridge and missing-view model |

## Identifiability corrections

### RDIB

1. `8 x 8 x 7 = 448` output dimensions cannot receive a meaningful
   pair-specific covariance estimate from two lineages. Any later version
   would need training-only pooled, block-diagonal, or shrinkage covariance.
2. If the role matrices are signed SVD coordinates, `A^T C B` is not a
   probability and BCE is invalid. The roles must be non-negative normalized
   memberships or the loss must be continuous.
3. BioLiP sequence equality is weaker than an exact construct. A valid
   construct contract needs mmCIF entity, mutations, missing residues,
   assembly, cofactors, water, altloc, protonation, and deposition-series
   rules.
4. Sequence plus SMILES maps to a distribution of possible conformations and
   contact states. Contact-difference supervision is identifiable only if
   cross-lineage signal exceeds PLIP version, hydrogenation, water, altloc,
   and conformational noise.
5. Stage 2 would need a frozen target-only plus ligand-only additive null,
   exact no-target-bypass interventions, and train-reference centering for
   unseen targets.

### PD-MVR

The proposed anchor loss is

`||W_C q_C(C_tl) - z_tl||^2 + retrieval`,

which contains no bridge affinity `Y_tl`. Consequently a bridge-affinity
permutation changes neither the objective nor its gradient. The loss can
align contact with a deployable sequence-ligand representation but cannot
identify a contact-affinity common space.

The proposed common-rank check likewise tests `C <-> X`, not `C <-> Y`.
A shared latent can split into orthogonal contact-private and
affinity-private subspaces while fitting every observed loss. Absolute
contacts, same-target contact differences, and affinity `2 x 2` double
differences also occupy different nuisance quotient spaces. A repaired
claim would require contrast-complete dual-view cycles or a full-rank
connected anchor design, followed by an affinity readout restricted to the
cross-fitted common subspace.

## Innovation boundary

The broad architectural claims are already covered:

- [LINKER](https://arxiv.org/abs/2509.03425) predicts a
  residue-by-functional-group-by-seven-interaction tensor from protein
  sequence and ligand SMILES using structure-derived supervision.
- [PBCNet2.0](https://www.nature.com/articles/s41589-026-02241-x) is a
  pairwise relative-affinity model trained on millions of complex pairs.
- [MVAE](https://papers.nips.cc/paper_files/paper/2018/hash/1102a326d5f7c9e04fc3c89d0ede88c9-Abstract.html)
  already trains across subsets of observed modalities.
- [Collective matrix completion](https://proceedings.mlr.press/v38/gunasekar15.html)
  already provides shared factors across relation matrices under explicit
  rank, incoherence, connectivity, and sampling assumptions.

RDIB therefore cannot claim the first interaction map, pairwise delta,
contact-to-affinity transfer, or interaction bottleneck. Its only plausible
novelty is the narrow, falsifiable combination of exact-construct directed
interaction changes, independently replicated provenance, strict
no-target-bypass deployment, and dual-cold evaluation.

PD-MVR's plausible contribution is a provenance-disjoint data-topology and
fail-stop protocol, not a new missing-modality architecture. That protocol
would still need real contrast-complete bridges; the current local sources
do not provide them.

## Authorized continuation

No formal model is authorized. The recurring-edit investigation subsequently
returned only 29 same-target cross-PubMed units with an optimistic ceiling of
20 and is now also stopped. The remaining public-data investigations are:

1. an official-schema and evidence-lineage audit of the BioLiP 2026
   LLM-rule affinity layer;
2. discovery of a materially new source with exact row mapping, endpoint,
   construct, provenance, and target-domain coverage.

If these do not introduce real observed information, the executable
successor remains the already frozen prospective cycle-closing acquisition
design. More adapters, completion losses, shared latents, or representation
changes cannot repair the identified observation deficit.

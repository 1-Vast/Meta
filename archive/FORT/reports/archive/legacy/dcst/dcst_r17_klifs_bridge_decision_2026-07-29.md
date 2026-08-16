# DCST-R17 KLIFS exact-overlap bridge decision

Date: 2026-07-29  
Decision: `STOP_KLIFS_BRIDGE_INADEQUATE`

## Result

The label-blind firewall started from 11,250 valid human KLIFS mechanisms and
retained 5,161 complexes, 181 accessions, 1,681 ligands, and 1,340 scaffolds.
No ChEMBL affinity column was loaded.

Two quality conditions passed:

- 42 development targets had an aligned KLIFS pocket and all had distant
  aligned-pocket identity support at least `0.25`, with zero `0.40` 4-mer
  firewall violations;
- repeated exact target-ligand structures had median IFP Jaccard `0.972953`.

The interaction and bridge conditions failed:

- iterative target-degree/ligand-degree 2-core: 393 edges, 92 accessions,
  120 recurring ligands, but only one connected component because promiscuous
  ligands joined the graph;
- exact allowed KLIFS-to-ChEMBL-train bridge: 300 pairs, 50 targets, and 45
  homology components, below the frozen 300/60/50 joint requirement;
- only `7.41%` of 20,000 ChEMBL-train ligands had maximum allowed-source
  Tanimoto at least `0.40`.

Runtime was `39.984 s` on CPU-only RDKit/graph preprocessing. Confirmation
metadata defined the firewall; confirmation features and labels were not
loaded.

## Interpretation

KLIFS solves the target-coordinate and repeated-structure quality problems,
but it does not solve broad chemical support or independent rectangular
topology. The exact bridge is concentrated in 50 kinase targets, and its
recurring-ligand graph is one provenance-sensitive component. Training the
proposed RIFB model would therefore repeat the prior kinase-only/MNI-0
limitations with a smaller affinity bridge.

The route stops before dataset materialization or affinity access. The next
source audit must add broad, independently curated protein and chemical
coverage with explicit complex identity, rather than increase model capacity.


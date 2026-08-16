# Post-UBSE-G1 multi-agent failure and route review

Date: 2026-07-29  
Status: completed independent implementation, failure-forensics,
identifiability, prior-art, and route-selection review

## Executive decision

UBSE-G1 is a valid empirical stop, not an optimization failure. The model
memorized the fit contact panels, did not transport to closed validation or
audit domains, and did not establish residue-position use. The next paper
route must add a new ligand-conditioned information source.

The retained program is:

1. learn and protect a target-marginal pocket anchor only as a proposal
   component;
2. obtain real residue-by-functional-group-by-interaction-type events from
   experimental complexes;
3. distill those events into a deployment-matched monomer-structure plus
   ligand optimal-transport student; and
4. use a fixed-measure ANOVA purification only as a nuisance-removal and
   exact-null protocol.

Neither pocket anchoring nor orthogonal decomposition alone supplies missing
pair information. ChEMBL affinity and Stage-2 remain locked.

## Independent G1 implementation audit

The accepted G1 ledger has 81 contrast audit panels times three seeds times
four controls, or 972 complete rows. An independent agent verified:

- no fit/validation/audit overlap in homology component, scaffold, or PubMed;
- no missing or duplicate `(seed, panel, control)` cell;
- identical null/cross trainable parameters, initialization, panel order,
  and update count;
- parameter-free cross multiplication;
- matched within-panel ligand cycling and mask-aware position destruction;
  and
- preregistered tie handling and panel metrics.

The original point check for S3-S5 subtracted cross-seed medians rather than
using the preregistered panel-over-seed delta. The bootstrap already used the
correct estimand. Corrected estimates were 0.0219, 0.0263, and -0.0209; every
gate remained false. The result therefore did not require retraining.

## Failure forensics

Re-evaluating the saved weights showed:

| Arm | Fit directional / cosine | Validation | Audit |
| --- | ---: | ---: | ---: |
| cross | 0.960-0.962 / 0.921 | 0.448-0.544 / -0.117-0.067 | 0.491-0.537 / 0.006-0.096 |
| additive null | 0.725-0.733 / 0.550-0.557 | 0.513-0.575 / -0.019-0.167 | 0.458-0.534 / -0.078-0.050 |

On the fit domain, ligand derangement reduced cross directional accuracy to
approximately 0.12 and the position-free control to 0.61-0.62. The model had
ample capacity to memorize panel-position interactions. Increasing width,
epochs, or backbone capacity is therefore not a justified rescue.

### Sparse residual and contact-count shortcut

Among the 81 contrast audit panels:

- 69 had two ligands;
- the median number of contrast residues was three;
- 53 panels had at most three contrast residues;
- 71 had at most five;
- 58/127 ligand pairs had nested contact sets; and
- 38/81 panels' differences were completely ordered by total contact count.

A label oracle that assigns each ligand the same score at every residue,
using only the ligand's contact-set size, obtained directional accuracy
0.814. Directional accuracy alone therefore does not certify residue
placement.

The random tie-aware assignment baseline was 0.4679. G1 cross obtained
0.5201, while the position-free control reached 0.5355.

Morgan similarity versus contact Jaccard had Spearman correlation 0.084
(`p=0.346`), and descriptor distance versus contact Jaccard had correlation
-0.058 (`p=0.514`). The selected two-dimensional covariates were weakly
aligned to the small positional residual.

### Double-centered exploratory diagnosis

For matrix \(Y\), define:

\[
D(Y)=Y-\overline Y_{l\cdot}-\overline Y_{\cdot i}+\overline Y.
\]

Only 45/81 audit, 36/57 validation, and 771/1,138 fit contrast panels retained
nonzero position-reallocation labels after this operation. Cross validation
double-centered cosine was -0.041, -0.056, and -0.137 across the three seeds.
The audit seed variation was too large to reinterpret as transport.

G1's exact null should therefore have included a position-free
target-by-ligand contact-burden scalar. The destruction control caught this
shortcut, but a future exact null must model it explicitly.

## Pocket-anchor diagnosis

For the target-marginal union pocket, the G1 additive arm retained:

- validation AP 0.195-0.213 and AUC 0.753-0.766;
- audit AP 0.208-0.231 and AUC 0.765-0.773; and
- audit oracle-size top-k recall 0.198-0.228.

Random prevalence was approximately 0.04. The cross arm reached fit AP about
0.59 and AUC about 0.96, but its validation and audit AP fell to roughly
0.13-0.16 and 0.12-0.14. Joint interaction training corrupted a modest
transferable anchor.

This supports large fit-side target-marginal pretraining with the anchor
frozen or otherwise protected. It does not show ligand specificity. P0A is
therefore an engineering proposal gate, not the new information claim.

## Orthogonal binding-state decomposition

For a complete target-specific ligand-by-position rectangle:

\[
g_t=\overline Y_{t\cdot\cdot},\quad
\mu_{ti}=\overline Y_{t\cdot i}-g_t,\quad
c_{tl}=\overline Y_{tl\cdot}-g_t,
\]

\[
r_{tli}=Y_{tli}-g_t-\mu_{ti}-c_{tl}.
\]

The residual satisfies zero row and column sums. This is classical two-way
ANOVA/double centering, not a new representation theorem. Relevant prior art
includes [Mandel's two-way table analysis](https://nvlpubs.nist.gov/nistpubs/jres/73b/jresv73bn4p309_a1b.pdf),
[Gabriel's biplot analysis](https://academic.oup.com/jrsssb/article/40/2/186/7027470),
and [functional-ANOVA purification](https://proceedings.mlr.press/v108/lengerich20a.html).

For unbalanced/masked event tensors, the measure must be frozen on fit:

\[
\langle X,Z\rangle_W=\sum_{li}w_{tli}X_{tli}Z_{tli}.
\]

Weighted row/column constraints must be solved as a projection. Test-batch
centering is forbidden because it would make a single prediction depend on
other test ligands. An audit target's marginal pocket and pair burden must
also be predicted from deployable inputs, never computed from its observed
holo contacts.

If \((g,\mu,c,r)\) are a deterministic invertible transformation of the same
coarse \(Y\), no information has been added. The decomposition is retained
only to:

- protect the pocket anchor;
- expose pair contact burden in the exact null;
- prevent contact-count evidence from being called spatial interaction; and
- define typed-event destruction controls.

Relative-contact learning itself is not novel; nearby examples include
[DeltaDelta](https://pmc.ncbi.nlm.nih.gov/articles/PMC7066671/) and
[PBCNet](https://www.nature.com/articles/s43588-023-00529-9).

## Route ranking

### A: typed three-dimensional event teacher plus deployable OT student

This is the primary route. The teacher target is:

\[
E_{rfgk}=
\text{residue }r
\times\text{ ligand functional group }fg
\times\text{ event type }k.
\]

The student may receive only monomer/predicted target structure and ligand
inputs available at deployment. Partial or unbalanced optimal transport with
a dustbin maps target surface patches to ligand functional groups. Holo
coordinates are teacher-only.

The innovation cannot be stated as generic structure pretraining or
interaction-map prediction. Those ideas are already covered by
[ATOMICA](https://pmc.ncbi.nlm.nih.gov/articles/PMC12026499/),
[LINKER](https://pubs.acs.org/doi/10.1021/acs.jcim.6c00527),
[NeuralPLexer](https://www.nature.com/articles/s42256-024-00792-z),
[DynamicBind](https://www.nature.com/articles/s41467-024-45461-2),
and large synthetic-structure affinity pretraining such as
[GatorAffinity](https://pmc.ncbi.nlm.nih.gov/articles/PMC12621780/).

The defensible combined contribution is:

> source-closed typed-event distillation, deployment-matched
> monomer-to-functional-group optimal transport, and a fixed-measure
> purification/exact-null protocol under strict dual-cold evaluation.

UBSE-A0 has prepared 3,467 complete coordinate locators over 2,833 PDB
entries, with zero PDB overlap among roles. Remote coverage is waiting on
network-enabled coordinate fetch.

### B: target-marginal pocket anchor

P0A uses 62,849 fit-side closed rows and 38,781 targets after removing held
homology, scaffolds, and PubMed sources. A pass can only freeze a target-only
proposal model for A. It cannot be cited as pair information.

### C: frozen pair-state or cofold probe

A frozen no-affinity pair-state model could cheaply test whether a pretrained
cofolder already contains a deployable interaction state. It must first bind
the checkpoint hash, training cutoff, and membership audit, and must beat
additive and partner-destruction controls. No admissible checkpoint is
currently local; shell network access is denied. This remains an upper-bound
probe, not an executed result.

## Required future controls

Before typed-event residual training:

1. exact null includes target pocket plus typed pair-burden \(C_{tlk}\);
2. wrong ligand, wrong protein, event-type shuffle, and structure-position
   destruction are matched;
3. event reliability is measured across PubMed/PDB repeats before model fit;
4. no test-batch centering or target support label is used;
5. audit event labels are extracted only after split and extractor freeze;
6. all Stage-1 checkpoints inherit source membership;
7. Ki and Kd remain separate in Stage-2; and
8. affinity, confirmation, and sealed outcomes remain inaccessible until the
   event/student semantic gates pass.

## Current authorization

- G1: stopped.
- OBD/ANOVA: retained only as a purification and exact-null protocol.
- P0A: authorized as a target-only CUDA engineering gate.
- A0 coordinates: locally addressable; remote fetch waiting.
- A1 typed-event extraction: not authorized until A0 remote coverage passes
  and an extractor contract is frozen.
- Stage-2 affinity: locked.


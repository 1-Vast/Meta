# UBSE, VBS, and AdaMBind multi-agent deep review

Date: 2026-07-29  
Review roles: mathematical identifiability, data topology, meta-learning task
definition, falsification controls, recent prior art, and executable gating

**Supersession note:** this review admitted UBSE-G1 prospectively. G1 later
returned `STOP_UBSE_G1_NO_DEPLOYABLE_INTERACTION_RESIDUAL`; the current
failure analysis and successor boundary are recorded in
`post_ubse_g1_multiagent_failure_and_route_review_2026-07-29.md`.

## Unified verdict

The three supplied suggestions must be separated:

1. **UBSE / conditional bound-state prediction:** survives, but only in the
   narrow form now admitted to G1. G0R and G0PB provide the first current
   public-data evidence for a reproducible, independently packable,
   ligand-conditioned residue signal.
2. **Virtual Binding System as a causal ligand intervention model:** stops.
   The current source observes bound complexes but no construct-matched
   unbound state in the same response space.
3. **AdaMBind as the primary strict dual-cold solution:** stops. It is a
   labeled few-shot adaptation method. It remains a conditional secondary
   `k=5` task only after an interaction representation passes.

No conclusion authorizes affinity decoding, coordinates, Stage-2 fitting, or
confirmation/sealed access.

## Why the UBSE question is different from the stopped RDIB routes

RDIB required replicated *differences* between two ligands on the same exact
target and source-independent bridge structure. That topology collapsed to
145 exact pair blocks, 29 recurring target-edit units, and ceilings below the
program floors.

UBSE-G0R instead tested an absolute ligand-conditioned observation with a
direct same-target wrong-ligand control. It found:

- 1,028 cross-PubMed/cross-PDB sequence/connectivity repeat units;
- 800 exact-sequence targets;
- median correct contact Jaccard `0.75`;
- hard wrong-ligand Jaccard `0.50`;
- target-bootstrap ligand-specific margin interval `[0.125, 0.200]`;
- optimistic observed-contact Recall@1 `0.5604` versus random `0.0846`.

This is enough to reject “the label is only a fixed target pocket” as a full
explanation. It is not enough to reject scaffold, ligand-identity, stereo,
construct, or additive two-tower shortcuts.

G0P therefore moved to same-target, same-PubMed, same-scaffold panels. Its
initial frozen run stopped on a `28.0397%` largest component and `6.7618%`
homology share. The single preregistered removal-only correction removed two
overrepresented homology blocks and passed all unchanged gates:

- 1,412 panels / 3,383 within-panel ligand contrasts;
- 492 homology components;
- largest conflict component `17.6346%`;
- conflict-free packing 450;
- frozen independent audit 88;
- closed residual training substrate 1,324 panels.

That is sufficient to run the affinity-blind G1 student, not to skip it.

## Two mathematical counterexamples

### Pair retrieval is not an interaction certificate

Let:

\[
z(t,l)=[e_t,e_l].
\]

This representation has no interaction term. It can nevertheless give
perfect pair retrieval, and both wrong-target and wrong-ligand replacements
can destroy its match. Protein-only, ligand-only, and centroid baselines do
not rule out this concatenated identity shortcut.

G1 therefore uses the exact additive logit null:

\[
\eta_A(t,l,i)=a(u_{ti})+b(v_l),
\]

and an identical-parameter interaction arm:

\[
\eta_X(t,l,i)
=\eta_A(t,l,i)+\frac{u_{ti}^{\mathsf T}v_l}{\sqrt{64}}.
\]

Only held-panel centered improvement over \(\eta_A\), plus ligand and
protein-position destruction, can admit the representation.

### VBS transition gauge

For:

\[
z^B_{t,l}=z_t^0+\Delta z_{t,l},
\]

any target-specific \(h(t)\) gives:

\[
z_t^{0\prime}=z_t^0+h(t),\qquad
\Delta z'_{t,l}=\Delta z_{t,l}-h(t),
\]

with the same bound output. Freezing a sequence encoder merely chooses a
convention; it does not identify an observed ligand-induced change.

The virtual-cell analogy is incomplete because methods such as
[CellOT](https://www.nature.com/articles/s41592-023-01969-x) observe control
and perturbed distributions in the same expression space. BioLiP does not
observe a ligand-free residue-by-functional-group contact tensor
\(C(t,0)\). The defensible name is *conditional bound-state predictor*, not
causal perturbation world model.

## AdaMBind boundary

[AdaMBind](https://www.nature.com/articles/s41467-026-70554-5) uses a
40%-sequence-identity novel-target split, so it should not be dismissed as a
mere random target split. But every test target receives 5 or 40 labeled
support pairs before query prediction. Its within-target support/query
construction does not jointly close ligand scaffold/chemical neighbourhood,
document/provenance, assay, and endpoint; no-support true zero-shot is left as
future work.

For the current project:

- primary `k=0` strict dual-cold: no AdaMBind adaptation;
- secondary `k=5`: same exact target and endpoint, but support/query scaffold,
  chemical-neighbour, document/assay, and provenance disjoint;
- mandatory controls: zero support, correct support, exact wrong-target
  support, shuffled labels, ligand derangement, protein-free adaptation, and
  mean-only calibration;
- report `k=0` and `k=5` on the identical query set.

The existing balanced rectangle audit already indicates that this legal
episode topology is scarce. G1 and a later source/topology gate must pass
before any MAML implementation.

## Prior-art boundary

The general interface-foundation-model and sequence-to-interaction-map ideas
are not new:

- [ATOMICA](https://pmc.ncbi.nlm.nih.gov/articles/PMC12026499/) trains
  atom/block/interface representations with geometric denoising and masking
  over 2,037,972 interaction complexes and constructs interface networks.
- [LINKER](https://pubs.acs.org/doi/10.1021/acs.jcim.6c00527) predicts a
  residue-by-functional-group-by-seven-interaction-type tensor from protein
  sequence and ligand SMILES, then freezes the interaction module for
  affinity prediction.
- [NeuralPLexer](https://www.nature.com/articles/s42256-024-00792-z) already
  provides ligand-conditioned complex/bound-state modeling.

Accordingly, “universal interface atlas,” masked interaction tokens,
sequence/SMILES interaction prediction, cross-attention, or a VBS name cannot
carry the paper's novelty alone.

The remaining defensible contribution is the combined identification
protocol:

> cross-provenance repeatability weighting + same-scaffold,
> same-publication centered contact residuals + explicit additive exact null
> + matched ligand/protein destruction + strict
> homology/scaffold/provenance dual-cold evaluation + a later Ki/Kd
> exact-null affinity increment.

This is a data-identification and falsification contribution, with a compact
multiplicative interaction module. It should not be marketed as inventing
interface foundation models or causal binding dynamics.

## Active execution boundary

UBSE-G1 is frozen to require:

- cross directional accuracy at least `0.60`;
- centered cosine at least `0.10`;
- cross-minus-additive directional increment at least `0.05` with bootstrap
  lower 95% bound above zero;
- at least `0.05` destruction by same-panel ligand derangement and by removing
  residue-specific protein information, each with positive lower bound;
- assignment Recall@1 at least `0.60` and at least `0.05` above every control;
- at least two of three seeds with positive exact-null and both destruction
  increments.

If G1 passes, it requests only a ChEMBL-TRAIN semantic-support gate. That gate
must still establish fold-specific provenance closure, approximately 423
independent predictive-scale support units, and noncollapsed cross residuals.
Only after that may a frozen shallow Stage-2 null test be registered:

\[
\widehat y=B0+\theta^\mathsf{T}r_{t,l},\qquad H_0:\theta=0.
\]

A deep affinity head, coordinate teacher, causal transition claim, AdaMBind
adaptation, and outcome access remain unauthorized.

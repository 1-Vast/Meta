# CIPF + TERM failure audit

## Decision

The corrected three-seed nested-k development run does not admit CIPF+TERM as
an effective cold-target few-shot mechanism. Most few-shot MSE improvement is
target-level calibration. TERM adds only 0.0008 MSE reduction at k=1 and
0.0112 at k=5; every component-bootstrap interval crosses zero. Cyclic support
label permutation is indistinguishable from correct binding. The current
candidate therefore cannot support SOTA, excellent-performance, or
important-mechanism-source claims.

## Causal diagnosis

The failure is a chain rather than a single numerical bug:

1. CIPF constructs all ligand-atom x full-protein-residue pairs without a
   pocket, contact, distance, or locality mask.
2. Learned primitive indices have no permutation/sign/scale identifiability or
   governed chemical anchors.
3. TERM receives `g=-r*phi`, but the router can bypass it through ligand-change,
   primitive-identity, and absolute primitive inputs.
4. The confidence entropy is detached before reliability and has no auxiliary
   loss. Consequently `protein_prior`, `confidence`, and `log_temperature`
   receive no gradient; an isolated backward check returns `grad=None` for the
   entire branch.
5. Regression, ranking, and centered-shape losses can all improve the shared
   zero-shot trunk. Only k>=2 receives explicit dependency loss, for about 144
   training episodes per seed.
6. Checkpoint selection minimizes full validation MSE, not TERM-vs-cut or
   real-vs-permuted evidence. It can select a strong level model with a dead
   mechanism path.

## Module findings

### Data and protocol

- CD-HIT40 component splitting, support/query cell isolation, unique ligand
  identities, and nested common queries are correctly enforced.
- The run uses only 120 steps x 2 episodes. Approximately 48 episodes are seen
  for each k in `{0,1,2,3,5}` while all encoders, CIPF, zero-shot, level, and
  TERM parameters are trained from scratch.
- The final population contains 42 targets but only six top-level components,
  so confidence intervals are wide. This does not explain the near-zero label
  permutation gap.

### Protein and ligand encoders

- Protein input is pooled/residue PLM state with no pocket localization or
  binding-site supervision.
- The governed 32-dimensional atom vector lacks explicit H-bond donor/acceptor,
  partial charge, ionization, and complementary residue-role contracts.
- Wrong-protein effects are small and uncertain, consistent with weak
  target-specific interaction use.

### BPSF / CIPF

- Pair construction and softmax pooling operate over every atom-residue pair.
  True pocket signal can be diluted by non-contact residues, and global
  softmax discards contact count/coverage.
- `tanh(atom_projection + residue_projection)` is an additive bias, not an
  explicit donor-acceptor, charge-complementarity, or hydrophobic cross term.
- Primitive response is an unbounded linear projection. Coefficient/primitive
  scale is not identifiable, while exact-evidence magnitude changes with that
  arbitrary scale.
- Query-only, episode-local mean/correlation regularization does not align
  primitive meaning across targets, supports, or random seeds.
- Endpoint and primitive heads share the pair trunk. The much stronger
  zero-shot/full objectives can train the trunk without producing reusable
  primitive functions.

### TERM evidence and router

- The exact-gradient identity is locally correct for a virtual squared loss
  around zero-shot `f0`; it does not establish that `phi` is a mechanism.
- Evidence uses `y-f0` but the correction is applied after scalar level
  adjustment, so target-level residual can be explained by both paths.
- `g` is only one router input. Once its global norm opens reliability, the
  router can emit coefficients from label-independent ligand/identity inputs.
- `sum/sqrt(k)` aggregation and a positive `log(1+k)` reliability term increase
  correction scale with k even when label binding is inconsistent.
- Confidence and signed-coefficient heads are nominally separate, but the
  confidence branch is untrainable because entropy is detached and no other
  confidence loss exists.

### Training losses and controls

- Full regression/ranking/shape losses do not isolate TERM; they also optimize
  the zero-shot query path.
- At a label-insensitive solution, real and permuted predictions and their
  shared-path errors coincide, making the margin contrast weak or stationary.
- Cross-component foreign support is an easy distribution intervention and can
  be detected without learning within-support ligand-label assignment.
- k=1 has no direct permutation or magnitude-matched wrong-label dependency
  supervision.

### Cartesian branch

- Production episode materialization supplies no coordinates or geometry
  edges, and the run has no geometry checkpoint. Cartesian code is not a
  performance source in this result.
- Future Cartesian slots also require explicit alignment with CIPF primitive
  indices before they can be added safely.

## Observability gap

The current result stores endpoint metrics but not primitive/TERM internals.
It cannot distinguish primitive collapse, coefficient collapse, reliability
suppression, or cancellation. A subsequent diagnostic run must record, by
seed/k/target:

- primitive mean, RMS, effective rank, and cross-seed slot alignment;
- atom-residue attention entropy/effective support and protein-length relation;
- exact-gradient norm, signed coefficient magnitude, reliability, and SAR
  magnitude;
- per-loss gradient norm/cosine for the pair trunk, primitive response, level,
  triad, and confidence modules;
- real/permuted/matched-wrong differences before endpoint aggregation.

## Required gates before another architecture expansion

1. Fix or separately supervise the dead confidence branch.
2. Run held-out-task synthetic tests comparing oracle-primitive TERM with
   learned-CIPF TERM; include level-only and private-mechanism abstention.
3. Add a label-binding objective that has useful gradient at the
   label-insensitive solution, including a k=1 matched-wrong control.
4. Run a matched A/B/C/D isolation: level, prior meta, CIPF+prior meta, and
   CIPF+TERM. Do not attribute a joint failure to either component without it.
5. Do not add Cartesian complexity until the sequence/2D primitive and router
   pass real-vs-permuted and TERM-vs-cut gates.

## Evidence files

- `report/meta_fewshot/cipf_term_corrected_nested_k01235_3seed_20260814/RESULT.json`
- `report/meta_fewshot/cipf_term_corrected_nested_k01235_3seed_20260814/EPISODES.json`
- `report/meta_fewshot/CIPF_TERM_DECISION_20260813.md`


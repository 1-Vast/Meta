# CIPF + TERM architecture decision

## Scope

The main cold-target method no longer treats pair-specific 3D-complex coverage
as a prerequisite.  It learns how complex models organize information using
sequence/residue embeddings and a 2D ligand graph.  Only two components are
candidate paper innovations:

1. **Complex-Inspired Interaction Primitive Field (CIPF)**;
2. **Triadic Evidence-Routed Meta-Learner (TERM)**.

The Cartesian encoder is an optional bias on the same primitive responses when
a legal common frame exists. It is not a separate required model and no search
for a pose sidecar for every BindingDB cell is part of the active plan.

## CIPF v0

Protein and ligand intra encoders remain shared.  Their rectangular
atom--residue field is read by globally indexed primitive queries. Raw governed
RDKit atom features and learned residue-region features add a weak semantic
compatibility bias. Each primitive emits a learned scalar response
`phi(P,L)[m]` and an embedding. The scalar response excludes the
dictionary-only identity baseline so it cannot become a pair-independent level
expert.

The current X32 atom contract provides element, degree, formal charge,
hybridization, aromatic/ring and chirality features. It does not provide a
complete H-bond donor/acceptor contract. Therefore v0 claims weakly
chemistry-anchored shared interaction functions, not physical energy terms or
explicit hydrogen-bond primitives. A future governed sequence/SMILES
physicochemical sidecar can strengthen these anchors without any 3D input.

## TERM

TERM defines virtual coefficients `a_m` only for the meta correction:

```
f(P,L;a) = f_level(P,L) + sum_m a_m phi_m(P,L).
```

At `a=0`, for an explicitly declared virtual squared loss around the zero-shot
endpoint `f0`, the exact coefficient gradient is `g_im = -r_i phi_im`, where
`r_i=y_i-f0(P,L_i)`. This is not the gradient of either the scalar-level
baseline or the trainer's complete SmoothL1/ranking objective. There is no
inner-loop coefficient update: a learned support-order-invariant router maps
the evidence to signed query-specific coefficients.

The router jointly consumes all three edges:

- support protein--ligand response `phi(P,L_i)`;
- a directional ligand-pair descriptor built from `(h_i,h_q)`;
- query protein--ligand response `phi(P,L_q)`.

Signed coefficients and routing confidence are separate heads. Entropy from
the confidence distribution is detached before the reliability calculation,
so the signed effect of a primitive is not confused with ambiguity.

Support labels never appear directly in query predictions. The reference
anchor, difference score/delta/blend, LOO reconstruction, sensitivity head and
internal BPSF adapters are absent from the active path.

## Zero/one-shot and level semantics

- k=0 hard returns the shared zero-shot endpoint; protein prior cannot create a
  correction without support.
- TERM evidence uses raw zero-shot residuals, retaining non-zero k=1 evidence.
- level calibration is a target-wise scalar shift and cannot change query
  ordering.
- primitive mean/decorrelation regularization discourages TERM from cheaply
  relearning target level; it does not guarantee identifiability.

## Training contract

Training remains one normal episodic end-to-end run. AdamW uses two parameter
groups in the same backward pass: a slower backbone/CIPF group and a faster
TERM group. Pairwise ranking has material weight because prior D-MEMT evidence
showed a larger ordering signal than regression signal. This is not staged
pretraining or target-specific fine-tuning.

## Construction checks and unexecuted falsification gates

The intended full synthetic falsification suite is:

- level-only tasks: mechanism correction must approach zero;
- source-shared primitive tasks: k=1/2 must recover a reusable primitive;
- private random mechanisms: transfer is impossible and the model must not
  claim support predictivity.

The repository currently passes only weaker arithmetic/construction smoke
checks. The shared-primitive check optimizes and scores the same synthetic
tasks and therefore does not demonstrate held-out-task one-shot transfer. The
private and level-only checks do not execute TERM end to end. The
pre-registered multi-seed sealed thresholds, held-out tasks, wrong-label arms,
bootstrap intervals and abstention checks were not executed. No synthetic
admission gate is claimed to have passed.

## Evidence status

A two-step GPU end-to-end smoke completed with finite losses and about 58 MiB
peak allocated memory at the tiny configuration. The next evidence gate is a
fixed-budget A/B/C/D development comparison:

- A: BPSF + level;
- B: prior D-MEMT;
- C: CIPF + D-MEMT;
- D: CIPF + TERM.

Only if D improves both representation and meta-learning contrasts should a
new sealed confirmation population be opened. The full matched A/B/C/D matrix
was not completed: D itself failed its TERM-cut admission gate, so the
pre-registered hard-stop prevented spending more compute on C or opening a
confirmation population.

## Completed BindingDB development result

The three-seed, 120-step, nested common-query k={0,1,2,3,5} run completed on
42 targets from six CD-HIT40 components. Peak allocated GPU memory was 3.57
GiB. The historical meta-test remains a consumed development population. This
is a pre-fix diagnostic: review subsequently found that TERM used
`mean/sqrt(k)` instead of `sum/sqrt(k)` and lacked a zero-evidence abstention
gate. The table is retained as failure evidence for that implementation and
must not be used to judge the corrected candidate.

| k | full MSE | scalar-level / TERM-cut MSE | full CI | full Spearman | cut - full (95% component CI) |
|---:|---:|---:|---:|---:|---:|
| 0 | 3.2834 | 3.2834 | 0.5315 | 0.0811 | 0 |
| 1 | 2.0555 | 2.0393 | 0.5241 | 0.0740 | -0.0162 [-0.1343, 0.1144] |
| 2 | 1.9367 | 1.9502 | 0.5240 | 0.0735 | +0.0135 [-0.10, 0.15] |
| 3 | 1.8260 | 1.8589 | 0.5243 | 0.0738 | +0.0328 [-0.05, 0.15] |
| 5 | 1.7433 | 1.7950 | 0.5263 | 0.0791 | +0.0517 [-0.03, 0.18] |

This fails the biological admission gate:

- k=1 is worse than scalar level;
- all TERM-vs-cut confidence intervals include zero;
- the k=5 relative MSE gain is about 2.9%, below the 5% importance target;
- full CI/Spearman are below the scalar-level ordering (0.5315/0.0811);
- cyclic label permutation produces virtually identical predictions for k>=2;
- the complete foreign intervention is strongly harmful, but it jointly
  changes protein, support labels and level and therefore cannot isolate
  correct support evidence or ligand--label binding.

The stored run's historical foreign arm also changed the protein and level,
and its `wrong_protein_state` arm replaced only transient TERM evidence.
Neither partial/mixed control is used for the admission decision. The current
code now (i) keeps the recipient protein/query fixed while replacing complete
support ligand/label pairs and recomputing their CIPF/base/evidence, and (ii)
replaces the protein throughout zero-shot, CIPF and TERM. Both corrected
controls require a new evaluation run before quantitative interpretation.

The result establishes failure of the combined learned CIPF+TERM path, but the
unfinished matched A/B/C/D matrix cannot distinguish primitive-identification
failure from router failure or joint-optimization failure. No sealed
confirmation split is opened and no SOTA, excellent-performance, or
important-mechanism-source claim is authorized.

## Corrected implementation rerun

After review, TERM aggregation was corrected from `mean/sqrt(k)` to
`sum/sqrt(k)`, a strict zero-evidence reliability gate was added, and the
foreign/wrong-protein interventions were separated. The same three-seed,
120-step development protocol was rerun from scratch on 2026-08-14.

| k | full MSE | scalar-level / TERM-cut MSE | cut - full (95% component CI) | permuted - full |
|---:|---:|---:|---:|---:|
| 0 | 3.2728 | 3.2728 | 0 | 0 |
| 1 | 2.0341 | 2.0349 | +0.0008 [-0.0069, 0.0105] | 0 |
| 2 | 1.9449 | 1.9486 | +0.0037 [-0.0083, 0.0195] | +0.0004 |
| 3 | 1.8526 | 1.8597 | +0.0071 [-0.0048, 0.0251] | +0.0004 |
| 5 | 1.7858 | 1.7969 | +0.0112 [-0.0037, 0.0346] | approximately 0 |

The corrected candidate still fails admission. Every TERM-vs-cut interval
includes zero; k=1 is practically unchanged; label permutation remains
indistinguishable; and k=5 relative MSE gain is only about 0.62%. Full
CI/Spearman at k=5 are 0.5124/0.0550, only marginally above the level path's
0.5096/0.0517 and not evidence of excellent cold-target ranking.

The corrected foreign-support intervention is highly harmful (about
2.1--2.4 MSE), but it replaces both support ligands and labels and therefore
still does not rescue the missing label-binding result. The full-path
wrong-protein contrast has intervals crossing zero. Gate authorization remains
`G2=false`, `G3a=false`, `G3b=false`; this consumed population remains
development-only.

## Sources

- PSICHIC: https://www.nature.com/articles/s42256-024-00847-1
- Interformer: https://www.nature.com/articles/s41467-024-54440-6
- GroupBind: https://openreview.net/forum?id=zDC3iCBxJb
- modular task generalisation: https://aclanthology.org/2023.eacl-main.49/
- Meta-Point: https://openaccess.thecvf.com/content/CVPR2024/html/Chen_Meta-Point_Learning_and_Refining_for_Category-Agnostic_Pose_Estimation_CVPR_2024_paper.html

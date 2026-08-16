# DCST-R11 affinity-destruction identification preregistration

Date: 2026-07-28  
Status: frozen before implementation and training

## Hypothesis

R11 keeps the R10 PMCE representation, R6 frozen teacher, exact-null matrix,
data, and certificate unchanged. It adds one source-train identification
objective so affinity supervision names the correct target–ligand pairing.

For each source train episode, let `r = y - base` and let `g`, `g_T`, and
`g_L` be PMCE scores for the true input, exact-target derangement, and the
registered within-episode ligand derangement. Define

```text
a(q,r) = cosine(q - mean(q), r - mean(r))
L_destroy = CE(
    [a(g,r), a(g_T,r), a(g_L,r)] / 0.10,
    class=true
).
```

The cosine bounds the objective and prevents increasing a destroyed score
without limit. Temperature `0.10` is inherited from the already registered
structural retrieval loss. `L_destroy` has weight `1.0`, also inherited from
that loss. There is no development tuning.

## Attribution control

PMCE-NoPriv receives the identical true/target-destroyed/ligand-destroyed
inputs, affinity labels, loss, temperature, weight, optimizer, and 4,000
steps. Its frozen upstream representation is the R6 no-privileged state.
Thus any certificate advantage must come from privileged information already
encoded in the pair-specific measure, not from access to destruction
augmentation.

## Frozen source gate

The pre-existing R10 gate remains:

1. privileged R6 segment mechanism pass;
2. at least one PMCE held-source certified band;
3. strictly more privileged than no-privileged certified bands.

The train report must include mean true/target/ligand alignment on the held
source certificate as before. The certificate itself is not relaxed or
trained on source-development labels.

## Conditional Stage 2

Only on source pass may R11 rerun ChEMBL train/development. The certified
source bridge, existing PMCE residual, arms, paired bootstrap, MDE, RMSE,
destruction-removal, negative-transfer, confirmation, and sealed rules remain
unchanged for the first R11 downstream test. This isolates the effect of
source counterfactual identification.


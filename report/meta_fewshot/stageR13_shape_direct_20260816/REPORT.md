# Stage R13: direct interaction-head shape family — gate-blocked at Stage 1

Numerical authority: `RESULT.json` in this directory, produced from the gate
suite `tests/test_shape_direct_synthetic.py`. No real-data training was run:
the family failed its own preregistered structural gate and is recorded, not
advanced.

**Gate count (corrected 2026-08-16).** The suite collects **18** gates:
**16 pass and 2 are recorded as `xfail`**. An earlier revision of this
report said "16 gates" and "15 of 16 algebraic gates pass", which is
inconsistent with itself (15 + 2 > 16) and with the suite. Verify with:

```bash
python -m pytest tests/test_shape_direct_synthetic.py --collect-only -q
```

## The hypothesis

The R9-R12 ladder localized the CI deficit to the supervision leak of the
relative-transport family: it supervises a bilinear potential
`delta(P,i,j)` while the deployed zero-shot ordering uses a *different*
quantity (anchor-mean shape differences). The direct-shape family
(`model/shape_direct.py`, `scripts/train_shape_direct.py`) removes the leak:
`shape(P,L) = s([e(P,L); u(P)]) - mean_m s([anchor_m; u(P)])` — a
protein-conditioned MLP readout of the interaction embedding, anchor-centered
for the no-constant guarantee — and the difference supervision targets
`s(e_i;P) - s(e_j;P)`, the exact deployed ordering quantity. The transport
stays the retained Tanimoto baseline; the closed gate family is not revived.

## Gate outcomes

**16 of 18 gates pass**: exact zero anchor mean, three-branch endpoint
identity, level constancy within a target, protein-blind ligand prior, k=0
endpoint identity, support permutation invariance, query permutation
equivariance, query independence, the label-residual contract, no
query-label input, geometry refused, no dead trainable branch at
k in {0,1,2,5} (four parametrizations), and the private-task abstention gate.

Two synthetic training gates are recorded as expected failures (xfail,
thresholds unmoved):

1. **interaction-branch gate**: three seeds, k=0-only training on the
   held-out bilinear task gives mean full CI 0.60 (gate: > 0.70) with a
   branch gap of 0.14 (gate: > 0.20) — probe diagnostics show the MLP shape
   branch **collapses to near-zero spread under the shape variance term**
   (shape std 0.00-0.06 in the standard recipe; the branch survives only
   with the variance removed, where it still underperforms the bilinear
   readout on this task). The synthetic task favors the bilinear inductive
   bias, but the collapse is a real measured weakness of the MLP readout.
2. **matched-wrong gate**: the k=1 transport (the retained Tanimoto
   baseline) is mildly harmful on the synthetic task (measured gap -0.030:
   correct 1.463 vs wrong 1.433), matching the recorded neutrality/harm
   signature of every query-specific channel tested in this project.

## Verdict

The direct-shape family is **gate-blocked at Stage 1** under its own
preregistered thresholds. It is retained as evidence (module, trainer, gate
suite). Per the protocol the thresholds are not moved.

## Where this leaves the program

R9-R13 complete a closed, evidence-consistent chain: the pair audit localized
the CI deficit; the cliff-weight dose response resolved it; the variance,
trunk-routing, margin-loss and direct-shape hypotheses were each falsified
or gate-blocked with measured causes. The consolidated boundary statement
lives in `report/BOUNDARY_20260816.md`.

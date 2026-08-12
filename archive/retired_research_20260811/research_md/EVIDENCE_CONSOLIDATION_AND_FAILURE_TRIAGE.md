# Evidence consolidation and failure triage

## Identifiability chain

```text
frozen probability-law operator                    PASS
correct-protein geometry and compatibility         PASS
exact-residue pocket localisation                  PASS (development)
teacher ligand-conditioned residue variation       PASS
gauge-free synthetic ordinal trainability          PASS
real ligand-conditioned residue direction          NOT IDENTIFIED
affinity direction                                 NOT TESTED
k<=5 support section                               NOT TESTED
biological z admission                             NOT AUTHORIZED
```

## Current failure

S3R passed every numerical and participation check but failed the first
biological Gate. Candidate gain over chance was `+0.01041`, below `+0.05`, and
the candidate did not separate from B5, foreign ligands, context corruption or
a trained permuted-label learner.

The failure is neither “no biological signal exists” nor “deep learning cannot
solve it.” Phase 2A showed signal in the labels, while S2R showed the estimator
can train. The unresolved bottleneck is the deployed representation tested in
S3R: a global mean of 41-D ligand atom features cannot preserve scaffold
topology or atom-local correspondence. That interpretation is a testable
hypothesis, not a proven causal mechanism.

## Triage

```text
current basis fails R1/R2
  -> stop optimizer retries and capacity scaling
  -> optionally preregister one graph-aware ligand-information audit

graph-aware single-axis audit fails
  -> pose-free sequence+2D residue-direction route lacks evidence
  -> consider a separately governed 3D/pose information route or stop

structural statistic passes independent confirmation
  -> only then register source-affinity correct>ligand and correct>wrong-protein

affinity Gate passes
  -> only then test support-rank sectioning and z admission
```

No raw pair map or variable-length residue score may enter `z`. The frozen
operator `A(F,z)=K(B(z)F(z))` is unchanged and does not certify the upstream
pairwise/AP experiments.

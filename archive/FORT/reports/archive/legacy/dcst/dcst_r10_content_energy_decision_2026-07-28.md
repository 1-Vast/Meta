# DCST-R10 continuous content-energy decision

Date: 2026-07-28  
Decision: `STOP_R10_ADD_HELD_SOURCE_MATCHED_COUNTERFACTUAL_IDENTIFICATION`

## Result

R10 reproduced the R6 privileged segment mechanism but certified `0/4`
PMCE bands versus `1/4` PMCE-NoPriv. No downstream affinity label was loaded.

The privileged candidate's best true band utility was `0.09896`, but the same
band had target-destroyed utility `0.08376` and ligand-destroyed utility
`0.18154`; it therefore failed certification. The no-privileged active band
had true utility `0.15573`, target-destroyed `0.07840`, ligand-destroyed
`-0.06162`, and confidence `0.78705`.

Wall time was `162.402 s`; peak allocated CUDA memory was `533.6 MiB`.

## Diagnosis

Continuous content removed absolute positions and arbitrary roles, but it
exposed a training-certificate mismatch. The source energy matrix is trained
only to rank affinity residuals for true inputs. Target- and ligand-destroyed
scores enter only after training, during certification. Nothing in the
objective requires the correct pair to be more affinity-aligned than either
destruction.

The privileged representation can consequently learn a useful true direction
whose shuffled-ligand version is even more correlated with affinity. This is
not an encoding failure. The next route must identify the true pair against
the same train-split destruction alternatives while retaining held-source
evaluation and an identical no-privileged control.


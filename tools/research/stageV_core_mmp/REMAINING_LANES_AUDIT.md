# Remaining-lanes audit — why the Phase 1 negative branch is complete

Machine-readable: `REMAINING_LANES_AUDIT.json`. This audit does not train, read
new labels or reopen any stage. It checks every legal route that could convert
the Phase 1 bounded negative into an admissible positive claim.

| route | status | decisive evidence |
|---|---|---|
| MSA / coevolution | **BLOCKED ON EXTERNAL ASSET** | M0 plan preregistered but no sidecar exists; no governed UniRef/UniClust snapshot anywhere in the workspace; historical `mmseqs40` directories are tiny local-sequence clustering outputs (largest local sequence JSON = 147 proteins), and network downloads are outside the local tool policy. |
| Functional-site annotations | **FALSIFIED AT IDENTIFIABILITY** | Stage P0: ProteinKG25 GO bags covered 313/387 targets and lost to the calibrated constant (2.27 vs 1.43); no positive ordering signal. |
| Structure/pocket priors | **NOT ADMISSIBLE + LEVEL-REJECTED** | Pocket descriptors are homology-transferred structure priors, not sequence-only inputs; even so Stage H0 rejected them for level (2.4398 vs 2.6179 constant, near shuffled control 2.4941). |
| Davis / KIBA external datasets | **PROMOTION-GATED** | Stage R plan frozen, NOT AUTHORIZED, NOT RUN; training authorized only after a BindingDB promotion that never happened; zero Davis/KIBA labels read. |
| Sealed confirmation split | **SEALED, 0 EVALUATIONS** | Physically isolated `meta_test` surface: 768 cells withheld, 0 evaluations. |
| Looser transformation equivalence classes | **UNREGISTERED; SCREEN-ONLY BY RULE** | Would reintroduce the measured `mu_tau` residual (median 0.269 pK); requires its own preregistration/cancellation analysis and could never satisfy the objective's strict core/context positive gate. |

## Conclusion

With current local assets and standing governance, **no remaining route can
change the Phase 1 bounded-negative verdict into an admissible positive
claim**. The negative branch is therefore complete. Reopening requires one of:
(a) a governed MSA/coevolution asset, (b) a corpus with cross-component
repeated complete transformations, or (c) a governance change authorizing
external datasets before a BindingDB promotion.

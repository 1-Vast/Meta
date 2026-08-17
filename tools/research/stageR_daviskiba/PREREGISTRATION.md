# Davis/KIBA boundary-check design — FROZEN, NOT AUTHORIZED, NOT RUN

Date: 2026-08-18 (night). Status: this document is a frozen plan only. The
standing governance authorizes exactly one dataset lineage for training
(the BindingDB-Ki double-cold protocol); Davis/KIBA training is authorized
only after a candidate passes the promotion gates. No candidate passed, so
nothing here has been executed and no Davis/KIBA label has been read.

## Purpose

Test the scope of the boundary (report/BOUNDARY_20260817_NIGHT.md) on the
sealed Davis assets: does the level wall reproduce when the protocol is
re-applied to a different public dataset? The Davis assets already carry
rich assay covariates (endpoint vocabulary Kd/Ki/IC50/functional-response/
score-like, DAVIS_mechanism_v2 manifest), which makes them the natural
independent test of the assay-history anatomy (D0_LEVEL_ANATOMY, D0b).

## Available sealed assets (inventory only, no labels read)

- dataset/sealed/DAVIS_mechanism_v2: source/metaval/recipient rows with
  ligand banks and biological context (63 metaval targets; endpoint
  vocabulary of 6 assay types).
- dataset/sealed/DAVIS_homology_v1: metaval-only sealed evaluation set.
- dataset/sealed/DAVIS_v1: carries DO_NOT_USE_FOR_STRICT_EVALUATION.txt
  (legacy) - excluded by design.
- KIBA: no local assets; acquisition and governance are prerequisites and
  are not planned here.

## What would be run (whenever authorized)

1. Split audit: verify the Davis mechanism split's component/document/
   ligand closures and seal state before anything else.
2. D0 anatomy replication: level decomposition per episode; within-
   document transfer R^2; endpoint/assay-covariate identifiability probes
   with component-fold selection.
3. Baseline replication: the T2 recipe (similarity_only trunk, Stage B
   loss, leak-free internal checkpoint selection) retrained from scratch on
   Davis source cells, evaluated on the frozen Davis metaval banks
   (k = 0/1/2/3/5), with the full counterfactual set.
4. Strongest-mechanism replication: K-REG (contrastive coembedding with
   regression alignment) under the same rules.
5. Gates mirror the BindingDB stages: G1 k=0 not degraded; G2 MSE/ranking
   joint; G3 no resolved ranking degradation; G4 controls; G5 cost. The
   boundary check itself passes if the level wall and the assay-history
   anatomy reproduce qualitatively (level share of k=0 error >= 50%,
   within-document transfer substantially above cross-document transfer).

## Governance

No Davis/KIBA label, row or statistic has been read for this plan; the
inventory above uses directory listings and the mechanism-v2 manifest's
schema/count fields only. Execution requires a written authorization
recorded in the stage artifact, exactly as the meta_test seal requires.

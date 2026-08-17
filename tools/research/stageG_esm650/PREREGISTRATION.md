# Stage G preregistration — ESM-650M residue-input trunk (external representation lane)

Frozen before the arm trained. Date: 2026-08-17 (night). Single seed,
development evidence; meta_test sealed and never constructed.

## Question

Does replacing the governed ESM-150M protein bank (640-dim pooled +
128-slot residues) with the local ESM-2 650M embeddings (1280-dim pooled +
128-slot residues, recorded with model/weights hashes in
tools/runtime/esm2_t33_650M_pooled/manifest.json) improve the incumbent
similarity_only recipe on the frozen meta_val banks? The D0 probes showed
ESM-650M pooled carries marginally more cross-component level signal than
150M (linear level MSE 1.6875 vs 1.7478); the residue-slot inputs are
untested end-to-end.

## Arm

G = the Stage D T2 arm byte-identical (1,200 steps, seed 20260815, Stage B
loss recipe, leak-free internal checkpoint selection, GPU verification) with
only the protein bank swapped. Baseline: the frozen T2 checkpoint.

## Gates

G1. k=0 MSE not degraded beyond +2% and its interval not resolved positive.
G2. At least two of k in {2,3,5} improved in MSE with resolved intervals.
G3. Spearman/CI not degraded at any k by a resolved interval.
G4. Correct-support dependence (permuted/matched-wrong above correct,
    resolved) preserved.
G5. Cost: parameters identical to T2; wall time and peak VRAM <= 1.5x T2.

Stop rules: S1 G1/G2 fail; S2 G3 fails; S3 any control inverts.

## Disclosure

The 650M bank is external data (a local model snapshot); if this lane were
ever promoted, the production pipeline must re-record the bank provenance
and the double-cold split remains the only authorized protocol. Davis/KIBA
would still be trained independently from scratch.

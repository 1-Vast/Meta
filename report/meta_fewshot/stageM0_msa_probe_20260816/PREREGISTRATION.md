# Stage M0 preregistration: MSA protein-side prior — diagnostic probes, no model training

Status: **preregistered, not run** (2026-08-16). This stage is the governed
entry point for the MSA/evolutionary-prior direction ratified in the
2026-08-16 adjudication. It produces diagnostics and named baselines only.
No production model is trained, no checkpoint is modified, `meta_test` is not
read and no MSA is built for a `meta_test` target. A small, fixed-budget
diagnostic probe may be optimized on `meta_train` only; it is not a candidate
model and cannot be promoted to deployment.

## Standing and non-goals

Established by adjudication and binding on this stage:

1. **No structural machinery.** 0/17,717 deployment cells have a common-frame
   complex (`task.md`, geometry is settled). Evoformer-style MSA-pair-structure
   exchange, Cartesian/equivariant interaction encoders and conformational
   routing stay closed. MSA enters, at most, as sequence-derived protein-side
   features.
2. **No target-similarity support transport.** Under the episode contract
   (`scripts/qpsmp_data.py::draw_episode`) support and query are ligands of
   one recipient target, so a protein-similarity support weight is constant
   within an episode. Any cross-target (donor) support arm is a protocol
   amendment and is not proposed here.
3. **No CD-HIT40-granularity family affinity prior.** The double-cold split
   guarantees zero protein-component overlap with training, so a component-
   level family mean is undefined for every evaluation target. A *coarse*
   family prior (below CD-HIT40 granularity, train-only construction) is
   permitted **as a named baseline only** (M0-B), never as a core mechanism.
4. **MSA is not a core-innovation candidate.** The R5 mandate allows at most
   two core innovations and requires the training innovation to be one of
   them; R14 (`stageR14_diagnostics_20260816/REPORT.md`) localizes the
   surviving lever to the within-target ordering coefficient `r`. M0 tests
   only the *calibration* half of `shape = Var(y)(1-r²) + amplitude`, i.e.
   whether MSA carries target-level information that the frozen ESM
   representation misses. Conservation is a **protein evolutionary prior**,
   not a binding-pocket map; no contact or pocket claim is licensed by
   anything in this stage.
5. **Retrieval remains a named baseline** (R0 falsification). M0-B is the
   MSA-specific instance of the protein-retrieval test, not a reopening of
   retrieval as a mechanism.

## D0 (prerequisite): MSA sidecar construction

- **Population:** the 346 `meta_train` targets and 41 `meta_val` targets of
  `bindingdb_ki_double_cold_v1` only. Sequences are read from the existing
  protein bank. `meta_test` targets are never queried.
- **Search:** MMseqs2 `easy-search` (or jackhmmer, if MMseqs2 is unavailable
  in the `drug` environment) against a local UniRef90 snapshot, E-value 1e-3,
  coverage ≥ 0.5. Tool version, database version and all parameters are
  recorded in the sidecar manifest. A target with fewer than 10 hits after
  filtering is recorded as `depth=0` and retained (depth is itself a
  stratification variable in M0-C, not an exclusion).
- **Hard prerequisite:** D0 cannot start until the executable, database path,
  database version/SHA256, query count, hit count and failed-query count are
  recorded. A historical `mmseqs40` directory is not a UniRef database.
- **Pooled target-level features (exactly eight, fixed before any probe):**
  1. `log_neff` — log10 of the effective alignment depth (N_eff);
  2. `gap_fraction` — mean gap/insertion rate over query positions;
  3. `mean_conservation` — mean per-position relative entropy against UniRef
     background frequencies;
  4. `top_decile_conservation` — mean conservation of the top 10% positions;
  5. `max_apc_mi` — maximum APC-corrected mutual information over residue
     pairs;
  6. `top1pct_apc_mi` — mean of the top 1% APC-corrected MI values;
  7. `mean_identity_top50` — mean sequence identity to the 50 closest hits;
  8. `length` — protein length.
- **Artifact:** `stageM0_msa_probe_20260816/sidecar/msa_features.parquet` plus
  `manifest.json`. Per-position vectors (conservation profile, APC-MI matrix
  row maxima) are additionally cached for a possible M1, but M0 itself
  consumes only the eight scalars.

## M0-A: frozen-representation supervised diagnostic

**Question.** Does MSA explain part of the incumbent's per-target k=0
calibration error that the frozen ESM representation does not?

- **Response.** Per-target residual `r_t = mean_cells(y − ŷ_k=0)` of each
  retained A0 seed separately. Do not average checkpoints before fitting or
  call the resulting ensemble the incumbent. Report the three paired seed
  deltas and their median; an ensemble is a separate named baseline.
- **Fitting.** Use a fixed low-capacity `torch.nn.Linear` diagnostic probe,
  standardized with `meta_train` statistics and optimized for a fixed 256
  AdamW steps inside leave-one-component-out folds on `meta_train`. No Ridge,
  matrix solve, pseudoinverse or closed-form adaptation is allowed, even in
  the diagnostic implementation. The probe is not a production model.
  One final evaluation is made on the 41 `meta_val` targets.
  Three feature arms: `ESM` (pooled ESM-2 embedding, the incumbent's own
  protein representation), `MSA` (the eight scalars), `ESM+MSA`.
- **Primary metric.** `meta_val` residual MSE against the constant-residual
  baseline (predicting the train residual mean), equal-target weighting;
  component-level bootstrap over the 19 `meta_val` components, 90% interval,
  10,000 draws. Reported as `Δ(MSE) = MSE(ESM+MSA) − MSE(ESM)` with interval.
- **Controls (all binding):**
  - *Residual permutation:* refit the `ESM+MSA` arm with train residuals
    shuffled across targets (100 draws). The val increment must collapse
    (≥80% of the unshuffled |Δ| destroyed); otherwise the increment is
    target-identity leakage, not signal.
  - *Depth partialling:* refit `MSA` with `log_neff` and `gap_fraction`
    removed. The `ESM+MSA` increment that survives this arm is
    non-depth signal; the difference is attributed to depth in M0-C.
  - *ESM sanity:* the `ESM` arm must beat the constant baseline on train
    LOCO before its val number is interpretable.
- **Preregistered gates:**
  - **PASS → M1 permitted:** Δ(MSE) point estimate ≤ −0.01 pK² *and* the
    component-bootstrap 90% upper bound < 0 *and* the permutation control
    destroys ≥80% of the increment.
  - **AMBIGUOUS → park:** point estimate ≤ 0 with interval crossing zero.
    No M1. Recorded as unresolved; the direction may be revisited only with
    new evidence.
  - **FAIL → close:** point estimate > 0. The MSA mainline is closed for
    the double-cold calibration target; the sidecar remains as recorded
    covariates.
- **Non-claims.** M0-A is a diagnostic of a frozen representation. It is
  not a DTA performance claim, not an input to any checkpoint, and confers
  no protein-conditioned language on any retrieval variant.

## M0-B: no-training protein-kernel and coarse-family baselines

The incumbent has **no active protein kernel** (its transport is
ligand-Tanimoto), so nothing is "replaced"; these are standalone named
baselines for target-level label-mean prediction on `meta_val`:

- `global_mean` — mean of train-target label means;
- `esm_kernel` — Gaussian kernel on train-centred pooled ESM embeddings;
- `msa_kernel` — Gaussian kernel on standardized MSA scalars;
- `esm_msa_kernel` — product of the two. Kernel bandwidth is selected only
  inside `meta_train` by fixed leave-one-component-out validation; no
  `meta_val` tuning is permitted;
- `coarse_family` — CD-HIT 60% and 80% clusterings built **on `meta_train`
  targets only**; a val target inherits the mean of its cluster's train
  label means, falling back to `global_mean` when its cluster is empty;
  family overlap fractions are reported per granularity;
- each kernel's **shuffled-target control** (train-target label means
  permuted), per the R0 finding that protein-conditioned retrieval lost to
  its own shuffle.

Metric: per-target label-mean MSE on `meta_val`, equal-target weighting,
component bootstrap 90% intervals. **Gate:** a kernel counts as carrying
protein-side information only if it beats both `global_mean` and its own
shuffled control with the interval excluding zero. Failure of every
protein-side baseline is a recorded strengthening of the R0 conclusion, not
a null result.

## M0-C: stratification and confound audit

All M0-A/B deltas are additionally reported by stratum with bootstrap
intervals: `log_neff` tertiles, `gap_fraction` tertiles,
`mean_identity_top50` tertiles, protein length tertiles, labels-per-target
tertiles, ligand novelty tier (exact-ligand-present vs exact-free, the R0
48.9% overlap split), and coarse-family overlap.

**Confound rule (preregistered):** if the M0-A increment is positive in the
top depth tertile and non-positive in the bottom tertile, the increment is
flagged **depth-confounded** (study-bias proxy) and M1 is *not* permitted on
that evidence alone; the park branch applies.

## Decision tree and successors (previews, each with its own preregistration)

- **M1 (only on M0-A PASS, unconfounded):** one local attention block over
  the existing protein slots with conservation/APC-MI scalar bias, trained
  **MSE-primary** (R14 Result 1: shape-first training trades `r` away; the
  G1 attribution trap), k=0 protein representation only, transport
  untouched, three seeds, matched budget with an A0 retrain control. Gates:
  k=0 MSE component-bootstrap vs A0, wrong-protein gap > 0 at k=0, no CI
  regression.
- **M2 (high-threshold candidate, low prior):** support-conditioned protein
  evolutionary-state routing. Enters a family with seven falsified
  query-specific channels; its only differentiating claim is that state
  weights derive from external evolutionary structure rather than free
  learned gates. Admission requires all six adjudicated conditions: k=1
  beats the scalar level shift; k=0 strictly non-regressing; support-label
  permutation destroys the effect; matched-wrong support significantly
  worse; state weights not reducible to ligand-Tanimoto; gain attributed to
  state selection, not calibration shift.
- **Standing constraint:** pairwise/cliff-aware counterfactual training
  remains the core-innovation track; MSA may not displace it.

## Records

Artifacts land in this directory: `msa_residual_diagnostic_meta_val.json`
(M0-A), `protein_kernel_baselines_meta_val.json` (M0-B),
`stratification_meta_val.json` (M0-C), `REPORT.md` and `RESULT.json` on
completion. No `RESULT.json` exists before the stage runs, so the record
audit's classification is untouched until then.

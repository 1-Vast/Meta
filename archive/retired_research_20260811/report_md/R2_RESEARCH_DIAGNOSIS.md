# R2 — Identifiability-regime diagnosis and repair for MetaSieve partner specificity

> **SUPERSEDED 2026-08-11.** The headline below is not accepted.  The registered
> E1 result (`gauge_ratio>1` on both consumed splits) falsified the proposed
> H0-regime; positive ridge and the learned population-coordinate term also
> invalidate the report's complete-predictor GL claim.  The corrected terminal
> diagnosis is [R2_MULTI_AGENT_RESOLUTION.md](R2_MULTI_AGENT_RESOLUTION.md).
> This file remains unchanged below this notice as an auditable Cowork record.

Date: 2026-08-11. Research report. **No Gate is opened by this document.**
Preregistration: `research/meta_fewshot/PREREG_R2_IDENTIFIABILITY_REGIME.md`.

```text
HEADLINE   IDENTIFIABILITY_REGIME_DEFECT_IDENTIFIED_REPAIR_DERIVED_EXECUTION_BLOCKED
```

---

## A. Research diagnosis

### A.1 What MetaSieve v0 actually learns

Reading `train_main_v0.py::MetaSectionRegressor`, the deployed predictor is

```
pop(L)      = Linear(fingerprint(L))                     # ligand only
m(P,L)      = B^T phi_288(P,L),   B in R^{288 x d}, B^T B = I
rho_i       = y_i - pop(L_i) - m_i . c_pop
w           = M^T (M M^T + lam I_k)^{-1} rho             # = (M^T M + lam I_d)^{-1} M^T rho
yhat_q      = pop(L_q) + m_q . (c_pop + w)
```

with the deployed operating point `d = 2`, `k = 5`, `lam = 1.0`.

**Proposition 1 (gauge invariance of the estimator).** Let `G ∈ GL(d)` act
episode-wise as `M → MG`, `m_q → m_q G`. If `rank(M) = d` (guaranteed
generically for `k ≥ d`) then at `lam → 0`

```
yhat_q  =  pop(L_q) + m_q (M^T M)^{-1} M^T rho
```

is **exactly** invariant under `G`. At `lam > 0` invariance is exact for
`G ∈ O(d)` and is broken only at `O(lam / sigma_min(M^T M))`.

*Proof.* `m_q G (G^T M^T M G)^{-1} G^T M^T rho = m_q (M^T M)^{-1} M^T rho`. For
orthogonal `G`, `G^T M^T M G + lam I = G^T (M^T M + lam I) G`. ∎

**Corollary 1 (the sufficient statistic).** All protein information reaching the
prediction is carried by the support-weight vector

```
u = m_q (M^T M + lam I)^{-1} M^T  in  R^k ,      yhat_q = pop(L_q) + u . rho
```

modulo an `O(lam)` shrinkage-toward-`c_pop` term. `u` lies in the `d`-dimensional
column space of `M`. **The biological frontend's entire causal role in v0 is to
choose how to weight five support labels.** It never contributes an additive
prediction of its own.

This is why replacing the protein on both sides recovers the correct-arm error.
Wrong/wrong recovery is a **theorem about the estimator**, not a fact about
biology, and it does not require any exact global orthogonal gauge. A0 tested
`O(288)` alignment of the 288D features and correctly found none; the invariance
that actually matters is the far weaker episode-wise `GL(2)` acting on two
coordinates. **A0 tested the wrong invariance group**, which is why its residuals
(0.445, 0.478, 0.729) look large while wrong/wrong nevertheless recovers.

### A.2 The regime error — biology is redundant *by theorem*

`FINAL_THEORY_TO_MODEL_HANDOFF.md` §6 states that the support supplies **at most
`k` continuous dimensions of member identity** (F20/CP-3), and §12 CI-A3 states
that auxiliary information is useful at `(D,x)` **iff it changes the joint
window** `T_{D∪{x}}`.

**Corollary 2 (regime redundancy).** If the model family restricted to the
covered region is exactly `d`-dimensional and the support design of size `k ≥ d`
is unisolvent for it, then `ω_{x,D}(0) = 0`, the section at `ε = 0` is a
singleton, and no fiber restriction can shrink it. Hence the certified auxiliary
information `Γ_c = ½(ω(2ε) − ω(2ε | c_b))` tends to 0 as `ε → 0`.

MetaSieve v0 runs at `d = 2`, `k = 5`. It is **2.5× over-supported**. Its own
frozen theory therefore predicts that the protein auxiliary contributes nothing
at first order, and that whatever it does contribute is an `O(ε²)` variance
effect. CI-A5's "harmless-for-all-data iff `T(c') = T(c_b)`" is satisfied
trivially, because the window is pinned by the support alone.

**MetaSieve selected an operating point at which its own theory forbids biology
from mattering, then measured that biology does not matter.**

### A.3 The published numbers already say this

From `MAIN_V0_REPORT.md`, target-macro MSE:

| Arm | MSE | Share of the total gain over `d=0` |
|---|---:|---:|
| Population `d=0` | 8.711 | — |
| **Full permuted** | **2.047** | **98.1 %** |
| Full correct | 1.916 | 100 % |
| Full wrong protein (= wrong/wrong) | 2.113 | 97.1 % |

Support-label permutation is a cyclic roll (`source_y.roll(1)`); it preserves
the label multiset exactly and destroys only the ligand↔label correspondence.
Therefore:

- **98.1 % of the v0 Meta-Section gain survives the destruction of all
  ligand-specific information.** The mechanism is a `k`-shot per-target
  intercept.
- Only **1.9 %** of the gain is ligand-specific, and only **2.9 %** is
  partner-specific.
- Within-target discrimination confirms it: Pearson is 0.082 at `d=0`, 0.097
  correct, and **0.103 with the wrong protein**. The section adds no reliable
  ranking signal, and the wrong protein ranks better than the right one.

Against the §12 acceptance criteria the published v0 result therefore already
fails *Ligand-specific prediction*, not only *Partner specificity*.

**The single most important missing control is absent from the v0 battery: the
one-parameter per-target intercept null** `yhat_q = pop(L_q) + mean_i(rho_i)`.
Because the permuted arm is a noisier version of that null and still reaches
98.1 %, the intercept null may well *beat* `full_correct`. Until it is run, the
v0 "REAL_BIOLOGICAL_META_SECTION_V0_PASS" verdict is unsupported. This is E0.

### A.4 Why T-BASIS fails — an independent, upstream defect

From `run_tbasis_radial.py` and `generate_tbasis_features.py`:

```
T[a,r,k] = (1/n_atom) * sum_{i,s} atom[i,a] * comp[s,r] * radial[i,s,k]
           8 atom channels x 6 residue CLASSES x 6 radial shells = 288
comp[s,r] = fraction of residue class r inside positional slot s of the sequence
```

The protein enters **only** through `comp[s,r]` — a 6-way chemical-class
composition profile over sequential positional slots — and through the predicted
contact weights, against which that profile is contracted and the residue axis
summed away. Consequences:

- **No residue identity.** Thr/Ser and Leu/Ile/Val are indistinguishable. Kinase
  gatekeeper chemistry, hinge motifs and the specific substitutions that create
  selectivity are largely invisible.
- **No pocket localisation.** Slots are sequence-index bins over the whole chain,
  not a binding site.
- **No conservation, no geometry, no structure.**
- Homologs at 40 % identity have nearly identical 6-class slot profiles.

This predicts A1's result exactly, and explains the otherwise puzzling A1 row in
which the *nuisance length/composition* baseline (0.839) scored better than
calibrated T-BASIS (0.926): **T-BASIS is, to first order, a composition
descriptor.** E2 measures this label-free, so it does not reopen A2.

### A.5 What the new solution must target

Two independent defects that must both be repaired, in this order:

1. **Regime (necessary).** Make the `k`-shot support *provably insufficient*, so
   the auxiliary is mathematically required. No representation repair can help
   while Proposition 1 holds — a perfect interaction descriptor would still be
   gauged away.
2. **Resolution (necessary).** Give the partner channel enough resolution to
   separate homologs — residue identity and pocket localisation.

Fixing 2 without 1 was V1. Fixing 1 without 2 gives a necessary but empty
channel. Both are required, and 1 is logically prior.

---

## B. Existing-data exploitation map

The A0/A1 report concluded that a dense assay-matched crossed cohort must be
*acquired*. Inspection of `dataset/processed/` shows a substantial part of it is
**already on disk**.

| Artifact | Content | Best scientific role |
|---|---|---|
| `crossed_panels_xp2/blk_metz_xp2.npz` | Metz kinase panel: 928 compounds × 147 kinases, **32,849 measured cells, density 0.241**, 258 scaffold components, KLIFS-mapped, 50,313 censored cells | **Primary ΔΔ interaction pre-training.** Interaction df ≤ 32,849 − (147+928−1) = **31,775**. Only 8 declared groups ⇒ training supply, *not* confirmation supply. |
| `multipanel/blk_bdb_panels.npz` | 85 assay-matched BindingDB panels, 6,363 cells, 129 targets, 2,845 ligands, **70 mmseqs40 clusters**, 179 scaffold components | **Confirmation supply for E5.** Multi-family, so it tests cross-dataset generality. Interaction df ≤ 3,390 (one-block bound). |
| `crossed_panels/pdsp_core.npz` + `esm2_t30_pdsp.npz` | PDSP Ki panel — GPCR/transporter, non-kinase | **Second, non-kinase confirmation supply.** Directly answers requirement 5. |
| `crossed_panels/esm2_t30_kinase_pocket85.npz` | ESM-2 embeddings restricted to the **85-residue KLIFS pocket** | **The repair for defect A.4.** Residue-resolved, pocket-localised partner descriptor. Already computed. |
| `crossed_panels/klifs_structures.json` | KLIFS pocket residue mapping | Pocket alignment across kinases; enables residue-indexed partner coordinates. |
| `crossed_panels/metz{60,70}_conformation_features.npz` | Conformational state features | DFG-in/out conditioning; secondary. |
| `meta_fewshot/bindingdb_ki_main_v0` (21,473 cells) | Sparse, heterogeneous BindingDB Ki | **Episodic Meta-Section training/eval only.** Its 1,820 selectivity groups closed into 21 components with 86.4 % in one — structurally unsuited to the A1 question. |
| `ssl_b2/teacher_dataset.npz`, `correspondence_router/c0_geometry.npz` | Contact/distance teacher, geometry | Geometry supervision for the interaction head; use only under a separate Gate. |
| `dta/kiba.tab`, `cache/protein_DAVIS.pt` | Davis / KIBA | External confirmation only. Do **not** merge endpoints into Ki regression. |
| `source_affinity/chembl37_sqlite_v1` | ChEMBL 37 | Assay-block/document context; matched molecular series; reserve. |
| `protein_benchmarks/LMDB/Contact.tar.gz`, `ProteinGym.tar.gz` | Contact maps, mutational scans | Representation pre-training only; ProteinGym gives *local mutation → function* supervision relevant to §6C. |

**The decisive structural point.** A1 asked whether partner selectivity is
recoverable, using a substrate with 21 dependency components and one component
holding 86.4 % of the mass. Metz supplies ~31,775 interaction degrees of freedom
under a single uniform assay protocol. **A1's negative result was obtained on the
wrong data for the question it asked**, and cannot be read as evidence that the
information is absent.

Honest counterweight, and the reason E3 is preregistered with two thresholds:
Metz's manifest declares only **8** dependency groups. Interaction df measures
*training* supply; dependency components bound *confirmatory* power. Metz is
ample for the former and inadequate for the latter. Confirmation must come from
BDB-panels (70 clusters) and PDSP.

---

## C. Literature synthesis

| Work | Mechanism | Supervision | What to borrow | What **not** to borrow |
|---|---|---|---|---|
| [ADKF-IFT (Chen et al., ICLR 2023)](https://arxiv.org/abs/2205.02708) | Bilevel deep-kernel GP; features meta-learned, GP hyperparameters task-fitted, solved by implicit differentiation | Molecular property, task = assay | The *split* between meta-learned and task-fitted parameters, and the principled inner solve | It adapts a kernel, not a reserved subspace; nothing forces a protein descriptor to be used. Adopting it wholesale reproduces the redundancy. |
| [MBP (Yan et al., Brief. Bioinform. 2024)](https://academic.oup.com/bib/article/25/1/bbad451/7469349) | Multi-task bioassay pre-training on ChEMBL-Dock; classifies **relative rankings within the same bioassay** | 300k affinities, 2.8M docked poses | The within-assay relative estimand — the closest published precedent to MetaSieve's ΔΔ. Removes assay offsets exactly as double-centring does. | Docked 3D poses (MetaSieve is sequence-only here); and MBP ranks *within one protein*, so it cancels the assay effect but **not** the protein main effect. MetaSieve needs the crossed, two-sided version. |
| [PSICHIC (Koh et al., Nat. Mach. Intell. 2024)](https://www.nature.com/articles/s42256-024-00847-1) | Physicochemical-constrained GNN producing residue-level interaction fingerprints from sequence alone | Paired affinity, no structure | Residue-resolved interaction fingerprints, sequence-only, with demonstrated **selectivity determinants** — the existence proof that A.4 is fixable without structure | It is a single-task predictor with no support-conditioned adaptation; importing it as an encoder without fixing the regime changes nothing. |
| [MetaDTA (ICLR 2022 workshop)](https://openreview.net/forum?id=yzlif16IASM) | Attentive Neural Process; support set defines the per-target regression function | BindingDB/Davis | The ANP episodic framing | Directly relevant caution: independent evaluations report MetaDTA performance tracks *task similarity* (r ≈ 0.49) and that plain **few-shot per-target calibration outperforms MAML/ProtoNet baselines** — i.e. the exact failure mode diagnosed in A.3. Do not adopt an architecture whose reported gains are calibration. |
| [AdaMBind (Nat. Commun. 2026)](https://www.nature.com/articles/s41467-026-70554-5) | MAML + adaptive task sampling weighted by query loss and support/query gradient similarity + label-noise strategy | Graph + sequence, three benchmarks | Task-value weighting and the label-noise treatment | Gradient-similarity task weighting is an *outcome-dependent* sampler — a leakage hazard under MetaSieve's discipline. MAML inner loops adapt everything, maximising redundancy. |
| [CAVIA (Zintgraf et al., ICML 2019)](https://proceedings.mlr.press/v97/zintgraf19a/zintgraf19a.pdf) / ANIL | Partition parameters into shared and *context* parameters; adapt only context | Regression/classification/RL | The formal precedent for **partial adaptation** — the nearest prior art to RFMS | Context parameters are still fitted from the support. RFMS's reserved block is fitted from the **auxiliary**, and the split is chosen by an identifiability argument (`d_S ≤ k`, `d_c ≥ 1`), not by parameter count. |
| Direct selectivity regression (kinase literature, e.g. [JCIM 2023](https://pubs.acs.org/doi/10.1021/acs.jcim.3c00347), [J. Med. Chem. 2017](https://pubs.acs.org/doi/abs/10.1021/acs.jmedchem.6b01611)) | Train directly on the affinity **difference** between two targets rather than on two absolute models | Kinase profiling panels | Confirms that `Δy` is a better-conditioned estimand than paired absolute predictions; validates §6B | Kinase-only; no meta-learning; no held-out-family discipline. Cannot be cited as evidence of cross-family transfer. |

**Novelty placement.** Partial adaptation exists (CAVIA/ANIL). Within-assay
relative supervision exists (MBP). Direct selectivity regression exists. Deep
kernels with a meta/task parameter split exist (ADKF-IFT). What does **not**
exist in this literature is a meta-training principle in which the adaptation
operator is *deliberately rank-deficient relative to the family dimension, with
the deficit assigned to the auxiliary channel and certified by the theory's own
`Γ_c`*. The novelty is **statistical/identifiability-level**, not architectural.
It should be claimed as such and no further.

---

## D. Candidate training mechanisms

### M1 — Reserved-Fiber Meta-Section (RFMS) — *derived, not proposed*

Split `R^d = V_S ⊕ V_C`, `dim V_S = d_S ≤ k`, `dim V_C = d_c ≥ 1`. The ligand
frame `ψ(L) ∈ R^d` takes **no protein argument**. The partner map `c0(P)` writes
only into `V_C`. The support cuts only `V_S`.

```
yhat(L) = psi(L) . ( c0(P) + Pi_S delta ),   delta = ridge fit on support
rho     = || y_s - Psi_s (c0 + Pi_S delta) ||_inf        (CI-A5 empty-fiber statistic)
```

Why partner specificity becomes **necessary rather than auxiliary**:

- there is no protein-dependent basis, so Proposition 1's group is trivial —
  substituting a wrong protein cannot be cancelled on the query side;
- `V_C` is outside the range of the support solve, so it is not re-fittable at
  any `k`;
- `Γ_c > 0` identically whenever `c0` varies across targets and `V_C` carries
  affinity variance. Wrong/wrong recovery becomes **mathematically impossible**;
  the size of the break *measures* `Γ_c` rather than merely passing a Gate.

Theory compliance: `d ≤ 5`; support cut `≤ k` dims; `c0 ∈ R^{d_c}`, `d_c ≤ 4` —
this is **not** a high-dimensional target embedding; permutation symmetry and
affine equivariance are preserved; `ρ` supplies the realizability flag the
handoff requires (§7 O4, tests T1/T3), which v0 does not emit at all. `c0` is
also the first legitimate candidate for a low-dimensional biological `z`.

### M2 — Crossed-interaction (ΔΔ) meta-pretraining on matched panels

Change the **information source**, not the network. Pre-train the interaction
head on real double differences inside assay-matched blocks,

```
DeltaDelta y = y_ai - y_aj - y_bi + y_bj
```

which removes protein and ligand main effects and the assay offset exactly, then
hand the resulting coordinate to RFMS. Feasibility rests on E3; the Metz
manifest alone implies ~3.2×10⁴ interaction df. Statistical care required:
censored cells (50,313 in Metz) must be handled as interval-censored, not
dropped and not imputed; unequal replication must be weighted; and ΔΔ variance
inflates measurement noise 4×, so the estimand needs the dense panel, not the
sparse BindingDB corpus.

### M3 — Pocket-resolved partner coordinate

Replace the 6-class positional-slot contraction with a residue-identity-resolved,
pocket-localised descriptor, using `esm2_t30_kinase_pocket85.npz` and the KLIFS
mapping already on disk, generalised to non-kinases via a pocket definition
derived from the local interaction teacher. This is representation repair; it is
*required* for M1 to have anything to say, and *insufficient* on its own.

### Ranking

| | Novelty | Bio plausibility | Theory compat | Data feasibility | Expected effect | Cost | Shortcut risk |
|---|---|---|---|---|---|---|---|
| **M1 RFMS** | High (identifiability-level) | Neutral | **Derived from §6/§12** | High — needs only a protein descriptor | Large on the *contrast*; unknown on absolute MSE | Low | Low: wrong/wrong break is structural, so it cannot be faked — but see E.2 |
| **M2 ΔΔ** | Moderate (MBP is close) | High | Compatible; changes the estimand not the operator | **Pending E3**; Metz suggests ample | Moderate | Medium | Medium: panel identity and scaffold series can leak |
| **M3 pocket** | Low (PSICHIC-adjacent) | High | Neutral | High — artifacts exist | Enabling, not sufficient | Medium | Medium: family memorisation |

Recommended order: **M1 → M3 → M2**, with M1's E0/E1/E2 diagnostics first
because they are label-free or descriptive and decide whether the other two are
worth funding.

---

## E. Adversarial review of M1

**E.1 "RFMS is just low-rank matrix factorisation with side information."**
Substantially true as an architecture — with `ψ` ligand-only and `c0` predicted
from sequence, the model is a `d`-rank factorisation of the protein×ligand matrix
whose protein factor is regressed from the partner. The defensible novelty is
*not* the factorisation; it is the **reserved block plus the identifiability
budget** `d_S ≤ k < d`, and the use of `Γ_c` and the empty-fiber residual as the
reported quantities. Claiming architectural novelty here would be wrong.

**E.2 "The wrong/wrong Gate now passes vacuously."** The strongest objection. If
`c0` collapses toward a constant, the Gate passes because the arms differ by a
fixed offset, not because biology was learned. Mitigation, preregistered:
report `c0_between_target_variance` and treat collapse as a **FAIL**; and require
the *ligand-specific* Gate (within-target Pearson/Spearman/CI) in addition to
MSE, which a constant `c0` cannot pass.

**E.3 "Reserving `V_C` will simply make absolute performance worse."** Likely,
and it must be reported. RFMS deliberately forfeits `d_c` dimensions of support
fit. If the biological channel is empty, RFMS is strictly worse than v0 — that is
the point: it converts a silent failure into a loud one. A worse MSE with a
large, cluster-macro-stable wrong/wrong break is scientifically more informative
than v0's better MSE with no partner content. It is **not** grounds for a
production migration.

**E.4 "Proposition 1 is defeated by `lam = 1.0`."** The ridge does break exact
`GL(d)` invariance, and A0's non-trivial Procrustes residuals show the correct
and wrong coordinate clouds are not literally aligned. But invariance under
`O(d)` is exact even with ridge, and the empirical wrong/wrong recovery
(`−0.033`, LCB `−0.090` cluster-macro in V1) is the direct measurement that the
residual non-invariance is negligible in practice. E1 settles this label-free by
measuring `u` directly rather than the coordinates.

**E.5 "The intercept-null claim over-reads the permuted arm."** Fair. A rolled
support label vector is not an intercept estimator; it is a linear fit to shuffled
targets, which is *noisier* than the mean. That makes the inference conservative
in one direction only: the true intercept null should be at least as good as
2.047, i.e. it may beat `full_correct` at 1.916. E0 must be run before the claim
is asserted as fact rather than as a strong prior.

**E.6 "This is one more architecture proposal."** The distinguishing test is that
M1's central claim is falsifiable *without any training*: E1 either shows `u` is
near-uniform and gauge-like, or it does not.

---

## F. Minimal implementation

Placed under `research/`; **no production `model/` or `scripts/` file touched.**

| File | Role | Cost | Labels opened |
|---|---|---|---|
| `research/meta_fewshot/PREREG_R2_IDENTIFIABILITY_REGIME.md` | Preregistration, decision rules, falsifiers, fresh-supply definition | — | — |
| `research/meta_fewshot/r2_e0_regime_audit.py` | **E0** intercept null + **E1** label-free sufficient-statistic audit, on existing v0 checkpoints | minutes | E0 only, already-consumed, descriptive |
| `research/meta_fewshot/r2_e2_tbasis_decomposition.py` | **E2** two-way ANOVA of T-BASIS, fixed-ligand partner dispersion, homolog resolution, corruption calibration | minutes | **none** |
| `research/meta_fewshot/r2_e3_crossed_census.py` | **E3** interaction df vs dependency components across Metz / BDB panels / PDSP | minutes | none (design only) |
| `research/meta_fewshot/r2_reserved_fiber_section.py` | **M1** RFMS module, control battery, Gate specification, collapse check | — | none until E4 |

Run order: `E1 → E2 → E3 → E0 → E4 → E5`.

---

## G. Preliminary experiment result — **not obtained**

**Execution was blocked for the whole session.** The Linux workspace failed to
start on every attempt with:

```
failed to set session disk path: session disk not found:
...\vm_bundles\claudevm.bundle\sessiondata.vhdx
```

No Python could be run, so E0–E3 were written but not executed, and no new
numbers appear in this report. Every quantitative statement above is derived
either from the frozen source code, from the published result tables in
`MAIN_V0_REPORT.md` / `V1_DEVELOPMENT_REPORT.md` / `BIOLOGICAL_GAUGE_AUDIT_REPORT.md`,
or from the dataset manifests — never from an unrun computation.

To execute, in the `drug` environment from the repository root:

```bash
python -m research.meta_fewshot.r2_e2_tbasis_decomposition
python -m research.meta_fewshot.r2_e3_crossed_census --list-keys   # confirm schema
python -m research.meta_fewshot.r2_e3_crossed_census
python -m research.meta_fewshot.r2_e0_regime_audit --split meta_val
```

The first three open no affinity label. E0 is descriptive on already-consumed
splits. None of them can open a Gate.

---

## H. Final recommendation

```text
IDENTIFIABILITY_REGIME_DEFECT_IDENTIFIED_REPAIR_DERIVED_EXECUTION_BLOCKED
```

Nearest listed verdict: **`BIOLOGICAL_INFORMATION_EXISTS_BUT_REPRESENTATION_REPAIR_REQUIRED`**,
qualified in two ways that the listed options do not express:

1. **The primary defect is not representational.** At `d = 2 < k = 5` the v0
   estimator is `GL(d)`-invariant and its own frozen theory (F20/CP-3, CI-A3)
   makes the auxiliary channel redundant. Representation repair alone cannot fix
   this, which is why V1 failed. The repair is the **reserved-fiber
   identifiability budget** (`d_S ≤ k < d`, `d_c ≥ 1`), which makes partner
   specificity mathematically necessary for the meta-training objective rather
   than an add-on — the outcome §16 asks for.
2. **`EXISTING_DATA_INSUFFICIENT_FOR_PARTNER_IDENTIFICATION` is not supported.**
   A1's negative was obtained on a substrate with 21 dependency components, 86.4 %
   in one. Metz (32,849 cells, density 0.241, ~31,775 interaction df) and
   pocket-85 ESM embeddings are already on disk. The A0/A1 "acquire data" clause
   should be amended to "use the panels already built".

Immediate next actions, in order:

1. Restore shell access and run **E1 and E2** — label-free, minutes, and they
   decide the diagnosis outright.
2. Run **E0**. If the intercept null matches or beats `full_correct`, withdraw
   the biological reading of `MAIN_V0_RESULT.json` and add the intercept arm to
   every future battery as a standing control.
3. Run **E3** to confirm the crossed supply, then amend
   `BIOLOGICAL_GAUGE_AUDIT_REPORT.md` §"A2 and next action".
4. Only then preregister **E4** (RFMS development on source + `meta_val`) and,
   separately, **E5** on the frozen BDB-panels/PDSP confirmation supply.

Do not connect CSMO, do not unfreeze the frontend, and do not migrate anything to
production on the strength of this document.

---

### Sources

- [Meta-learning Adaptive Deep Kernel Gaussian Processes for Molecular Property Prediction (ADKF-IFT)](https://arxiv.org/abs/2205.02708)
- [Multi-task bioassay pre-training for protein–ligand binding affinity prediction (MBP)](https://academic.oup.com/bib/article/25/1/bbad451/7469349)
- [Physicochemical graph neural network for learning protein–ligand interaction fingerprints from sequence data (PSICHIC)](https://www.nature.com/articles/s42256-024-00847-1)
- [MetaDTA: Meta-learning-based drug-target binding affinity prediction](https://openreview.net/forum?id=yzlif16IASM)
- [A meta learning and task adaptive approach for drug target affinity prediction (AdaMBind)](https://www.nature.com/articles/s41467-026-70554-5)
- [Fast Context Adaptation via Meta-Learning (CAVIA)](https://proceedings.mlr.press/v97/zintgraf19a/zintgraf19a.pdf)
- [A Hybrid Structure-Based Machine Learning Approach for Predicting Kinase Inhibition by Small Molecules](https://pubs.acs.org/doi/10.1021/acs.jcim.3c00347)
- [Profiling Prediction of Kinase Inhibitors: Toward the Virtual Assay](https://pubs.acs.org/doi/abs/10.1021/acs.jmedchem.6b01611)

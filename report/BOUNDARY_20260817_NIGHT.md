# Boundary, 2026-08-17 (night): the k=0 level wall and the scope of the conclusion

This replaces the boundary statements in report/BOUNDARY_20260816.md and
task.md for the questions below. All numbers are development evidence on the
governed BindingDB-Ki double-cold protocol; meta_test remains sealed
(logical exclusion after parsing, 768 cells withheld, never evaluated).

## 1. What the new measurements establish

### The level term is assay history, and the protocol removes the transferable part

- Within meta_train, a target's level transfers between different proteins
  **inside one document** with R^2 = +0.451 (leave-one-target-out, MSE 1.010
  vs grand-mean 1.912; D0b_DOC_TRANSFER.json).
- The double-cold split has document_overlap = 0, so exactly this signal is
  unavailable at inference. What remains transferable: panel composition
  (held-out R^2 +0.239), protein sequence embedding (+0.119), jointly
  +0.259 (D0_LEVEL_ANATOMY.json). At most ~26% of level variance is
  predictable from legal inputs — before any model or training method.
- Best legitimate level predictors on frozen meta_val banks: ESM-650M linear
  probe 1.6875, panel-feature MLP 1.8868, trained panel head 1.438
  (panel-shuffled control 1.539); meta_train-only constant 2.15-2.17;
  meta_val-calibrated REFERENCE 1.3471. The k=0 budget is level MSE <= 0.1239
  at the measured shape term — more than 10x beyond the best measured
  predictor.

### The decomposition re-audit

- The Stage C level/shape split is per episode (per drawn panel), not per
  canonical target; the drawn-panel contribution to level^2 is small
  (0.013-0.034 pK^2), so the level term is genuinely between-target.
- The calibrated constant (1.3471) reads meta_val labels (disclosed
  REFERENCE); the honest meta_train-only constant is 2.15-2.17, which the
  tested features do beat.

### The shape term is representation-limited

- The occupancy ordering signal (r +0.203 meta_val) survives scaffold
  novelty (+0.154) and low ligand recall (+0.221), but its full exploitation
  is worth ~0.04 pK^2.
- No arm in any stage has ever produced a resolved shape/centered or ranking
  improvement in a multi-seed comparison. The orthogonal level/shape routing
  (Stage E) moved centered MSE by at most 0.036 and never below the
  baseline's 0.865 at k=0.

## 2. The falsification ledger (multi-family, all leak-free protocols)

Training/framework families falsified: analytic and legacy operators
(QPSMP/LIRMS, HyperSAR, D-MEMT/DORM, CIPF/TERM, K3/ELMT); BPSF/CIPF pair
trunk; contact-grammar trunk; moment-form support adaptation (A2 exact
operator); inner/outer-loop meta-learning (Stage A/B); centered-objective
protein conditioning (Stage P); panel-set level head + orthogonal
level/shape routing (Stage D/E, with loss-only and framework-only
ablations); pairwise learned transport (Stage F, with framework-only
ablation); ranking-loss substitutions (R9-R14); ESM-150M and ESM-650M
residue-input trunks (Stage G/G2, 3 seeds each arm).

External representations measured: ESM-150M pooled, ESM-650M pooled,
panel composition (handcrafted and learned), assay covariates (endpoint,
panel/replicate counts, document identity), and structure-derived pocket
priors from homologous holo complexes (209/387 targets covered at >=30%
identity + >=50% query coverage; pocket-probe level MSE 2.4398 vs 2.6179
constant, shuffled-pocket control 2.4941 — rejected at the identifiability
gate, H0_POCKET_IDENTIFIABILITY.json). End-to-end LM conditioning (Stage I:
live ESM-150M with LoRA adapters, single-stage training) produced two
resolved RANKING gains (k=2 Spearman vs its frozen live control, k=3 Pearson
vs T2) but no resolved MSE movement and slightly worse level — rejected by
its own gates. Assay provenance (Stage J: journal/publisher codes parsed
from the legal panel_ids metadata, D0c probe level MSE 1.619 vs 2.155
constant, shuffle 2.522) is the strongest single covariate family found, and
the trained assay-aware level head reaches level^2 ~1.30 at k=0 (the best on
record) — but coupling it to the k>=1 transport degrades ordering with
RESOLVED intervals (k=2 Spearman -0.0624, k=3 -0.0598) and it is rejected by
its own gates. The legal covariate space for the k=0 level is now
essentially exhausted; MSA remains blocked on a governed UniRef snapshot.
The contrastive coembedding family (Stage K/K2) was then tested: the K-REG
configuration produced the first ALL-k resolved MSE improvement across three
seeds (k=0 -0.1118 [-0.1851, -0.0490], k=1 -0.0480, k=2 -0.0273, k=3
-0.0218, k=5 -0.0122; ranking preserved; zero control inversions), but its
k=0 centered gain did not survive pooling, so it is NOT confirmed and
nothing is promoted. It is read as calibration-consistency regularization,
not a new information source: per-seed k=0 MSE stays >= 2.44.
The last composition was then tested (Stage L): a support-gated version of
the assay-aware level head (active at k=0 only). Its gate kept k>=1 MSE
statistically indistinguishable from T2 and produced the best k=0
calibration in the record (MSE 2.0997, level^2 1.2151), but ordering
degraded with RESOLVED intervals at k=2/3/5 and k=0 CI — training the head
reshapes the shared trunk, so zero-shot level and within-target ordering
conflict on the same representation no matter how they are routed. Three
compositions (E, J, L) have now failed this gate; a separate frozen-feature
calibrator has a measured ceiling below L. A fourth composition (Stage Q:
a level head over FROZEN features — ESM bank vector, handcrafted panel
statistics, journal table — gated to k=0) was then falsified: the best
frozen joint probe on record (1.3416) does not survive training, and the
trunk's ranking still degrades with resolved intervals at k=0/2/3. The
conflict between zero-shot level and within-target ordering on one shared
trunk is therefore fundamental to single-stage end-to-end training: the only
escape is a separately trained inference calibrator, which the governing
contract excludes as a multi-stage regime. Independent external validation:
Nelen et al. (J Cheminform 17:8, 2025, PMID 39833966) measured on ChEMBL
that absolute values from different assays are rarely comparable while
potency differences between matched pairs are robust — the same
level/ordering structure this boundary established on BindingDB. The
protocol-level conclusion stands as the programme's final state. The
ligand-side language-model family
was then measured and falsified (Stage M0: frozen ChemBERTa-77M ligand
embeddings; within-target ordering r +0.147 [-0.026, +0.318] below the
occupancy record, level probe collapses to the grand mean exactly). The
external-representation ledger now covers every locally testable legal
input family. The protein-function-annotation family was then measured and
falsified (Stage P0: ProteinKG25 GO bags, 313/387 targets matched; level
probe 2.27 vs 1.43 constant on the covered subset), so the bounded
conclusion is final for this machine and protocol.

## 3. The conclusion, and exactly where it applies

**Under the governed double-cold protocol with sequence + 2D ligand inputs,
ordinary end-to-end training, and the legal-input families above, k=0 MSE
<= 1.00 pK^2 is not achievable**: the best trained model's level^2 alone
(1.52-1.56) exceeds the entire 1.00 budget, and the strongest legitimate
level predictor signal (within-document assay history) is structurally
unavailable at inference by the split's document closure. This is a
protocol-and-input conditional conclusion, not a claim about all conceivable
models or datasets.

What remains reachable in development: k=5 at 0.944-1.007 across three seeds
with clean label/protein controls; k=3 at 1.178-1.251; k=2 at 1.272-1.363.

What would change the verdict: (a) a governed MSA lane; (b) a structure or
pocket lane covering the corpus; (c) end-to-end ESM fine-tuning; (d) a
contrastive coembedding framework; (e) restating the k=0 target as
centered-MSE/ranking or admitting per-document calibration — all outside the
current protocol and none authorized by this document.

## 4. Hygiene

Every stage in this cycle preregistered before training; single-seed screens
gated by stop rules; the one promising lane (G) went to a preregistered
three-seed confirmation and was not confirmed; no sealed meta_test label
entered any fitting, selection or reported metric; query labels stayed
loss-and-metric-only; no closed-form solvers, no
cross-dataset support, no transduction. Davis/KIBA independent training was
not authorized because no candidate passed the promotion gates.

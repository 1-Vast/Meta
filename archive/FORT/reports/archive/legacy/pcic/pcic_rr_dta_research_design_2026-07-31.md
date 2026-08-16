# PCIC-RR-DTA: research design, falsification protocol, and current decision

Date: 2026-07-31

## Executive decision

`NEEDS_NEW_PUBLIC_INFORMATION_SOURCE`

Do not implement or train CITA-DTA, SSNM-DTA, CRISP-Mamba-DTA, or the design
below on the current FORT labels. The software and CUDA stack are ready, but
the scientific admission test is not. Existing TRAIN-only evidence has already
triggered the cheapest STOP condition:

- ChEMBL-37 supplies many algebraic rectangles but no verified same-assay
  rectangles, weak difference-in-differences signal relative to q90 replicate
  noise, and strong provenance-family concentration.
- Papyrus 05.7++ is aggregated to one parent-target row and cannot recover
  independent document-resolved replication.
- The minimal correct-protein probe is worse than its protein-free control.
- Earlier support-memory and trainable interaction probes did not distinguish
  correct from wrong protein.

The conditional main design is **PCIC-RR-DTA**: Provenance-Closed Interaction
Contrast with Residual Reordering. It has exactly two core innovations:

1. Provenance-Closed Cross-Target Interaction Contrast (PCIC).
2. Counterfactual Evidence-Gated Contrast Posterior (CEGCP).

Frozen ESM-2 features, Morgan/physicochemical ligand features, the small
interaction operator, and the ligand baseline are infrastructure, not claimed
innovations. A standard first-order MAML adapter on the same admitted PCIC
features is the only backup. SSNM multi-scale Mamba and CRISP recurrence remain
capacity-matched controls after, and only after, a new information source
passes D0 and P0.

## Assumptions and verifiable success criteria

Assumptions:

- The scientific estimand is target-conditioned ligand reordering, not merely
  lower global affinity error.
- The primary endpoint is pKi and the primary support size is k=5. Every
  admitted model is also evaluated at `k in {0,1,3,5,10}` with frozen support
  selection. pKd is a separately powered secondary endpoint. pIC50 is never
  pooled with either.
- A public record is admissible only when its values, version, license, schema,
  provenance, and checksum are recoverable. A readable paper is not enough.
- Confirmation, Davis target-conditioned labels, and other sealed labels stay
  unread. Structures, poses, pockets, private teachers, and future measurements
  remain closed.

Success means that a new public source first passes the model-free D0 topology
gate, a frozen low-capacity P0 probe then establishes an attributable
correct-protein effect, and a one-seed M0 model exceeds both statistical MDE80
and material-effect floors on error and ranking. No architecture is admitted
by aggregate RMSE alone.

## 1. Root-cause diagnosis

### 1.1 What the current model does and does not do

The current architecture is not a pure late-concatenation model. It first
encodes protein residues, then a ligand projection controls residue-attention
weights, and the pair feature contains pooled protein, ligand, product, and
absolute-difference terms. The accurate diagnosis is narrower:

> The current interaction is a single, one-way ligand-to-residue pooling step;
> the residue state is not updated by the query ligand, and ordinary affinity
> supervision has not made the resulting representation protein-specific.

Changing that operator may be useful only after the labels identify a
protein-specific effect. The present failure cannot be attributed to the
absence of any ligand-protein interaction code.

### 1.2 Decisive observed behavior

The TRAIN-only IDG-RBP probe used 58 strict pKi gate episodes, 446 fit targets,
413 fit homology components, and only 256 bilinear parameters.

| Arm | RMSE | Spearman | Pairwise accuracy |
| --- | ---: | ---: | ---: |
| B0 ligand-only | 1.3598 | 0.0865 | 0.5334 |
| Protein-free calibration | 1.3550 | 0.0865 | 0.5334 |
| Correct protein | 1.4266 | 0.0520 | 0.5182 |
| Random wrong protein | 1.3965 | 0.0946 | 0.5359 |
| Matched wrong protein | 2.1465 | -0.4382 | 0.3075 |

Correct-minus-protein-free RMSE gain is `-0.0709`, with component-bootstrap
95% CI `[-0.1221, -0.0200]`. The matched-wrong arm covers only five targets in
two components and cannot support attribution. The machine-readable decision
is `promote=false`. The earlier AnchorDelta correct and wrong-protein outputs
were nearly identical, and the existing memory gate is `STOP`.

### 1.3 Causal chain of the failure

1. Global ligand potency is highly predictive, so ordinary affinity loss can
   be minimized by B0(d) and a target intercept or calibration slope.
2. With k=5, removing an intercept and the B0 calibration direction leaves at
   most three identifiable support contrasts. A high-dimensional prompt,
   fast-weight memory, or full MAML update is therefore underidentified.
3. The current public measurement graph does not provide verified
   protocol-comparable crossed target-ligand observations. Raw rectangle count
   greatly overstates independent information.
4. A scheduler can preferentially sample informative tasks, but it cannot
   create protein-specific information when no admitted task contains it.
5. More Mamba scales, recurrence depth, graph tokens, cross-attention width, or
   epochs change capacity and inductive bias. They do not repair steps 2-3.

Therefore the missing prerequisite is **new public, provenance-comparable
interaction information plus a task construction that isolates reordering**.
PCIC addresses the task construction; CEGCP limits adaptation to directions
that the support set can identify. They are conditional on the new information.

### 1.4 Decision on the three supplied proposals

| Proposal | Engineering feasibility | Scientific feasibility now | Reason |
| --- | --- | --- | --- |
| CITA-DTA | Feasible | STOP at D0/P0 | The contrast idea targets the right estimand, but current crossed labels lack protocol closure and the frozen correct-protein probe fails. |
| SSNM-DTA | Feasible with nontrivial Mamba scan changes | STOP | Multi-Delta protein dynamics and fast memory add capacity without admitted protein-specific supervision; the proposed rank bound is not automatic. |
| CRISP-Mamba-DTA | Feasible after adding ligand tokens | STOP | Recurrent cross-modal refinement is an unvalidated transfer and would simultaneously change ligand representation, fusion, recurrence, and head. |

## 2. Public-data audit, roles, overlap, and leakage

### 2.1 Frozen source ledger

| Source | Version and checksum | License | Schema and provenance | Role | Overlap/leakage risk and status |
| --- | --- | --- | --- | --- | --- |
| ChEMBL | Release 37; dual-cold registry SHA-256 `0e754f73f5d75913d61791d6ccd08e05662cd8015fc608ba370d4ee303e6b784` | CC BY-SA 3.0 | `chembl-dualcold-v1`; exact pKi/pKd rows with target, ligand parent, scaffold, document, assay, source row, split, and homology component | TRAIN/development estimator engineering only | 343,211 registered rows; 274,400 audited TRAIN rows. Same database supplies historical panels, so document, assay, target, connectivity, scaffold, and component closure are mandatory. pIC50 excluded. |
| Metz/ChEMBL document | ChEMBL-37 doc `CHEMBL1201862`; registry SHA-256 `94da6bb5a59c2911672fde982530c8dd6a673c194b2b2d7b4638df7768c8173e` | CC BY-SA 3.0 | `chembl-kinase-panel-dualcold-v1`; single-document pKi kinase panel with target, ligand, scaffold, assay, homology, and frozen split | Spent development only | One campaign and kinase family; 64 accessions and 618 ligands overlap ChEMBL TRAIN. Development outcomes were already used, so it is neither new pretraining information nor confirmation. |
| Papyrus++ | 05.7, README dated 2024-10-24; activity SHA-256 `8004e0d1027a760f205b45264386f792e7d49658da39f77f52e660a6f19760dd`; target SHA-256 `832e564fb82daea0e4da79abcb44834d10104229382874e79915a1288d80783c` | **Missing locally**: README points to an absent `LICENSE.txt` archive | Aggregated parent-target table with source/doc/assay lists, endpoint flags, aggregate pChEMBL values and protein metadata | Provenance/topology audit only | Includes ChEMBL 34 and other sources; not independent of ChEMBL. Semicolon provenance cannot be split into pseudo-records. Zero document-replicated parent-target cells after strict filtering. STOP. |
| BindingDB native articles | 202607 ZIP; SHA-256 `d2584d1519318d00ab5f46289da5ab3549affe732d598a5072f8777b6b3b5262` | CC BY 3.0 US | Article-resolved exact human single-chain Ki/Kd subset; target and ligand identifiers present; no source-native assay identifier | Source recovery only | 93,712 rows, of which 93,023 are BindingDB-curated and 429 imported from ChEMBL. Only 38 pKi and 6 pKd targets have an article block with at least 40 ligands. Target-shallow and not powered; assay comparability unresolved. |
| Reinecke Kinobeads | 2024, DOI `10.1038/s41589-023-01459-3`; all five raw workbook hashes are frozen in `raw/sha256.txt`; registry SHA-256 `92ca4981efb87ee42ae03fd4f6d2b2a0fc11ee90a678921a1c67d6b0c2768405` | CC BY 4.0 | `reinecke-2024-dualcold-development-v1`; `pKd_app = 9-log10(Kd_app,nM)`, target/accession, ligand parent/scaffold, homology, document, reliability, and split | Development only; already observed | One kinase campaign, not independent confirmation. pKd_app is not pooled with biochemical pKd. Development labels were consumed historically; confirmation labels remain unread. |
| Novartis SPD | 2023, Zenodo `10.5281/zenodo.8103950`; primary activity SHA-256 `7132723f85e746de2f8387d01dcde6ffff703c92561fda9751cbd6753e900240`; all file hashes in manifest | Data CC BY 4.0; code MIT | `fort-source-manifest-v1`; summarized IC50/AC50, censoring, assay/group, gene map, tested negatives, and provenance | Source-shape/censoring audit only | 91.3% censored inactive; median 14 compounds per gene-assay; marketed drugs overlap ChEMBL heavily. Not exact pKi/pKd and not an affinity-training source. |
| Davis/ChEMBL document | Davis 2011 pKd, ChEMBL doc `CHEMBL1908390`; registry SHA-256 `f15daa5478f63a648a07d52d76aee588e4dc6d7275444fc50d204a774a3499fe` | ChEMBL redistribution under CC BY-SA 3.0; original-paper rights are not used to broaden this role | pKd-only, self-contained anchor/query split; 116 targets, 102 homology components, scaffold/connectivity-disjoint query | Registered single-use confirmation only | Underpowered: median 12 query ligands and target-conditioned MDE80 0.1596. `consumed=false`; target-conditioned labels remain sealed and unscored. |
| KIBA benchmark | Exact local release/version and checksum **not registered** | **Not verified for a local training artifact** | Composite KIBA score combines heterogeneous Ki, Kd, and IC50-derived evidence; row-level protocol provenance is insufficient for this estimand | Historical comparability only; no label use | Overlaps common kinase/ligand benchmarks and violates endpoint separation for the primary task. It is not an admitted source. |
| Binder2030 | Paper reports 3,384 ligands by 400 membrane proteins; exact downloadable release/version and checksum **unresolved** | **Unresolved** | Full value matrix, endpoint details, assay/protocol fields, provenance blocks, and row identifiers **unverified** | Candidate reopening source only | Paper-level claims do not establish public download, rights, completeness, independence, or training admissibility. No labels may be used until all missing fields close. |

KIBA and Metz provide no new independent information under their registered
roles.

### 2.2 Current ChEMBL topology result

The registry-closed audit is TRAIN-only and keeps pKi and pKd separate.

| Quantity | pKi | pKd | Combined where meaningful |
| --- | ---: | ---: | ---: |
| Raw endpoint rows | 248,775 | 25,625 | 274,400 |
| Target-pair/document units | 12,661 | 22,640 | 35,301 |
| Algebraic rectangles | 15,142,136 | 983,933 | 16,126,069 |
| Same-assay rectangles | 0 | 0 | 0 |
| Unit median absolute DD | 0.2000 | 0.5916 | 0.3932 |
| DD/noise(q90) | 0.802 | 0.296 | 0.633 |
| Unit reversal fraction | 0.360 | 0.413 | 0.394 |
| Largest document-family fraction | 63.3% | 72.5% | 69.2% |

Zero same-assay rectangles can partly reflect target-specific ChEMBL assay IDs;
it is not a biological null. It nevertheless means that protocol equivalence
has not been demonstrated. The 16.1 million rectangles are dependent algebraic
combinations, not 16.1 million experiments. The independent count is at most
the 35,301 target-pair/document units and is smaller after provenance and shared
component clustering.

Papyrus 05.7++ has 707,461 aggregate parent-target rows and 147,434 strict
retained rows, but exactly zero document-resolved replicated parent-target
cells. It cannot repair the missing provenance by volume.

### 2.3 Required split and source firewall

For every endpoint, construct splits in this order and freeze their hashes:

1. Normalize parent connectivity and exact target accession/construct.
2. Assign endpoint without conversion across Ki, Kd, IC50, or apparent Kd.
3. Union duplicate source rows and mark all cross-database lineage.
4. Build target homology components before target split.
5. Build ligand connectivity, Bemis-Murcko scaffold, and chemical-neighbor
   components before support/query split.
6. Close document, assay/protocol, campaign, and source families across fit,
   gate, and evaluation roles.
7. Freeze episode membership, support draw, query rows, source-row IDs, and
   checksums before any model score.

## 3. Literature and transfer-evidence matrix

Evidence classes are: **direct** (the paper validates the stated method in DTA
or the exact diagnostic), **structural transfer** (the mathematical structure
is established elsewhere), and **unvalidated hypothesis** (the proposed DTA
use is new and must be falsified).

| Module or claim | Primary paper | Class | What the paper supports | Failure boundary |
| --- | --- | --- | --- | --- |
| Joint ligand-target modeling | Lapinsh et al., Bioinformatics 2005, DOI `10.1093/bioinformatics/bti703` | Direct | Proteochemometric modeling across ligands and targets | Joint input does not prove protein contribution under cold splits. |
| Protein permutation gate | Avdiunina et al., JCIM 2025, DOI `10.1021/acs.jcim.5c00395` | Direct | Rigorous PCM evaluation and small protein-embedding contribution under permutation | Does not supply a remedy; motivates wrong/protein-free controls. |
| Public aggregate source | Bequignon et al., J. Cheminformatics 2023, DOI `10.1186/s13321-022-00672-x` | Direct | Papyrus curation and aggregation | Aggregation cannot reconstruct independent raw provenance. |
| Local ligand edits | Hussain and Rea, JCIM 2010, DOI `10.1021/ci900450m`; Hu et al., JCIM 2012, DOI `10.1021/ci3001138` | Direct | Matched molecular pairs and MMP cliffs | Does not validate cross-target affinity DD pretraining. |
| Cliff-specific evaluation | van Tilborg et al., JCIM 2022, DOI `10.1021/acs.jcim.2c01073` | Direct | Standard models fail on activity cliffs | Cliffs may be assay noise or source effects without provenance closure. |
| Four-cell biochemical contrast | Hudspith et al., OBC 2019, DOI `10.1039/C9OB01558B` | Structural transfer | Mixed double-mutant cycles use protein variant by ligand edit double differences | Small controlled biochemical cycles do not validate database-scale sequence transfer. |
| Difference-in-differences | Callaway and Sant'Anna, J. Econometrics 2021, DOI `10.1016/j.jeconom.2020.12.001` | Structural transfer | Contrast algebra and dependence-aware inference | Heterogeneous assays violate causal comparability assumptions. PCIC is not called causal. |
| PCIC objective | This proposal | Unvalidated hypothesis | Residualized cross-target DD may force target-conditioned ligand reordering | Fails if same-protocol units, SNR, or independent components are insufficient. |
| Sequence interaction operator | Koh et al., Nature Machine Intelligence 2024, DOI `10.1038/s42256-024-00847-1` | Direct for sequence DTA; transfer for this implementation | PSICHIC learns physicochemical sequence-ligand interaction fingerprints | Fingerprints are not experimental contacts, and training-lineage closure must be audited. |
| Ligand shortcut risk | Mastropietro et al., Nature Machine Intelligence 2023, DOI `10.1038/s42256-023-00756-9` | Direct | Affinity models can memorize ligands and show limited protein dependence | Supports controls, not a new encoder. |
| Gradient meta-learning | Finn et al., ICML 2017, arXiv `1703.03400` | Structural transfer | MAML episodic adaptation | k=5 cannot identify unrestricted updates. |
| Few-shot DTA and task scheduling | Wan et al., Nature Communications 2026, DOI `10.1038/s41467-026-70554-5` | Direct | AdaMBind uses protein-as-task, MAML, query loss, and support/query gradient similarity | Its reported protocol does not establish FORT's scaffold/document/assay/component closure or wrong-protein attribution. |
| Counterfactual task utility | Lin et al., Nature Communications 2025, DOI `10.1038/s41467-025-66915-1` | Structural transfer | TAPB uses intervention-style debiasing in binary DTI | Continuous DTA wrong-protein evidence gating is unvalidated. |
| Closed-form few-shot posterior | Patacchiola et al., NeurIPS 2020, arXiv `2008.05414` | Structural transfer | Deep-kernel Bayesian adaptation | Does not prove the support residual directions are protein-specific. |
| CEGCP | This proposal | Unvalidated hypothesis | A support-only Bayes-evidence gate may prevent unsupported reordering | Fails when correct and wrong proteins have indistinguishable support evidence. |
| Multi-scale SSM | Karadag et al., Neurocomputing 2026, DOI `10.1016/j.neucom.2026.133226` | Structural transfer | ms-Mamba validates multiple Delta scales in time-series forecasting | Residue index is not physical time; no few-shot DTA evidence. |
| Fast weights | Schlag et al., ICML 2021, PMLR 139 | Structural transfer | Delta-rule associative memory | Does not establish that support labels contain target-specific information. |
| Protein Mamba feasibility | Peng et al., Nature Methods 2025, DOI `10.1038/s41592-025-02656-9`; Sgarbossa et al., Bioinformatics 2025, DOI `10.1093/bioinformatics/btaf348` | Direct for protein sequence modeling | Bidirectional/long-context Mamba is feasible on biological sequences | Neither validates multi-scale ligand-conditioned few-shot affinity. |
| Recurrent cross-modal fusion | Zou et al., Medical Image Analysis 2025, DOI `10.1016/j.media.2025.103549`; Wang et al., AAAI 2025, DOI `10.1609/aaai.v39i8.32879` | Structural transfer | Cross-modal MRI fusion and multimodal prompting | Aligned images/global prompts are not k-labeled protein-ligand tasks. |
| State caution | Merrill et al., ICML 2024, PMLR 235 | Direct theoretical counterevidence | An SSM state does not imply successful state tracking | Repeating a static pair cannot be assumed to create biochemical interaction state. |

The novelty claim is deliberately narrow. Four-cell algebra, Bayesian linear
updates, sequence encoders, and counterfactual controls all have precedents.
The potentially novel contribution is their provenance-closed use to identify
and safely adapt only a target-dependent reordering signal in strict dual-cold
DTA.

## 4. Main design and one backup

### 4.1 Main: PCIC-RR-DTA

**Innovation 1 - PCIC.** Pretrain a small interaction residual using
endpoint-specific, same-protocol crossed target-ligand units after cross-fitted
removal of ligand potency, target intercept, and target-specific B0 slope.
Train on the four-cell residual contrast, not on an unqualified union of
affinity rows.

**Innovation 2 - CEGCP.** At a new target, fit at most three support-contrast
directions with a closed-form posterior. Enable the ranking correction only
when support-only marginal evidence favors the correct protein over calibration
and deterministic wrong proteins. The gate is exactly zero when evidence,
rank, or numerical stability is insufficient.

Non-innovation infrastructure is intentionally small:

- frozen ESM-2 residue embeddings;
- current Morgan-1024 plus 10 physicochemical descriptors;
- the current ligand-conditioned residue pooling, reduced to a low-rank
  interaction feature;
- an explicit frozen ligand-only B0;
- endpoint-specific heads and noise models.

A fixed D-MPNN may replace Morgan only in a preregistered representation
control. It is not part of the core claim.

### 4.2 Backup: PCIC-FOMAML

Use the same admitted PCIC representation, B0, episodes, parameter count, and
compute budget. Replace the analytic posterior only with first-order MAML on a
rank-at-most-three residual head and target calibration. Freeze the encoders.
This tests whether gradient adaptation helps after the information gate passes.
AdaMBind-style scheduling is an equal-budget control around this backup, not a
third innovation and not a substitute for PCIC.

### 4.3 Explicit exclusions

- Do not combine PCIC, SSNM memory, CRISP recurrence, graph ligand tokens, and
  adaptive scheduling in one model.
- Do not attribute a gain to architecture if parameters, optimizer steps, or
  negative/counterfactual evaluations are unmatched.
- Do not reopen from a failed gate by increasing width, depth, epoch count,
  support size, or relaxing closure.

## 5. Mathematics, tensors, gradients, and inference

### 5.1 Observation and closure notation

For endpoint `e` in `{pKi, pKd}`, an admitted observation is

```text
o = (t, d, y, e, a, v, s, h_t, c_d)
```

where `a` is assay/protocol, `v` is document/provenance family, `s` is source
row lineage, `h_t` is target homology component, and `c_d` is chemical
component. Endpoint datasets and heads are disjoint. A crossed unit contains
two targets and two exact ligand parents within a verified comparable protocol:

```text
U = (t_a, t_b, d_i, d_j, e, v, a)
```

No fit/gate/evaluation role shares target homology, ligand connectivity,
scaffold/chemical component, document, assay/protocol, or source-row lineage
where that axis is declared cold.

### 5.2 Ligand baseline and cross-fitted nuisance removal

Generate out-of-fold B0(d) predictions on TRAIN under the same held chemical,
document, and provenance components used by PCIC. For every training fold, fit
both B0 and the nuisance calibration without the held component:

```text
y(t,d) = alpha_t + beta_t B0(d) + epsilon(t,d)
r(t,d) = y(t,d) - alpha_hat_t^(-fold) - beta_hat_t^(-fold) B0(d)
```

Only after representation and model selection are frozen may B0 be refit on
all admissible TRAIN rows for deployment. This removal is necessary because a
target-specific potency slope `beta_t B0(d)` can create a nonzero raw double
difference without learning a specific ligand-protein interaction. All
nuisance estimates are cross-fitted; no row predicts its own residual.

### 5.3 Interaction feature

For a batch of proteins and ligands:

```text
H_t       in R[B, L, d_p]       frozen ESM-2 residue features
x_d       in R[B, 1034]         Morgan plus physicochemical features
P_t       = H_t W_p             in R[B, L, d]
q_d       = W_d x_d             in R[B, d]
A_td      = softmax(P_t q_d / sqrt(d), residue axis)
p_td      = sum_l A_td[l] P_t[l] in R[B, d]
z_theta   = W_z [q_d, p_td, q_d * p_td, abs(q_d-p_td)] in R[B, d_f]
```

For a fixed, label-free set `R_episode` of M length/composition/family-matched
proteins outside every correct or counterfactual homology component in the
episode, define the target-contrast feature:

```text
phi_theta(t,d) = z_theta(t,d) - (1/M) sum_{u in R_episode} z_theta(u,d)
```

Any component of `z` that is purely ligand-only cancels exactly. A target-only
constant can remain but is removed by support contrast projection below. The
reference-set construction and weights are frozen from TRAIN metadata and
never use query labels. Correct, matched-wrong, random-wrong, and shuffled arms
within one episode use the identical `R_episode`; otherwise the intervention
would change both the protein and its reference population.

### 5.4 PCIC target and loss

For an admitted unit U:

```text
DDy = [r(t_a,d_i)-r(t_a,d_j)] - [r(t_b,d_i)-r(t_b,d_j)]
q(t,d) = w^T phi_theta(t,d)
DDs = [q(t_a,d_i)-q(t_a,d_j)] - [q(t_b,d_i)-q(t_b,d_j)]
```

Because `DDs` is constructed from a scalar interaction score, it is exactly
antisymmetric under swapping targets or ligands. With reliability weight
`omega_U` frozen from replicate/protocol metadata:

```text
L_PCIC = sum_U omega_U Huber(DDs_U - DDy_U)
L_rev  = sum_{U: abs(DDy_U)>tau_noise} omega_U
         log(1 + exp(-sign(DDy_U) DDs_U))
```

Only TRAIN episodes may add the counterfactual outer loss:

```text
L_anchor = sum_(t,d) omega_(t,d) Huber(q(t,d) - r(t,d))
L_Q(protein) = mean_Q Huber(yhat(protein)-y)
               + lambda_Qrank L_Qrank(protein)
L_outer = L_Q(correct protein)
L_cf = max(0, m - [L_Q(wrong protein) - L_Q(correct protein)])
L_train = L_PCIC + lambda_rev L_rev
          + lambda_anchor L_anchor
          + lambda_outer L_outer + lambda_cf L_cf
```

`L_anchor` fixes the additive gauge left unidentified by four-cell contrasts
and makes the optional k=0 residual score numerically meaningful. It uses the
same cross-fitted nuisance residual, component-balanced weights, endpoint, and
closure as PCIC. It is an anchoring term, not a third innovation. An
absolute-affinity-only model with the same term and operator is a required
control. `L_outer` is the positive episodic objective that trains `A`, the
posterior prior, and the support-to-query path; its Huber delta, ranking loss,
and weights are frozen before M0. NLL is a secondary calibration metric, not an
alternative selected after results. `L_cf` is only an attribution margin and
cannot replace the positive correct-query loss.

To prevent a ligand-deep target-pair block from dominating via quadratically
many pairs, `omega_U` is normalized to sum to one within every target-pair by
protocol/provenance block. A fixed maximum number of ligand pairs per block is
sampled each epoch. Component counts, not rectangle counts, set block weight.

The wrong-protein arm recomputes reference features, support features,
posterior, gate surrogate, and query features. No correct-protein prompt or
feature may survive that replacement.

### 5.5 Support calibration and contrast posterior

For a new target with support `S={(d_i,y_i)}_(i=1..k)`, define

```text
b_i = B0(d_i)
u_i = y_i - b_i - g0_abs q0_centered(t,d_i)
```

where the centered frozen zero-shot PCIC head is enabled only if the separate
no-support absolute gate in section 7.2 passed. At k>=2, let
`C=[1,b] in R[k,2]`; at k=1, `C=[1]`. Estimate endpoint-specific calibration
with an unpenalized intercept and a fixed B0-slope ridge:

```text
P_c = diag(0, lambda_c)               if k>=2
P_c = [0]                             if k=1
gamma_hat = (C^T W C + P_c)^-1 C^T W u
```

Let the columns of `Q` be an orthonormal basis for the null space of
`(W^(1/2) C)^T`, and let `A in R[d_f,3]` be a TRAIN-learned, frozen low-rank
projection. With `F_i = phi_theta(t,d_i)^T A`:

```text
y_tilde = Q^T W^(1/2) u
F_tilde = Q^T W^(1/2) F
r_k     = min(3, k-rank(C), rank(F_tilde))
V_r     = first r_k right singular vectors of F_tilde
F_r     = F_tilde V_r
A_r     = A V_r
Lambda0_r = V_r^T Lambda0 V_r
Sigma_r = (Lambda0_r + F_r^T F_r / sigma_e^2)^-1
mu_r    = Sigma_r F_r^T y_tilde / sigma_e^2
```

For every `k>=3`, a ranking update additionally requires numerical
`rank(W^(1/2)C)=2`; tied or nearly tied B0 support values make the gate exactly
zero. The failure rate is reported separately at k=3, 5, and 10. Only the
`r_k` largest singular directions under a frozen tolerance and tie
rule are retained. Thus k=5 cannot silently write more than three independent
ranking directions.
`Lambda0`, `sigma_e`, `lambda_c`, and singular-value tolerances are learned or
selected on TRAIN components, then frozen.

### 5.6 Support-only evidence gate

For correct and deterministic wrong proteins, compute the contrast marginal
likelihood without query labels:

```text
ell(F_r) = log N(y_tilde; 0,
                 sigma_e^2 I + F_r Lambda0_r^-1 F_r^T)
F_null = 0                              # calibration-only contrast model
E_model = ell(F_correct) - ell(F_null)
E_protein = ell(F_correct) - max_j ell(F_wrong_j)
```

The deployment gate is

```text
g(S,t) = 1[ r_k>=1
            and rank(W^(1/2)C)=2
            and s_min(F_r)>=s_min_threshold
            and cond(F_r^T F_r / sigma_e^2 + Lambda0_r)<=kappa_max
            and E_model>=tau_model
            and E_protein>=tau_protein ]
         * clip(E_protein / tau_saturation, 0, 1)
```

Require `tau_saturation>=tau_protein>0`.
All thresholds are frozen using TRAIN-only nested components. During training,
a smooth sigmoid surrogate permits gradients; deployment uses the hard gate,
which is exactly zero below threshold. TRAIN query labels may optimize an outer
loss and calibrate thresholds, but they are never inference inputs.

### 5.7 Query prediction

For query ligand q:

```text
c_q = [1, B0(q)]                     if k>=2
c_q = [1]                            if k=1
yhat(t,q) = B0(q)
             + g0_abs q0_centered(t,q)
             + c_q^T gamma_hat
             + g(S,t) phi_theta(t,q)^T A_r mu_r
```

The zero-shot head is centered over a frozen TRAIN ligand reference set for
each target, fixing the target-additive gauge:

```text
q0_centered(t,q) = q0(t,q) - mean_(d in D_ref) q0(t,d)
```

`D_ref` is a fixed, query-independent TRAIN ligand reference set whose
connectivity, scaffold, and chemical components obey the evaluation firewall;
its labels are not used in centering.
`g0_abs` is a globally frozen boolean from a separate no-support P0 comparison
against B0, not a per-query choice. It is independent of the support-ranking
gate `g(S,t)`.

| k | Calibration | Adaptation rank | Allowed prediction path |
| ---: | --- | ---: | --- |
| 0 | None | 0 | B0 plus centered zero-shot PCIC only if `g0_abs` passed; under current evidence B0 only. |
| 1 | Unpenalized intercept only | 0 | B0, optional admitted centered zero-shot PCIC, and one exact support offset; no target-specific reordering update. |
| 3 | Intercept plus B0 slope | at most 1 | One gated contrast direction. |
| 5 | Intercept plus B0 slope | at most 3 | Primary three-direction gated posterior. |
| 10 | Intercept plus B0 slope | capped at 3 | Same capacity as k=5; extra supports improve evidence and variance, not rank. |

### 5.8 Gradient paths

The training graph is:

```text
TRAIN labels -> residual/DD targets -----------------------> losses
protein sequence -> frozen ESM -> small interaction theta -+
ligand features -> frozen B0 and small ligand projection ---+
support labels -> calibration/projection -> Cholesky solve -+
TRAIN query labels ----------------------------------------> outer loss only
```

Gradients update only `W_d`, `W_p`, `W_z`, `w`, `A`, and optionally the three
prior precisions. They flow through correct and wrong feature construction and
through the Cholesky solve for `(Sigma_r,mu_r)`. They do not update ESM-2, B0, split
membership, provenance weights, support selection, endpoint transforms, or
hard deployment thresholds. At inference there is no gradient step and no
query label; only a small linear solve is performed.

## 6. Formal shortcut analysis

Let a generic crossed score be

```text
DDs = [s(t_a,d_i)-s(t_a,d_j)] - [s(t_b,d_i)-s(t_b,d_j)].
```

**Ligand-only:** if `s(t,d)=L(d)`, then both within-target differences equal
`L(d_i)-L(d_j)`, so `DDs=0`.

**Target bias:** if `s(t,d)=T(t)`, each within-target ligand difference is zero,
so `DDs=0`.

**Additive ligand plus target:** if `s(t,d)=L(d)+T(t)`, both terms cancel and
`DDs=0`.

**Support offset:** any `c(S,t)` added equally to all query ligands cancels in
within-target differences. It may lower RMSE but cannot change ranking.

**Global scaffold/potency shortcut:** a target-independent function
`G(scaffold(d),B0(d))` also cancels in crossed DD. At adaptation time, projection
off `[1,B0]` prevents an intercept or B0 slope from being counted as ranking
evidence.

**Important non-impossibility:** `s(t,d)=beta_t B0(d)` can yield nonzero DD.
Therefore raw DD alone does not prove interaction. PCIC cross-fits and removes
target-specific B0 slopes, includes a calibration-only arm, and requires
source/scaffold-bin invariance. Likewise, assay or source labels correlated
with target and scaffold can mimic interaction; only protocol/provenance
closure and source-exclusive evaluation address that failure.

A richer shortcut `s(t,d)=a(sequence_family(t)) f(B0(d),scaffold(d))` can also
produce nonzero DD and may beat an unmatched wrong protein. Linear nuisance
removal cannot make this class impossible. The formal claim is therefore
limited to **protein-sequence-conditioned reordering**, not biochemical contact
or causal molecular interaction. Attribution additionally requires strictly
family-matched wrong proteins and a parameter-matched nonlinear
sequence/family-covariate by ligand null. Only an independent interaction label
could support the stronger biochemical-interaction interpretation.

The task is not solved merely because the algebraic DD is nonzero. It is solved
only if the residual contrast generalizes to held homology and chemical
components, correct protein beats wrong/protein-free controls, and correct
support changes query ranking beyond calibration.

## 7. Cheapest audits before implementation

### 7.1 D0: model-free topology and noise audit

For each endpoint independently:

1. Enumerate exact-ligand 2x2 target-ligand rectangles.
2. Require verified same-protocol/campaign comparability; an assay ID alone may
   be target-specific, so map protocols from source records rather than assume
   equality or inequality.
3. Collapse rectangles to target-pair by provenance units; cluster dependence
   through shared homology components, ligand components, and provenance
   families.
4. Report reversal fraction, residual absolute DD, replicate-noise
   distribution, DD/noise at q50/q90/q95, source concentration, and effective
   component count.
5. Run a power calculation for the registered P0 contrasts after all target,
   scaffold, chemical-component, document, assay, and source closures.

D0 PASS requires all of the following for one endpoint:

- downloadable public values, license, version, schema, raw provenance, and
  checksums are complete;
- at least two independent provenance families, with no family over 50% of
  units;
- at least 100 independent held homology/provenance components and sufficient
  post-firewall query depth to make each primary P0 MDE no larger than its
  material floor;
- verified protocol-comparable crossed units exist;
- residual `median_abs_DD / propagated_noise_q90 > 1`;
- pKi, pKd, pIC50, and pKd_app remain separate.

The current ChEMBL result fails protocol closure, q90 SNR, and source
concentration. Papyrus fails raw replication. BindingDB fails target breadth and
assay schema. Therefore current D0 is STOP.

### 7.2 P0: frozen low-rank protein-specificity probe

Only after D0 PASS, run two separately gated probes:

- Freeze B0, ESM features, ligand features, episodes, and all normalizers.
- At k=0, test the centered absolute q0 head against B0 and its wrong,
  shuffled, and protein-free versions. Set `g0_abs=1` only if q0 improves B0 on
  error and both ranking metrics beyond their MDE/material floors.
- At k=5, fit one rank-16 or smaller antisymmetric bilinear support-contrast
  probe on TRAIN components. This separately admits the ranking representation
  and does not set `g0_abs`.
- Compare B0, calibration/protein-free, correct protein, length/family-matched
  wrong protein, random wrong protein, and sequence-shuffled protein.
- Recompute every support/query feature under the replacement protein.
- Use a single seed and paired component bootstrap. Promote only if correct
  protein improves RMSE, Spearman, and pairwise accuracy over protein-free and
  matched wrong protein, with positive 95% lower bounds and effects above the
  registered MDE/material floors.
- Test exact constant-shift equivariance at k=1, 3, 5, and 10. Support-label
  permutation is required only where a ranking update is identifiable; it is
  inapplicable at k=0 and k=1.

The existing P0 consumed 126.8 seconds and 1.86 GiB peak GPU memory on the RTX
4060 Laptop GPU and failed. This is already cheaper and more decisive than a
multi-hour Mamba training run.

The existing TRAIN-only diagnostics can be reproduced without touching their
authoritative outputs by writing to `tmp`:

```powershell
& 'D:\anaconda\envs\drug\python.exe' -m research.interaction_identifiability_audit `
  --registry-closed --range-policy clip --noise-source registry `
  --output tmp/interaction_identifiability_repro.json

& 'D:\anaconda\envs\drug\python.exe' -m research.residual_bilinear_probe `
  --endpoint pKi --seed 1729 --support 5 --ligand-dim 16 `
  --protein-dim 16 --ridge 10 --out tmp/idg_rbp_repro.json
```

Module invocation (`python -m`) is required so repository imports resolve. No
PCIC/SSNM/CRISP training command is authorized while D0/P0 remain STOP.

## 8. Control and ablation matrix

Every arm uses the same endpoint rows, episode hashes, B0, support/query draws,
optimizer-step budget, and evaluation code. Counterfactual arms are charged for
their extra forward passes; compute-matched controls receive the same budget.

| Control | Exact intervention | What it falsifies | Required full-model behavior |
| --- | --- | --- | --- |
| Ligand-only B0 | Remove all target and support paths | Global potency/scaffold shortcut | Full improves error and ranking beyond floor/MDE. |
| Calibration-only | B0 plus support intercept/B0 slope | Support offset and potency rescaling | Full improves within-target ranking; RMSE-only gain is STOP. |
| Protein-free | Remove protein from feature, prompt, gate, posterior, and query | Hidden protein bypass or unnecessary protein | Correct protein must win. |
| Wrong protein | Replace protein everywhere and recompute support posterior/gate | Protein identity attribution | Correct must beat matched and random wrong proteins. |
| Sequence-shuffled protein | Shuffle residues with length/composition preserved | Sequence-order dependence | Correct sequence must win. |
| Wrong support | Use a matched other target's supports and labels; recompute all adaptation | Target-specific support value | Correct support must win. |
| Support-label permutation | Permute labels within the correct support set | Label write versus support chemistry | Correct labels must win; prediction should revert toward gated base. |
| Ligand/support shuffle | Break ligand-label correspondence | Ligand identity leakage | Full gain must disappear. |
| Standard FOMAML | Same rank-3 head, no scheduler, equal compute | Need for analytic evidence gating | Main must be no worse and must have better failure abstention/calibration. |
| AdaMBind-style scheduler | Query-loss/gradient-similarity scheduling, same representation | Whether scheduling alone explains gain | Cannot pass without protein attribution; full must exceed it. |
| Ordinary affinity pretraining | Same operator and rows, absolute affinity loss only | Value of PCIC supervision | PCIC must improve correct-vs-wrong gap and ranking. |
| No PCIC | Random/equally trained interaction feature | Innovation-1 necessity | Expected loss of P0 and reversal gains. |
| No counterfactual loss | Keep posterior/gate, remove TRAIN wrong-protein margin | Counterfactual training contribution | Expected smaller correct-wrong evidence gap. |
| No evidence gate | Always set g=1 | Innovation-2 safety value | More negative-transfer episodes and worse calibration are expected. |
| Gate without wrong protein | Use only model-vs-null evidence | Protein-specific part of gate | Should admit more ligand/calibration shortcuts. |
| Parameter-matched additive model | Replace interaction with widened additive/random feature head | Capacity explanation | Must not reproduce the gain. |
| Compute-matched ensemble | Give baseline equal forward-pass/parameter/step budget | Extra compute explanation | Must not reproduce the gain. |
| Nonlinear family-by-ligand null | Parameter-matched sequence/family covariates interacting with B0, scaffold, and ligand features, without residue interaction features | Family-conditioned ligand shortcut | Full must beat it under strictly family-matched wrong proteins. |
| Source/scaffold probes | Predict source/assay; stratify by scaffold similarity and source-exclusive folds | Provenance/chemical shortcut | Effect must survive low-similarity and source-exclusive strata. |

SSNM's multi-scale Delta and CRISP's recurrent fusion may be added only as
separate `single-scale vs multi-scale` and `single-pass vs recurrent`
parameter-matched controls after PCIC P0 PASS. They may not change ligand
tokens, pooling, recurrence, and parameter count simultaneously.

## 9. Statistics, seed policy, and PASS/STOP rules

### 9.1 Statistical unit and intervals

- For prediction, the biological unit is the held target homology component,
  not a row, ligand pair, support draw, query, or seed.
- For PCIC topology, first collapse to target-pair by protocol/provenance unit,
  then cluster units sharing either target homology component or provenance
  family. Raw rectangles are never resampled as independent observations.
- Use paired resampling: every bootstrap replicate samples components once and
  evaluates all arms on the same sampled components.
- Report 10,000 component-bootstrap replicates, two-sided 95% intervals, the
  number of components, target/query depth, and source-family sensitivity.
- Seeds measure optimization stability, not biological sample size; never
  multiply the component count by the number of seeds.

Primary metrics at pKi, k=5 are target-macro RMSE, within-target Spearman, and
pairwise ranking accuracy. MAE, concordance, reversal accuracy, NLL, interval
coverage, and calibration error are secondary. pKd and other k values are
separate secondary families.

For a paired component contrast `Delta_c`, estimate

```text
MDE80 = (z_0.975 + z_0.80) * sd(Delta_c) / sqrt(C_eff)
      approximately 2.80 * sd(Delta_c) / sqrt(C_eff).
```

Use the cluster-aware effective component count. Before scoring M0, freeze the
larger requirement for every primary metric:

```text
required_effect = max(MDE80, material_floor)
```

Material floors are RMSE improvement 0.10 pK, Spearman improvement 0.05, and
pairwise-accuracy improvement 0.03. Define RMSE improvement as
`RMSE_control-RMSE_model`; ranking improvements use `model-control`.

Core-ablation gates are frozen separately. PCIC versus ordinary-affinity
pretraining must increase the correct-minus-matched-wrong Spearman gap by at
least 0.05 and pairwise gap by at least 0.03, both with positive lower CIs.
CEGCP versus the always-on matched posterior must improve RMSE by at least 0.10
pK and reduce by at least 10 percentage points the fraction of components on
which adaptation is worse than calibration-only; Spearman and pairwise lower
CIs must be nonnegative. Each threshold is also raised to MDE80 when MDE80 is
larger.

### 9.2 Sequential compute policy

1. D0 is CPU/model-free. Stop on data admission failure.
2. P0 is one seed and a low-rank frozen probe. Stop on attribution failure.
3. M0 is one seed of main, backup, and required controls. Stop unless every
   primary mechanism criterion passes.
4. Only after M0 PASS run five total seeds. Require consistent effect sign and
   report seed dispersion without treating seeds as independent biology.

On the verified RTX 4060 8 GiB environment, a future admitted PCIC frozen probe
is expected to fit within the observed 1.86 GiB. A provisional 20,000-unit
PCIC four-arm one-seed experiment is budgeted at 1-4 GPU hours; M0 at 8-24 GPU
hours. These are engineering estimates, not authorization. Memory and wall
time must be measured in a smoke run before the full allocation.

### 9.3 Stage PASS rules

**D0 PASS:** every requirement in section 7.1 passes for an endpoint.

**P0 PASS:** the k=0 absolute and k=5 support-contrast gates are decided
separately. For either gate, correct protein beats protein-free, strictly
family-matched wrong, and shuffled protein on RMSE, Spearman, and pairwise
accuracy; all paired component-bootstrap lower bounds are positive; each effect
exceeds its frozen requirement; source and scaffold sensitivities preserve the
sign. Failure of k=0 sets `g0_abs=0` but does not by itself reject a k>=3
support-contrast route. Failure of the support-contrast gate stops adaptation.

**M0 PASS:** all of the following hold at pKi/k=5:

1. Full PCIC-RR beats ligand-only, calibration-only, standard FOMAML, and the
   AdaMBind-style scheduler on all three primary metrics.
2. Each primary effect exceeds `max(MDE80, material floor)` and has a positive
   paired component-bootstrap lower bound.
3. Correct protein beats matched wrong and protein-free; correct support beats
   wrong support and permuted labels.
4. Within-target Spearman and pairwise accuracy improve together. RMSE-only
   improvement is insufficient.
5. PCIC versus ordinary-affinity pretraining, and CEGCP versus the always-on
   matched posterior, each pass their own preregistered paired component
   contrast: effect above a TRAIN-nested MDE/material threshold and positive
   95% lower CI on their mechanism-primary metric. Arbitrarily small ablation
   differences do not support either innovation claim.
6. A parameter/compute-matched additive model does not reproduce the gain.
7. Effects survive source-exclusive, low-scaffold-similarity, homology, assay,
   and provenance-family analyses.

**Five-seed PASS:** M0 criteria remain satisfied on the prespecified pooled
component analysis, effect signs are consistent across seeds, and no result is
driven by one seed or one provenance family.

### 9.4 Immediate STOP rules

Stop the route if any of these occurs:

- public version/license/schema/checksum/provenance is missing;
- no powered, protocol-comparable crossed units remain after closure;
- correct protein does not beat protein-free and matched wrong protein;
- correct support does not beat wrong/permuted support;
- ligand-only or calibration-only matches the full model;
- only RMSE improves while both ranking metrics do not;
- the effect is below MDE80/material floor or its lower CI is nonpositive;
- performance is confined to near-neighbor scaffolds or one source/assay;
- scheduler utility selects easy tasks but not protein-specific tasks;
- endpoint/source leakage, target intercept, or B0 slope explains the result;
- either claimed core innovation fails its quantitative ablation MDE/CI gate.

After STOP, do not rescue by adding epochs, width, Mamba scales, recurrence,
structure inputs, private labels, or weaker splits.

## 10. Final selection and executable decision

### 10.1 Selection

`NEEDS_NEW_PUBLIC_INFORMATION_SOURCE`

The preferred future route is PCIC-RR-DTA because it first changes the
information-bearing training target and then constrains few-shot adaptation to
the three or fewer directions k=5 can identify. CITA contains the useful seed
of this route. SSNM and CRISP do not currently meet the prerequisite because
they modify capacity before protein-specific information has been admitted.

Current scores, stated separately:

| Dimension | Judgment |
| --- | --- |
| Conditional methodological novelty | Moderate-high: about 6/10; the combination may be new, but each algebraic or architectural ingredient has precedent. |
| Current-data identifiability | Low: about 2/10. |
| Current-data probability of attributable gain | Low: about 2.5/10. |
| Engineering readiness | High; CUDA, Mamba, and tests are healthy. |

### 10.2 What was and was not executed

- Verified environment: `D:\anaconda\envs\drug`, PyTorch `2.6.0+cu124`, CUDA
  12.4, `mamba-ssm 2.2.4`, NVIDIA RTX 4060 Laptop GPU with 8,188 MiB.
- Existing full test result: `76 passed`.
- Existing D0 topology audit and P0 IDG-RBP results were reused; no new sealed
  or confirmation target-conditioned labels were read in this review.
- No model code, dataset, registry, episode, checkpoint, or training entry
  point was modified.
- No new training was launched. Doing so would violate the proposals' own
  cheapest-first falsification contract.

### 10.3 Reopening trigger

Reopen only when one of these is registered:

1. a public, downloadable, source/assay/protocol-comparable crossed affinity
   panel with independent provenance blocks; or
2. a public and authorized non-structural interaction/selectivity label
   satisfying the same target, homology, scaffold, chemical-component,
   document, assay, and source closure.

The current structure, pose, pocket, and privileged-teacher routes remain
closed; this reopening rule does not authorize them.

Binder2030 is the most obvious source to investigate next, but only metadata and
rights recovery is authorized now. Its full downloadable data, license,
checksum, matrix completeness, endpoints, assay/protocol schema, and provenance
must be resolved before reading values or rerunning D0. If it or another source
passes D0, rerun P0 first; do not jump directly to PCIC pretraining.

## Audit trail

This report is grounded in:

- `reports/active/innovation_gate_decision_2026-07-31.md`
- `reports/active/interaction_identifiability_audit_registry_closed_v2.json`
- `reports/active/residual_bilinear_probe_pKi_seed1729.json`
- `reports/active/papyrus_f0.json`
- `reports/active/memorygate.v1.json`
- `reports/active/open_data_only_amendment.md`
- `reports/active/panel_davis_registration.md`
- `reports/active/pretraining_anchor_three_candidate_failure_report_2026-07-27.md`
- `reports/active/cita_ssnm_crisp_multiagent_feasibility_2026-07-31.md`

The three supplied proposal hashes are recorded in the multi-agent feasibility
report. All pre-existing worktree changes remain untouched.

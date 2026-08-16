# FORT Research History: Evidence, Dilemmas, and Decision Record

> **Dense evidence ledger, updated 2026-07-30.**
>
> This document replaces an append-only execution log of about 8,073 lines and 621,557 bytes.
> The pre-rewrite file had SHA-256
> `89e3fbaf762c739fdd211ae46158deae2d4ec4b92239c5f6a46c9a529f59c4c9`.
> Repeated plans, environment setup, transient debugging, superseded interpretations, and duplicate
> conclusions were removed. The ledger retains task definitions, immutable statistics, destructive
> controls, failure reasons, firewall events, formal decisions, reopening conditions, and
> authoritative artifact paths. Committed Git history and frozen files under `reports/active/`
> retain lower-level execution detail.
>
> `task.md` is the authority for current state and the next allowed action. This file explains why
> that state is justified and what evidence any future positive claim must exceed.

---

## 0. Current State and Reading Rules

### 0.0 User-directed two-stage model amendment (2026-07-28)

The active executable program is now `DCST-P0`: PLINDER-to-ChEMBL
destruction-certified spectral transfer under
`reports/active/dcst_two_stage_preregistration_2026-07-28.md`. This amendment
authorizes only the registered filtered PLINDER train/development computation,
ChEMBL train optimization, and ChEMBL development scoring. It does not reopen
confirmation or sealed evaluation and does not alter the historical
substitution-geometry conclusion below. `OMUT-X7` is parked unexecuted.
The current Stage-1 privileged objective and its source-only counterfactual
correction are frozen in
`reports/active/dcst_counterfactual_training_amendment_2026-07-28.md`.
P0's 4,000-step source run certified two of four privileged spectral bands
versus one of four without privileged supervision, but failed the centered
joint mechanism gate. The separately preregistered `DCST-R1` route tests the
identified resolution failure with 32 frozen ESM-2 residue intervals before
any further downstream-label load.
R1 also failed. Its audit established that 64.6% of source rows belong to
multi-accession PLINDER split clusters, only 50.1% contain the arbitrarily
selected cluster-representative accession, and the prior PDB-residue
projection changed bins for 58.3%/74.4% of contacts at 8/32 segments. The
active R2 successor separates the split cluster from the exact pocket-UniProt
model key and maps PDB entity coordinates through RCSB SIFTS.
R2 mapped 13,329 contacts with only 133 unmapped, certified two privileged
bands versus zero without privileged supervision, and made wrong-target
cross-entropy substantially worse. It still failed the centered mechanism
gate because the wrong-ligand margin was only 0.0417. R3 therefore replaces
the single global ligand vector in the structural head with active-Morgan
substructure tokens while preserving the exact entity alignment.
R3 improved structural cross-entropy but certified one band in both privileged
and no-privileged arms. R4 keeps its representation and adds a symmetric
within-target interaction-map retrieval objective to prevent the aggregate
centered loss from hiding ligand-to-ligand confusion.
R4 passed the ligand structural-destruction margin and missed the target
margin by 0.00112, but its privileged affinity certificate count fell to
zero. R5 treats this as multi-task gradient conflict: a structure-only teacher
is frozen before fitting a separate bilinear affinity readout, with matched
random-frozen and no-privileged controls.

### 0.0a Two-stage DCST execution ledger through R18 (2026-07-29)

The paragraph above is the historical entry point, not the current stage.
The complete PLINDER-to-ChEMBL execution now reaches R14. The frozen R6 source
teacher is real source-mechanism evidence, but no tested absolute-state
transport improved strict ChEMBL dual-cold prediction.

| Route | Tested change | Decisive evidence | Decision |
| --- | --- | --- | --- |
| `DCST-P0` | Initial privileged spectral transfer | Privileged certificate `2/4` versus NoPriv `1/4`, but the centered joint structural mechanism gate failed | Stop |
| `DCST-R1` | 32 ESM sequence segments | Found that PLINDER split clusters had been reused as protein identities and resolved PDB indices had been projected onto unrelated representative sequences; the resulting coordinates were invalid | Stop and repair entity identity |
| `DCST-R2` | Exact pocket UniProt plus SIFTS alignment | Mapped 13,329 contacts with 133 unmapped; certified `2/4` versus `0/4`, but wrong-ligand margin was `0.0417 < 0.05` | Stop |
| `DCST-R3` | Active Morgan environments as ligand tokens | Better structural fit, but privileged and NoPriv each certified `1/4` | Stop |
| `DCST-R4` | Bidirectional interaction-map retrieval | Passed ligand destruction and missed target destruction by `0.00112`, while privileged affinity certificates collapsed to `0/4` | Stop; shared-encoder gradient conflict |
| `DCST-R5/R5c` | Structure-only teacher, frozen representation, separate affinity fit | The first run was invalid because a zero-valued autograd anchor exposed teacher paths to Adam weight decay. The corrected run had nonzero weights but `0/4` privileged bands and failed structural attribution | Stop |
| `DCST-R6` | Shared mechanism bottleneck: affinity reads only the predicted `32 x 8` interaction map | Source passed: privileged `2/4` versus NoPriv `0/4`; centered destruction margins `0.06893/0.10704`. Stage 2 scored `0.095585` versus B0 `0.098237`; paired effect `-0.0025 [-0.0110, 0.0038]` | Retain source teacher; stop direct Stage-2 transfer |
| `DCST-R7` | Learned content-addressed protein roles | Router was nearly uniform (normalized entropy `0.9847`, effective roles `7.85`); privileged `0/4` versus NoPriv `2/4` | Stop |
| `DCST-R8` | Frozen ESM role atlas | Privileged `1/4` versus NoPriv `0/4`, but joint atlas-affinity gradients broke the upstream mechanism margins | Stop |
| `DCST-R9` | Frozen R6 teacher plus atlas energy | R6 mechanism reproduced, but privileged and NoPriv each certified `1/4`; NoPriv confidence was much larger | Stop atlas family |
| `DCST-R10` | Continuous position-free content energy | R6 mechanism reproduced; privileged `0/4` versus NoPriv `1/4`, with destroyed-ligand utility exceeding true-pair utility | Stop |
| `DCST-R11` | Source affinity destruction-identification loss | Privileged remained `0/4` versus NoPriv `1/4`; the held-source true pair was not identified against both destructions | Stop source-affinity-energy transfer |
| `DCST-R12` | Frozen structural interaction-moment prompt only | B0 `0.0982`; Priv `0.0966`; NoPriv `0.0989`; Uniform `0.0982`. Priv minus B0 `-0.0017 [-0.0048, 0.0009]` | Stop prompt-only sufficiency |
| `DCST-R13` | Structural prompt fused with a trainable CFRI-compatible direct branch | B0 and Priv both `0.098237`; grouped difference exactly `0.0000 [0.0000, 0.0000]`. Target and ligand destruction left predictions unchanged | Stop raw privileged prompt |
| `DCST-R14` | Label-blind target, ligand, and mechanism transport audit | Pair responsiveness ratio `0.946823` passed, but target support was `5.0089%`, ligand support `9.8500%`, and domain AUCs were `0.908641/0.882452/0.941130` | Stop absolute PLINDER transfer |
| `DCST-R15 / SISMT` | Intersect common source/target covariance directions with privileged R6 certificate directions | Priv had 143 support-compatible directions, but only one met mechanism overlap (`0.050474`) and its 20-bootstrap median subspace projection was exactly `0.0`; stable retained dimension `0` | Stop before partial transport or affinity |
| `DCST-R16 / DTIOD T1` | Transfer local segment-by-active-Morgan mixed finite differences instead of absolute states | Priv exceeded NoPriv by `6480.53x`, random segment by `2.1309x`, and random bit by `12.8637x`; wrong target removed `80.13%`, but wrong ligand increased response by `10.89%` | Stop before student, support gate, or affinity |
| `DCST-R17 / KLIFS bridge` | Replace PLINDER with aligned KLIFS `85 x 7` IFPs and an exact ChEMBL-train overlap bridge | 5,161 firewalled complexes and IFP Jaccard `0.973`, but rectangular core only 393 edges in one giant component; bridge 300 pairs/50 targets/45 homology components; ChEMBL-train ligand support `7.41%` | Stop before affinity or model construction |
| `DCST-R18 / BindingDB PDB source` | Audit an independent native-article, single-chain, PDB-linked source without decoding numeric affinity | Firewall passed, but only 1,267 rows/159 accessions/340 ligands remained; bridge 348 pairs/58 targets/55 homology components; target support `12.70%`, ligand support `7.145%` | Stop; do not read affinity or download PDB coordinates |
| `DCST-R19 / ERIP` | Build a same-domain high-confidence Stage-1 objective from exact `2 x 2` target-ligand affinity contrasts | 19,712 reliable pairs and 527,654 rectangles passed scale/topology/provenance gates, but one target pair contributed `38.7527%` of all rectangles | Stop raw rectangle population; preregister balanced successor |
| `DCST-R20 / balanced ERIP` | Require same endpoint, prohibit one-document four-cell support, and cap target/homology-pair blocks | 324,444 valid rectangles collapsed to 11,871 after caps, spanning 233 target pairs/145 targets; only 256 were pKd and the largest block remained `2.1565%` | Stop exact-rectangle family; no manifest or affinity load |
| `DCST-R21 / HCRR` | Proposed high-confidence residual-rank Stage-1 curriculum | New source-identifiability audit established that affinity reliability alone does not supply exact paired-complex structural deltas or repeated directed edits | Withdrawn before implementation and before any label load; not an empirical failure |
| `PBCNet2.0-D0` | Audit the new open 2026 Zenodo/GitHub release before any package download | Pair member paths and labels are separable, but the 8.6M generated-pair manifest has no original BindingDB row, assay, document, DOI, patent, or source lineage; repository also lacks a license file despite README MIT claim | Stop before 38.4 MB manifest and 5.44 GB archive download |
| `BioLiP2-D1` | Build an affinity-blind exact-complex/contact registry and test source/support topology | 562,794 firewalled complexes, 426,410 same-target ligand pairs, target/ligand support `75.13%/49.82%`; however largest 2-core component was `97.47%` and largest PubMed `2.00535%` | Stop D1; preregister generic-ligand/provenance closure |
| `BioLiP2-D1C` | Remove degree>50 ligands and collapse `(PubMed,target,ligand)` structure duplicates | Closed source retained 66,660 rows/391,262 ligand pairs and target/ligand support `70.30%/49.79%`; PubMed concentration passed, but largest 2-core component remained `79.98%` | Stop as independent RBSDD source; retain as exact-complex/contact index |
| `PSICHIC-G0` | Audit public code/weights, complete XL membership, row-level lineage, and label separation before inference | Apache-2.0, pinned weights, and separated regression/binary/functional labels pass; repository has no XL membership file, linked XL folder exposes only `test.csv`, `finetuning.ipynb`, and `degree.pt`, and documented rows have no source dataset or stable source-record ID | Stop before any dataset row or weight byte |
| `RBSDD source round` | Recover one independent source with exact complexes, real ligand-conditioned deltas, affinity, and provenance | PBCNet2.0 lacks lineage, BioLiP2 retains a `79.9767%` giant 2-core after closure, and PSICHIC lacks public training membership/lineage | Close current public-source round; no RBSDD construction or affinity load |
| `PCIC original proposal` | Enumerate provenance-aware circuits in the nuisance null space as a new interaction-training substrate | `P0-Cycle-A` already projected raw TRAIN exact cells onto the full `ker(Z^T)` for target/ligand/assay/document; pKi projected SD was `0.35596`, pKd top-document energy was `0.93100`, and the verdict was `P0_CYCLE_A_BIOLOGICAL_FAIL_STOP` | Stop the renamed null-space/circuit route; retain only the new label-blind `P_Z_perp X` operator-rank and query-span delta audit |
| `PB-CEC current-source proposal` | Identify shared physical target/ligand factors by coupled completion across PLINDER/KLIFS/BindingDB/ChEMBL-like relations | The completed source round supplies no exact independent row-lineage bridge: PBCNet2.0 lacks source lineage, BioLiP2 fails independent topology, and PSICHIC lacks membership/lineage | Stop current implementation before model or affinity access; retain only as prospective acquisition and bounded sensitivity framework |
| `PCIC-O0-P` | Test whether safe ChEMBL-37 exact cells have enough provenance-separated homology-scaffold-lineage topology to justify operator-rank computation | Outcome-safe projection reproduced 231,090 pKi and 28,749 pKd cells with zero protected values decoded, but pKi/pKd retained only 37/65 joint components and largest components of `99.4690%/95.5129%`; pKi exact-lineage cell coverage was `91.1692%` | `STOP_PCIC_O0_PROVENANCE_OR_TOPOLOGY_INADEQUATE`; skip O0-I and all affinity/model training, then design prospective cycle-closing blocks |
| `RDIB / PD-MVR G0-A` | Test whether BioLiP2 contains cross-PubMed replicated same-construct exact ligand-pair contact-difference units before coordinates or affinity | 133,352 source-specific pairs collapsed to 145 replicated exact differences/59 targets; optimistic conflict-free ceiling 56, largest PubMed share `17.2414%`, and only 5 blocks contained both ligands against the same ChEMBL TRAIN target | `STOP_RDIB_PDMVR_G0A_REPLICATION_TOPOLOGY_INADEQUATE`; stop exact-pair RDIB and PD-MVR structural packing before coordinates/PLIP/model |
| `RDIB-Edit-G0` | Test whether a frozen single-cut directed-edit vocabulary rescues cross-PubMed replication for the same exact target | 730 global edit classes and 200 globally repeated edits collapsed under the correct `(target_key, edit, >=2 PubMeds)` unit to 29 units/17 targets; only 6 repeated an exact ligand pair, and the optimistic resource ceiling was 20 | `STOP_RDIB_EDIT_G0_REPLICATION_TOPOLOGY_INADEQUATE`; strict stereo/construct/provenance closure can only reduce the failed upper bound |
| `PD-MVR-B0 exact bridge` | Replace R18 identity overlap with an exact PDB-sequence-accession-ligand BindingDB–BioLiP2 bridge audit | 348 R18 identity links reduced to 23 exact candidates/5 targets/6 ligands/6 PDBs; optimistic provenance-conflict ceiling 5 and only 2 known TRAIN homology components | `STOP_PDMVR_B0_EXACT_BRIDGE_UPPER_BOUND_INADEQUATE`; no coordinate, bridge-rank, affinity, or missing-view training |
| `PD-MVR current formulation` | Determine whether its exact bridge identifies a contact-affinity common latent | The bridge loss contains contact and deployable pair representation but no bridge affinity `Y`; permuting bridge `Y` leaves objective and gradients unchanged, and the proposed rank checks `C <-> X` rather than `C <-> Y` | Stop on mathematical non-identifiability independently of the already failed bridge scale |
| `BioLiP26-D0` | Audit the official 2026 LLM/rule layer as a strict new-matrix source while byte-skipping four affinity fields | 23,502 rows yielded 17,679 strict UniProt-InChIKey-PubMed triples but only 258 new exact TRAIN-matrix edges; 18 components, largest `88.76%`, resource ceiling 91, greedy packing 59, and no row evidence/model/rule version | `STOP_BL26_D0_STRICT_MATRIX_ADOPTION`; retain metadata only and do not decode affinity |
| `AdaMBind primary applicability` | Determine whether target-wise MAML addresses strict zero-support target-and-scaffold cold prediction | AdaMBind gives each test target 5 or 40 labeled support pairs; its 40%-identity target split does not close scaffold, chemical-neighbour, provenance, assay, or endpoint axes, and its authors leave no-support true zero-shot as future work | Stop as the primary task; retain only a separately reported future `k=5` adaptation branch |
| `Virtual Binding System causal claim` | Interpret `z0(target)+delta(target,ligand)` as an observed ligand intervention | Only bound complexes are observed; for any target shift `h`, `z0'=z0+h` and `delta'=delta-h` produce the same bound output, unlike control/treated virtual-cell data | Stop causal/perturbation-world-model wording; a conditional bound-state encoder remains testable |
| `UBSE-G0R` | Test whether existing BioLiP binding-residue labels are independently repeatable and ligand-conditioned before model or coordinates | 1,028 cross-PubMed/cross-PDB sequence/connectivity units over 800 targets; correct contact Jaccard `0.75`, hard wrong-ligand `0.50`, margin bootstrap interval `[0.125,0.200]`, contact-label Recall@1 `0.5604` vs random `0.0846` | Pass only to stricter same-scaffold/additive-null topology; observed-label retrieval is not a deployable student |
| `UBSE-S0P initial whole-component split` | Force all same-scaffold panels into five whole homology-scaffold-PubMed conflict-component folds | 3,281 raw panels reduced to 2,170 after generic-scaffold closure, but the largest component held `58.3410%`, fit retained only 904 panels, and largest homology share was `5.3456%` | `STOP_UBSE_S0P_SAME_SCAFFOLD_INDEPENDENT_TOPOLOGY_INADEQUATE`; do not claim a population-wide balanced five-fold source |
| `UBSE-G0P` | Enumerate homology-scaffold-PubMed-closed same-scaffold multi-ligand panels without reading contacts | 1,612 panels/4,691 contrasts and packing 452 passed scale, but largest conflict component `28.0397%` and largest homology share `6.7618%` failed frozen caps | Formal STOP; permit one removal-only overrepresented-homology correction with unchanged gates |
| `UBSE-G0PB` | Remove once every pre-existing homology block above 5% and test a conflict-free pilot audit manifest | Removing 2 blocks/200 panels left 1,412 panels, 492 homology components, largest component `17.6346%`, packing 450, frozen audit 88, and closed residual train 1,324 panels | Pass only the selected 88-panel pilot to an affinity-blind G1 student with identical-parameter additive exact null; S0P's population-wide five-fold claim remains stopped |
| `UBSE-G1` | Predict held-domain within-panel residue-contact residuals from frozen ESM2 sequence and 2D ligand covariates beyond an exact additive null | Substrate passed with 1,138/57/81 contrast fit/validation/audit panels, but cross directional accuracy was `0.5107`, cosine `0.0689`, assignment `0.5201`, cross-minus-null delta `0.0219` (95% interval `[-0.0513,0.0952]`), and the protein-free-position control was stronger | `STOP_UBSE_G1_NO_DEPLOYABLE_INTERACTION_RESIDUAL`; no G2, affinity load, Stage-2 fitting, confirmation, or sealed access |
| `UBSE-A0 / A0C 3D event source` | Build and then independently correct a label- and affinity-blind coordinate locator manifest for a real residue-functional-group event teacher | A0 found 3,467 BioLiP file instances over 2,833 PDBs with zero cross-role PDB overlap, but incorrectly equated filename serial with mmCIF `auth_seq_id`; A0C byte-skipped affinity columns, uniquely recovered BioLiP column 20 for 3,467/3,467 rows, found only 39 serial/sequence-ID equalities, and produced zero duplicate corrected locators | Original mmCIF-locator claim withdrawn; `FREEZE_UBSE_A0C_LOCATOR_V2_KEEP_A0_REMOTE_WAIT`. Corrected metadata passes, but remote coordinates/events and A1 remain locked |
| `UBSE-A0W PLINDER weak-teacher audit` | Determine whether local PLINDER interaction strings can replace the pending residue-by-functional-group-by-event 3D teacher | Any-candidate coverage was 2,897/123/168 fit/validation/audit rows, but only 2,587/98/140 had a non-empty residue-by-event signature shared by every candidate; some rows had 93 candidates, and the stored field has no ligand atom, functional-group, exact BioLiP ligand-instance, or PubMed axis | `STOP_UBSE_A0W_AS_A1_REPLACEMENT`; retain only a separately preregistered A0-fit-only residue-event auxiliary probe under new source closure, with no affinity unlock |
| `UBSE-P0A execution attempt 1` | Run the frozen three-seed CUDA target-marginal anchor gate | Process remained at 100% GPU after 57.3 minutes, then disappeared before the next check with no checkpoint, ledger, result, traceback, or Windows crash event; the foreground execution session was also unavailable | `ABORT_UBSE_P0A_ATTEMPT1_EXECUTION_SESSION_LOSS_NO_SCIENTIFIC_DECISION`; rerun identical code detached with persistent logs, without changing any scientific gate |
| `UBSE-P0A corrected decision` | Learn a source-closed protein-only target-marginal pocket proposal before the typed 3D event student | On 62,849 rows/38,781 targets, corrected validation AP/AUROC/top-k recall were `0.3159/0.8398/0.2841`; AP exceeded propensity by `0.2147` (LCB95 `0.1604`) and fixed-seed shuffle by `0.2551` (LCB95 `0.1971`); three-seed AP range `0.0040`; CUDA utilization mean `98.40%` | `FREEZE_UBSE_P0A_FOR_A1_POCKET_PROPOSAL_ONLY`; ranking proposal only, with no ligand-conditioned evidence, calibration claim, affinity load, confirmation, or sealed access |
| `UBSE-A1-v1 coupling design` | Use unbalanced OT over target and ligand marginals to identify typed residue-FG placement | UOT can change pair-conditioned row/column marginals and win without coupling; the v1 rectangle collapses the FG axis, its dustbin does not exactly close, and the previously read G1 audit is not an independent confirmation set | `REVISE_UBSE_A1_BEFORE_PREREGISTRATION`; require a pair-conditioned rank-one null, explicit dustbin, within-complex FG checkerboard, stronger split/membership closure, and a fresh confirmation role |
| `UBSE-A0D remote availability` | Test every corrected A0C v2 official RCSB coordinate URL with HEAD-only metadata requests | 2,833/2,833 unique URLs returned HTTP 200 `application/gzip`: fit 2,496/2,496, validation 140/140, audit 197/197; one URL needed a second attempt; zero response-body bytes and zero origin violations | `REQUEST_UBSE_A1V2_SOURCE_AND_TOPOLOGY_PREREGISTRATION`; remove only the remote-addressability WAIT, while coordinate parsing/events/affinity remain locked |
| `UBSE interaction-space / episodic advice refresh` | Reassess UCE-style universal interaction embeddings, AdaMBind, and virtual-cell world-model language against the completed UBSE evidence | Pair-level interaction tokens and strict episodes are compatible with A1-v2, but bound-only BioLiP does not identify an unbound-to-bound causal transition; AdaMBind requires 5 or 40 same-target labels; universal/pooled embeddings and meta-learning cannot create missing residue-FG correspondence, source independence, or checkerboard topology | Retain only a conditional bound-state interaction representation plus training-time dual-cold episodes. Keep the primary task at `k=0`; reserve constrained `k=5` adaptation for a later secondary task after coupling passes; no affinity, confirmation, or Stage-2 unlock |
| `UBSE-A1-v2 source-role certificate` | Freeze label-blind A1-R/A1-S/A1-C roles and independently bind retained homology, chemical-neighbour, locator, dependency, and inherited-model membership | A1-R 153 targets/459 instances; A1-C 512+64; A1-S 1260/59/81 panels. Retained A1-S cross-role containment maxima `0.3806/0.3511/0.0781` and ECFP4 maxima `0.4694/0.4932/0.4107` have zero conflicts. Locator identities are unique, but A1-R forms 124 PDB/PubMed/physical-ligand components with worst-case Kish `n_eff=98.77`; current P0A overlaps 153/153 and 576/576 targets | Freeze deterministic metadata certificate and component-bootstrap correction. `SR-4` remains failed; retrain P0A-v2 and keep event/coordinate-body access locked |
| `UBSE-A1-v2 locator / H0V` | Recover distinct BioLiP filename serial and mmCIF auth-sequence IDs, then verify official URL availability without coordinate access | 1,035/1,035 strict unique joins, 997 unique URLs, and 997/997 HTTP 200 `application/gzip`. Corrected H0V followed zero redirects, never iterated a body stream, and recorded zero actual downloaded bytes | `FREEZE_A1V2_H0V_STRICT_HEAD_FIREWALL_KEEP_BODIES_LOCKED`; locator half passes, extractor/FG/assembly half of SR-5 remains pending |
| `UBSE-P0A-v2 prelaunch review` | Decide whether the A1-C-closed target-marginal proposal is ready for a 5.5-hour three-seed CUDA run | Safe counts reproduce 54,868 rows/32,769 targets/11,126,109 residues/33,216 windows and 13 focused tests pass, but three independent reviews found false case-sensitive PDB closure (655 case-fold matches), hard-coded homology zero, count-only label membership, no real preregistration binding, nonrecoverable seed state, incomplete ledger/heartbeat gates, and no worst-batch CUDA smoke; old peak memory was 7,945/8,188 MiB | `NO_GO_UBSE_P0A_V2_LONG_TRAINING_PRELAUNCH_AUDIT`; no v2 training product exists. Reopen only after all closure, manifest, recovery, telemetry, and memory gates pass and A1-direct makes the proposal scientifically necessary |
| `UBSE-A1-v2 finite DTA bridge` | Prevent A1 from becoming an open-ended structural side project by specifying its exact path into strict dual-cold affinity | Four dependent stages are frozen: construction/topology, teacher coupling beyond a pair-conditioned rank-one null, a deployment-side no-P0A student with a 224-D coupling-only coordinate, and a minimal frozen-student head `B0 + theta^T z_int` with exact null `theta=0`; no raw target bypass is permitted | `NO_DECISION_A1V2_STAGE1_EXTRACTOR_TOPOLOGY_PENDING`; only the Stage-1 extractor/FG/assembly/component-power/coordinate-body preregistration is next. Events, student, affinity, confirmation, and sealed outcomes remain locked |
| `UBSE-A1-v2 flexible-kernel amendment` | Test whether the LOCK/CLOCK principle supplies a bounded estimator improvement without pretending that it creates A1 information | Official arXiv v1 confirms LOCK reverts from a local nonlinear term to an extra linear predictor and CLOCK learns structure-conditioned amino-acid correlations across mutation landscapes. A1 has no substitution-response landscape; sample-wise `B B^T` does not define a cross-sample PSD kernel; a linear-plus-local residual does not return to B0 at distance | `REVISE_AND_DEFER_A1V2_FLEXIBLE_KERNEL_TO_STAGE4B`; preserve Stage-4A `B0 + theta^T z_int`. Only after Stages 1-3, an affinity-blind PSD/semantic gate, and Stage-4A PASS may a frozen Nyström local block be executed. Do not claim ligand-conditioned CLOCK, calibrated OOD uncertainty, or new independent samples |
| `Flexible-kernel review firewall incident` | Preserve independence of the affinity-blind design review | One computational reviewer used an overbroad local search and received snippets from a historical affinity-bearing JSON artifact. Its numerical output was not used, the review vote was discarded, and a no-context replacement reviewer used a strict document whitelist | The contaminated review cannot support a preregistration or independent-audit claim. No further value access is authorized; confirmation and sealed permissions remain unchanged |

These failures do not support the claim that protein structure is generally
uninformative. They show that the current PLINDER absolute mechanism state is
outside most of the ChEMBL strict-dual-cold support and that increasingly
expressive fusion heads learn to suppress it. Reopening requires a different
information object or a new source, not another router/attention block.

The subsequent RDIB/PD-MVR review tested the requested exact-complex
alternatives rather than inventing another representation. Exact replicated
differences and exact dual-modal bridges both failed optimistic upper bounds.
The recurring-directed-edit audit also failed after correcting the unit from
globally recurring edits to same-target cross-PubMed replication. The
official BioLiP-2026 LLM/rule source subsequently failed strict new-matrix
topology and evidence-version auditability.

The later UBSE advice introduced a genuinely different observation object:
the absolute ligand-conditioned binding-residue state rather than a replicated
ligand-pair difference. Its cross-publication reliability gate passed, and a
label-blind same-scaffold panel program survived after one preregistered
removal-only homology correction. This does not reverse the RDIB failures:
UBSE-G1 then showed that sequence and 2D chemistry did not predict centered
contact residuals beyond an identical-parameter additive null on the frozen
audit. The exact-null, ligand-destruction, protein-position-destruction, and
pair-assignment gates all failed. The program therefore returns to new
information-object/source research or the frozen prospective cycle-closing
acquisition design; another post-hoc fusion-head rescue is not authorized.

The two independent, user-requested label-blind follow-up routes were then
executed:

- `SISMT`: test whether the intersection of target-supported covariance
  directions and privileged destruction-certified source directions has
  nonzero stable dimension and adequate effective sample size.
- `DTIOD`: test whether local segment-by-substructure mixed finite-difference
  responses have stronger privileged semantics and target support than the
  absolute R6 moment.

Neither route is allowed to fit ChEMBL affinity until its own feasibility gate
passes. All GPU-capable teacher inference, spectral algebra, distillation,
and model fitting use the CUDA-enabled `drug` environment.

SISMT had no stable retained direction. DTIOD found a real privileged local
response, but it failed the mandatory ligand-pair destruction: the response
was protein-local rather than an identified target-ligand operator. The
subsequent KLIFS/R17 bridge also failed its rectangular support gate. The
BindingDB/R18 replacement then failed four of five frozen source gates despite
6,561 PDB identifiers, because those structures covered only 340 independent
ligands. The current PLINDER/KLIFS/BindingDB source program is therefore
stopped for absolute-state, tangent, and exact-bridge transfer. Repeated
structures cannot substitute for target-domain chemical and target support.
The subsequent R19/R20 ChEMBL rectangle audits also failed after concentration,
endpoint, and provenance control. The binding successor is now source
recovery for real same-target paired-complex deltas (RBSDD), beginning with
PBCNet2.0-D0.

### 0.1 Current decision

| Item | Current state |
| --- | --- |
| Final category | **3: current data do not identify transferable amino-acid substitution geometry; a new source-resolved substrate or prospective measurements are required** |
| Active program | OpenMut data recovery and mutation-reordering identifiability |
| Formal gate | `OPENMUT_SOURCE_AND_POWER_AUDIT__NO_TRAIN` |
| `OMUT-D0` | **Complete 2026-07-28**, `OMUT_D0_SOURCE_FREEZE_COMPLETE`; 16 sources frozen, 25 reads, 0 rows materialized, 0 firewall violations; anchor `1953cff4c1d7301c51d1ef934c0c5c913f7c154022ad515c53e416bbce8f82f9` |
| `OMUT-F0` | **Complete 2026-07-28**, `OMUT_F0_DAVIS_ROLE_FROZEN__PRESERVE_CONFIRMATION`; `panel_davis` stays sealed and unconsumed, DAVIS-Complete WT-mutant contrast foreclosed, 0 rows read |
| `OMUT-D1` | **Complete 2026-07-28**, `OMUT_D1_TOPOLOGY_ADEQUATE`; 62 `k=4` candidate mutation components, 17 accessions, concentrated in HIV Gag-Pol/ABL1/EGFR |
| `OMUT-X0` | **Active** |
| Later stages | Blocked |
| Tau transfer | Conceptual supervision-density principle passes; bilateral intervention topology and teacher admission both stop; no training |
| Open-evidence pretraining | Source recovery is admissible as D0 metadata work; real interaction-pretraining hypotheses remain untested and unauthorized |
| `MEDIP-S0` | `MEDIP_S0_ENGINEERING_CALIBRATION_STOP`; synthetic recovery/calibration passed, but metadata and selectivity destruction did not |
| New affinity values in this round | None read |
| Davis target-conditioned confirmation | Not consumed; a historical arm-blind power audit did read base-arm power labels |
| ChEMBL confirmation | Five rows were displayed during an earlier CROSSDOC schema inspection; permanently record `confirmation_labels_read=true`, although those rows entered no experiment |
| Sealed test | Never consumed |
| Currently authorized real-outcome model training | **None** |
| Runtime | The `drug` environment is operational; hardware is not the current scientific bottleneck |

Only D0 metadata work is currently authorized on real sources. MEDIP-S0 was an independent
generated-data estimator falsification and has stopped; it grants no permission to open outcomes.
Historical TRAIN or development analyses in this ledger are evidence about prior work, not
permission to reopen their labels or train a predictor.

### 0.2 The four claims must be proved in order

| Level | Claim | Current evidence | State |
| ---: | --- | --- | --- |
| 1 | Real signal exists in a local protein landscape or biological information source | Local LOCK/CLOCK results in their original task; residue teacher, GO/compatibility pretexts, and KirHub mutation double differences here | Established only in those local tasks |
| 2 | Substitution geometry is identifiable in this DTA setting beyond identity, family, composition, and random geometry | Fixed LOCK did not beat aligned identity, BLOSUM-label permutation, or matched random PSD | Failed |
| 3 | The coordinate predicts target-specific ligand reordering | KirHub and CROSSDOC show that some reordering exists; WTPAIR and LOCK semantic gates failed | Biological indication exists; coordinate prediction is not established |
| 4 | The mechanism improves strict target-cold plus ligand/scaffold-cold DTA | No compliant real substrate and no authorized predictive training | Untested and not claimable |

Passing one level never proves the next. In particular:

1. LOCK or CLOCK working on local mutation landscapes does not establish BLOSUM or CLOCK in
   dual-cold DTA.
2. Correct synthetic estimator calibration does not establish that the real mechanism exists.
3. Observed target-specific reordering does not identify a protein coordinate that transfers to
   unseen targets and unseen scaffolds.

### 0.3 Evidence language

| Term | Strict meaning in this project |
| --- | --- |
| `PASS` | Passes only the preregistered gate at that stage; it does not authorize a higher claim |
| `STOP` | Stops the route under the tested information conditions; only explicitly new information can reopen it |
| `NO_DECISION` | The estimator or execution conditions do not identify the effect; it is not evidence of a zero effect |
| Engineering pass | Software, data access, or synthetic calibration works; no biological or predictive credit |
| Mechanism evidence | Correct semantics beat matched destructive controls |
| Performance evidence | Under a strict firewall, the method beats the strongest nested null by a material, adequately powered effect |
| Independent `n` | Paired complete-case biological units after joint, transitive blocking over base protein, homology/family, and provenance, or an explicitly preregistered multiway-clustered unit; rows, pairs, quartets, folds, seeds, and technical replicates are descriptive only |

---

## 1. Core Dilemma: No Real Substrate Yet Identifies Dual-Cold Interaction Geometry

### 1.1 This is not ordinary affinity regression

An observed value can be decomposed as:

```text
y_obs(t,l,e,a,d,s)
  = mu
  + alpha_target(t)
  + beta_ligand(l)
  + h_interaction(t,l)
  + q_endpoint(e)
  + q_assay(a)
  + q_document(d)
  + q_source(s)
  + epsilon
```

Random splits can exploit ligand potency, target identity, family, document, and assay recurrence.
Strict dual-cold prediction instead asks for `h_interaction(t*,l*)` when both target `t*` and
ligand/scaffold `l*` are unseen. This is extrapolation of an interaction function, not random
completion of a partly observed matrix.

Earlier cold-target work exposed the first obstacle: measurement-noise correction can improve the
measurements without predicting an unseen target's affinity baseline. Dual-cold prediction is harder
because unseen ligand chemistry and its non-additive interaction with an unseen target must both
transfer. Cleaner labels, a deeper network, or a successful pretext task is therefore insufficient.

### 1.2 The estimand is ligand reordering, not absolute potency or taxonomy

For a wild type, a single mutant `m`, and shared ligands `a,b`, the central estimand is the four-cell
mixed difference:

```text
C(m;a,b)
  = y(mut,a) - y(mut,b) - y(WT,a) + y(WT,b)
  = [y(mut,a)-y(WT,a)] - [y(mut,b)-y(WT,b)].
```

Primary rectangles require one endpoint and the same assay/context and source conditions across all
four cells. Under that restriction, target and ligand main effects and common measurement offsets
cancel, leaving mutation-specific ligand reordering. Cross-assay rectangles are sensitivity-only.

The algebra is valid but not new. REWIRE, MISO/MISO-OR, DICE/AXIS, KirHub DD, WTPAIR, and R-MAON
already cover mixed differences, additive-nuisance removal, antisymmetric bilinear carriers, or a
direct centered operator. The genuinely untested information object is narrower:

> A canonical WT-to-single-substitution operation jointly acting with a recurring directed ligand
> chemical change, measured in continuous Ki or Kd complete rectangles across independent sources.

Renaming the same four cells as a quartet, listwise loss, Hodge decomposition, functional ANOVA,
Delta2Rank, or OMRO creates neither information nor independent samples.

### 1.3 Transfer requires an identified transferable coordinate

The minimum transferable hypothesis is:

```text
h_interaction(t,l) ~= u(t)^T Theta phi(l).
```

The hard part is not fitting `Theta`; it is demonstrating that:

- `u(t)` is not merely family, identity, composition, study depth, or source;
- `phi(l)` is not merely generic potency, scaffold identity, or library composition;
- the correct `u(t)` and `phi(l)` pairing beats wrong-target, wrong-mutation, ligand-only,
  source-only, and matched-random controls;
- the advantage survives joint isolation of target, homology, family, ligand, scaffold,
  chemical-neighbor, assay, document, and provenance components.

Many functions fit a kinase-only, single-source, or ligand-warm graph but extrapolate in opposite
directions on an unseen family. If two data-generating processes agree on every observed kinase unit
but imply opposite unseen-family ligand orderings, the training likelihood cannot choose between
them without cross-family measurements or an independently validated target coordinate. This is an
identifiability limit, not an optimizer limit.

### 1.4 Complete rectangles are not many independent observations

With `L` shared ligands, one target or mutation has only `L-1` centered-response degrees of freedom.
Expanding them into `L(L-1)/2` ligand pairs or more quartets is correlated algebraic replication:

- `L=92` gives exactly 4,186 unordered pairs or 8,372 directed pairs, not tens of thousands of
  independent target units, and it cannot turn 34 genes into 4,186 genes;
- about 2.81 million source labels can expand into about 8.6 million PBCNet2.0 pairs without adding
  biological independence;
- folds test split generalization and seeds test training stability; neither adds an experimental
  source;
- BindingDB, ChEMBL, PLATINUM, MdrDB, and supplements that trace to one DOI or experiment form one
  provenance lineage.

Power is calculated on a paired component contrast:

```text
MDE80 = (z_0.975 + z_0.80) * SD(delta_component) / sqrt(n_component)
      = 2.8016 * SD(delta_component) / sqrt(n_component).
```

At the optimistic planning value `SD=0.10`:

| Independent units | Exact/rounded MDE80 | Interpretation |
| ---: | ---: | --- |
| 6 WT-mutant base-protein pairs | `0.1143748` (`0.114`) | The 384-cell A0 can estimate reliability and variance only |
| 10 / 11 base proteins | `0.089 / 0.085` | DAVIS-Complete alone is underpowered; these are projected, not empirical, MDEs |
| 22 families | `0.060` | KirHub family-level inference cannot stably resolve `0.03` |
| 34 genes | `0.048` | KirHub gene-level inference remains insufficient |
| 88 components | `0.0298652` (`0.030`) | Optimistic mechanism floor before clustering, missingness, multiplicity, or model fitting |

The approximately 423 independent multi-family components and at least 40 randomized,
scaffold-diverse query ligands per target belong to strict predictive T1 planning. They are not the
initial mechanism gate and must not be conflated with the optimistic 88-unit mechanism floor.

### 1.5 Model capacity cannot repair missing identification

A larger model cannot:

- manufacture missing protein families, independent provenance lineages, or complete rectangles;
- determine whether BLOSUM semantics, identity, or matched random PSD is the load-bearing geometry;
- turn duplicated database copies of one paper into independent replication;
- make IC50, percent inhibition, Kd, and Ki one interchangeable endpoint;
- turn pair or quartet expansion into new biological units;
- repair a factorization whose direction parameters are unidentified under the null;
- prove that a synthetically generated mechanism exists in real data.

Capacity instead increases the ability to exploit taxonomy, ligand potency, assay/document identity,
known poses, and dataset membership. The project therefore uses falsification first: identify data
and coordinate information before authorizing the smallest model. A failed low-cost gate cannot be
rescued by bandwidth, rank, backbone, seed count, or training duration.

### 1.6 Exact statement of the bottleneck

> The binding constraint is neither "biological information is useless" nor "the model is too
> small." It is the absence of an adequately powered, source-resolved, multi-family,
> endpoint-consistent, complete or correctly weighted factorial measurement substrate that can
> identify target-specific interaction geometry when both target and ligand are cold.

Three non-substitutable gates follow:

1. **Data:** legal, independent, sufficiently deep target-by-ligand rectangles exist.
2. **Semantics:** the correct substitution coordinate beats family, identity, composition,
   permuted BLOSUM, matched random PSD, and wrong-target controls.
3. **Prediction:** the coordinate that passes semantics beats the strongest nested null in strict
   dual-cold, cross-family, cross-provenance evaluation.

No success at one gate compensates for failure at another.

---

## 2. Detailed Data, Representation, Estimation, and Evaluation Dilemmas

### 2.1 Target extrapolation and the affinity baseline

Measurement-process de-noising reduced held-out measurement error from `0.606` to `0.460`; real
assay offsets improved it by `24%` and a shuffle control by only `5%`. Downstream cold-target
macro-Spearman nevertheless changed by `-0.051`, worse in all three seeds. Other examples show the
same estimand mismatch:

- pair-compatibility retrieval AUC was `0.635` versus `0.512` for broken pairs, but downstream `k0`
  RMSE was `+0.0234` worse than random initialization in all three seeds;
- residue-teacher KL improved from `2.5366` to `0.0742`, with shuffle `0.2609`, while RMSE at
  `k=0/4/16` worsened by `+0.0651/+0.0476/+0.0412`;
- true GO hierarchy beat shuffled labels but lost to random initialization in multi-seed evaluation;
- a NARD single-seed improvement did not reproduce: three-seed RMSE difference
  `+0.0035 +/- 0.081` with mixed signs.

These tasks contain signal, but not necessarily the signal required for strict affinity
extrapolation.

### 2.2 Real substrates collapse under rectangle requirements

A compliant primary rectangle requires:

- explicit, verifiable WT and single-mutant constructs;
- the same ligand parent measured against WT and mutant;
- one endpoint, with exact Ki and exact Kd analyzed separately;
- the same assay/context for all four primary cells;
- traceability to assay, document, original paper, institution, and experimental lineage;
- enough shared, scaffold-diverse ligands;
- multiple independent base proteins after joint homology/family/provenance blocking.

DAVIS-Complete is ligand-dense but has only about 10 canonical WT base proteins, or 11 by the paper's
broader naming convention. KirHub has 222 eligible constructs but only 34 genes, 22 families, and
one source. CROSSDOC's strongest strict signal has only 11 homology components. Database row count
is not effective sample size.

### 2.3 Power, pseudoreplication, and unstable statistics

Repeated failure modes include treating pairs, quartets, folds, seeds, or technical replicates as
independent; a ligand giant component joining many scaffolds/documents; very large design effects
from homology and provenance clustering; and target-level ratios dominated by near-zero
denominators.

The PARC/HQ-GBMA `inside/signal` containment ratio illustrates the last problem. The smallest positive
denominator was `0.00309`, the median was `0.1473`, a roughly 48-fold span, and one component reached
`56.09`; deleting one target could reverse the mean contrast. Robust reanalysis established only:

- pooled ESM minus shared-global: median `-0.1165`, sign-test `p=0.0013`, Wilcoxon `p=0.0003`;
- PARC minus shared-global: median `-0.0243`, sign-test `p=0.068`, Wilcoxon `p=0.079`, unresolved;
- PARC minus random positions: median `+0.0317`, Wilcoxon `p=0.0171`, but sign-test `p=0.254` and
  Holm-adjusted `p=0.137`, therefore not evidence;
- PARC minus pocket composition: median `+0.0079`, Wilcoxon `p=0.857`.

The original mean is not a valid effect-size estimate. Future ratio statistics require
ratio-of-sums, denominator reporting, and component-blocked inference. For this unusually
high-variance containment estimand, detecting `0.03` needs about 618 independent components and
detecting `0.02` about 1,391; only 77 were evaluable. Those figures are specific to this ratio and do
not replace the general OpenMut 88-unit optimistic mechanism floor or approximately 423-unit strict
prediction plan.

Query depth also matters. Metz's same-arm retraining MDE was `0.0181`, but its cross-arm
heterogeneity reference was `0.0614`. Davis had about 102 components and median query depth 12; its
raw retraining MDE was `0.0688`, and the frozen `2.32` multiplier gave
`0.0688 x 2.32 = 0.1596`. Reinecke had median query depth 5 and MDE80 `0.0668`. Retraining variance
within one model arm cannot replace real paired-arm dispersion.

### 2.4 Family, document, assay, and provenance leakage

In WTPAIR, the KLIFS-group centroid scored `0.0911`, far above pooled-ESM WTPAIR at `0.0323`.
Taxonomy can therefore predict an average homologous-kinase profile without identifying
within-family substitution geometry.

RECRO initially found cross-document residual correlation `+0.3345`, but:

- `79.75%` of co-measured target-ligand cells were nearly exact numeric duplicates;
- `91.6%` of nominal cross-document comparisons remained in one provenance family;
- document-family-disjoint residual correlation fell to `+0.0901 [-0.0557,+0.2358]`;
- correct-target minus matched-wrong-target was `-0.0744`.

The decision `RECRO_SIGNAL_EXPLAINED_BY_PROVENANCE` rejects ChEMBL document IDs as independent
experimental environments. Database, assay, plate, and table identifiers are not independent
sources by default.

### 2.5 Endpoint incompatibility, censoring, and selection

- Primary labels are direct biochemical exact Ki or exact Kd, analyzed in separate strata.
- IC50, censored values, multi-mutations, indels, phosphorylation, and cross-assay rectangles are
  sensitivity-only.
- Cell viability, clinical outcome, qualitative resistance, docking, and simulated scores are
  excluded.
- Inactive, censored, failed-curve, BQL/AQL, and out-of-quantification outcomes remain in the
  registry; hit-only selection is forbidden.
- Any outcome-dependent rectangle sampling requires frozen, nonzero inclusion probabilities and
  weighting before value access.

KirHub's aggregate 1-uM residual activity supports an ordinal reordering check but is not a Ki/Kd
training label. MdrDB mixes direct binding, cellular resistance, predicted structures, and
simulation, so it can index papers but cannot supply bulk labels. Papyrus, BindingDB, and ChEMBL
aggregations cannot support confirmation without reconstruction of endpoint, construct, and
provenance.

### 2.6 Ligand-only, scaffold, and chemical-edit shortcuts

Generic ligand potency is the only signal that remained strong under provenance-family isolation:
RECRO raw family-disjoint correlation was about `+0.700`. Every protein-conditioned model must beat
both ligand-only `B0` and an identical shared-global interaction null, not merely random prediction.

Local chemical edits are not automatically informative:

- MMP-dense public kinase data such as PKIS are often percent inhibition, while continuous-pK panels
  are MMP-sparse;
- historical MMP/local-edit features lost to whole-molecule differences in chemistry-cold tests;
- DICE RAW-CHEM was only `+0.00119 [-0.00815,+0.01006]`, and target-disjoint transfer failed;
- the label-free conformer topology had a ligand giant component containing `74.37%` of ligands and
  `80.82%` of rows and lacked provenance-family metadata.

MMP, Morgan difference, physicochemical delta, whole-molecule difference, and conformer ensembles
must therefore be separate coordinate arms or controls. Concatenating them and attributing a
capacity gain to mechanism is prohibited.

### 2.7 Protein coordinates and semantic non-identification

#### Pooled ESM-2

Pooled ESM-2 is substrate-dependent, not universally null. In corrected KirHub A1,
component-macro Spearman was `0.0719`, above ligand-only `0.0429`, shuffle `0.0261`, and matched
random protein `0.0200`. However, `ESM - ligand-only = +0.0290` missed the frozen `+0.030` threshold,
and `ESM - group centroid = -0.0110`. In WTPAIR, pooled ESM trained directly on protein-pair by
ligand-pair mixed differences also failed strong controls. In Metz/PARC, it lost to shared-global
and was indistinguishable from exposure-matched derangement (`p=0.82`).

The supported diagnosis is that whole-sequence pooling expresses coarse homology, family, or global
semantics but does not reliably carry localized ligand-conditional reordering on these substrates.
This does not mean ESM has no protein information or that every PLM is invalid.

#### Aligned identity, position, and composition

Aligned KLIFS pocket identity beat pocket shuffle, matched random target, and pooled ESM in a
ligand-warm kinase oracle, so correct position and target assignment contain information. PARC,
however, added only median `+0.0079` over composition with Wilcoxon `p=0.857` in 77 evaluable
components. The position-specific increment is unidentified; pocket biology is not refuted.

KLIFS is kinase-only and can encode family, structural-template quality, and target identity.
Cross-family use requires an independently frozen alignment/structure map plus position shuffle,
structure shuffle, within-family wrong-target, and composition-only controls.

#### Fixed BLOSUM-derived LOCK

The fixed LOCK label-free geometry is valid and nontrivial:

| Label-free quantity | Result | Gate |
| --- | ---: | --- |
| Normalized kernel minimum eigenvalue | `2.07e-16` | PSD pass |
| Family/composition residual energy | `0.3851` | Non-redundancy pass |
| CKA with pooled ESM-2 | `0.1633` | Not an ESM copy |
| Within-family nonconstant pair fraction | `0.9992` | Nonconstant pass |

Its effective rank was about `289.36`, while the top 16 dimensions retained only `0.2340` energy, so
it is not an established low-dimensional coordinate. On 353 activity-eligible genes, 303 frozen
target components, and 92 shared ligands, all coordinates entered one fixed top-eight
squared-similarity estimator:

| Contrast | Mean | 95% interval | Interpretation |
| --- | ---: | --- | --- |
| LOCK - group centroid | `+0.02997` | `[+0.01223,+0.04799]` | Missed the `0.030` gate by `0.0000257` |
| LOCK - aligned identity | `-0.01370` | `[-0.02304,-0.00459]` | Exact identity was stronger |
| LOCK - composition | `+0.03565` | `[+0.02081,+0.05087]` | Beat composition |
| LOCK - pooled ESM-2 | `+0.04548` | `[+0.02826,+0.06271]` | Beat pooled ESM |
| LOCK - position shuffle | `+0.08327` | `[+0.06315,+0.10300]` | Alignment position mattered |
| LOCK - sequence shuffle | `+0.09184` | `[+0.07168,+0.11210]` | Correct sequence mattered |
| LOCK - BLOSUM permutation | `-0.00161` | `[-0.01138,+0.00840]` | BLOSUM token semantics were not identified |
| LOCK - matched random PSD | `-0.00444` | `[-0.01553,+0.00684]` | Did not beat matched random geometry |
| LOCK - matched wrong target | `+0.09403` | `[+0.07372,+0.11534]` | Correct target assignment mattered |

The precise decision is that alignment, correct sequence, and correct target carry information, but
fixed BLOSUM substitution semantics are not identified as the cause. Rounding `+0.02997` to a pass
or generalizing the failure to substitution biology is prohibited.

#### Conservation_LOCK and true CLOCK

`conservation_LOCK - fixed LOCK` was `+0.00664 [+0.00125,+0.01219]`, and
`conservation_LOCK - composition` was `+0.04229 [+0.02703,+0.05756]`; its group-mode mean
`0.45374` remained below aligned identity `0.46079`. This adds a label-free positional conservation
prior, not a structure embedding or local mutation-landscape likelihood. It is not CLOCK, does not
establish low dimension, and cannot rescue fixed LOCK under the preregistration.

True CLOCK requires an externally frozen position-specific, structure-conditioned map,
leave-family-out mapping, and a parameter-matched structure-shuffled control. Current targets
usually have one WT sequence and no local mutation landscape. CLOCK is untested, not failed.

### 2.8 Estimator and null identifiability

#### Overparameterized target maps

HQ-GBMA/PARC `ProteinGrassmann(32 -> 64, r=6)` had 27,072 parameters for about 89 training targets,
or roughly 304 parameters per target. A shared-global Grassmann subspace had about 348 manifold
degrees of freedom. A large target map can memorize family, noise, or target-specific coefficients
without extrapolating. Any unlocked target-specific channel is therefore restricted to a directly
centered operator or a very low-degree correction.

#### The K-LBP factorized null is nonregular

K-LBP used:

```text
Theta = gamma * a * c^T.
```

At `gamma=0`, every `a,c` pair defines the same shared-global null:

```text
d mu / d a = 0
d mu / d c = 0
d mu / d gamma depends on arbitrary a,c.
```

Direction parameters exist only under the alternative, so Fisher information is singular at the
null. Unit-norm constraints remove scale exchange away from zero but do not identify null
directions. R3 stopped at the first `S1_gamma_0.0` replicate; Variant E converged in only 3/5 folds
within the frozen 200 iterations. The decision was
`R3_ESTIMATOR_INSENSITIVE_NO_DECISION`. This is not evidence for `gamma=0`, and more iterations or
seeds do not repair the null geometry.

#### R-MAON repairs the estimator, not the missing substrate

R-MAON directly parameterized and centered `Theta`, making `Theta=0` the unique regular null. G0-B
used synthetic outcomes with empirical Metz TRAIN noise:

| Regime | Rejection/power | Median recovered/truth |
| --- | ---: | ---: |
| `S1_null` | `0.055` | - |
| `S1_active` | `1.000` | `1.0027` |
| `S2_heteroscedastic` | `1.000` | `0.9993` |
| `S3_degenerate_signal` | `0.995` | `0.9610` |
| `S5_null_coordinate` | `0.045` | - |

This shows that the estimator recovers a registered synthetic mechanism with the intended error
rate. It does not show that the mechanism exists. G0-A stopped because no prospective manifest
existed: `RMAON_G0_TOPOLOGY_OR_POWER_STOP`. The assay-monotone multi-fidelity proposal remains
untested, with no real Spearman, RMSE, calibration, or paired gain.

### 2.9 Structure, pose, and privileged external teachers

The physical-structure route stopped because inputs and estimands failed, not because of GPU limits:

- deployment docking rank correlation was about `0.005`;
- experimental holo versus AlphaFold increment on eight matched targets was `-0.019`;
- native-pose joint model scored `0.257`, below ligand marginal `0.394` and heavy-atom control
  `0.429`;
- adding physics reduced a descriptor baseline from about `0.381` to `0.125`, gain LCB `-0.395`.

These results close adding docking/contact/pose to the same observational graph as a rescue. They do
not reject structure coordinates on a genuinely new, compliant substrate.

| Candidate | Useful information | Why it is not an admissible current source |
| --- | --- | --- |
| PBCNet2.0 | Same-pocket ligand-pair relative-affinity teacher | About 8.6M pairs expand about 2.81M mixed-endpoint labels; requires a similar co-crystal, Glide/MCS poses, and restricted exact training membership |
| FLOWR.ROOT | Dense project-specific measurement plus LoRA adaptation | Zero-shot failed on four unseen proprietary SAR projects; useful adaptation depended on private dense labels and structural inputs |
| NISE/Boltz | Local protein design or compatibility confidence | Confidence or binary binding is not cross-target quantitative affinity or reordering |
| Neosurface | Ligand-conditioned protein surfaces are biologically meaningful | Partner retrieval/design with ligand-bound PDB input, not small-molecule DTA |
| ViT affinity | Known-complex structural baseline | PDBbind random 90/10; even the hardest subset allowed structural neighbors; not strict dual-cold |
| EnsembleEGNN | Conformer-set ligand encoder | Cyclic-peptide permeability without protein input; current data have a ligand giant component |
| Hybrid UQ | Calibration, rejection, and aleatoric/epistemic reporting | Released multitask path is ligand features plus target-specific outputs without target/provenance-disjoint validation |

Pose, support labels, future states, known complexes, or membership unavailable at deployment are
training-only privileged information. They require contamination audits and can never be relabeled
as measured affinity.

### 2.10 Tau future-supervision transfer

The paper *tau: Learning Touch-Augmented Vision-Language-Action Models from Future Visual
Supervision* was audited from frozen arXiv PDF `2607.24485v1`, SHA-256
`b9dd13ee0fc69ce091c6f10b738caaca05aef95608b317e06b2c8021d4b2be4b`.

Its actual data comprise four robot tasks with 100 trajectories per task, collected under one robot,
site, and acquisition system. Current tactile tokens plus a horizon of 32 actions predict changes in
future pretrained VLA visual latents. Future tactile-change magnitude weights the auxiliary loss,
which is removed at deployment. Reported mean full-task success was:

| Arm | Mean full-task success |
| --- | ---: |
| Full tau model | `71.25%` |
| Without action conditioning | `58.75%` |
| Without predictive self-supervision | `51.25%` |
| Without tactile encoder/adapter | `28.75%` |

Each model-task condition had 20 physical trials. The paper did not report task- or seed-cluster
confidence intervals. Generalization was limited to two of the four tasks, with two unseen objects
and two distractor sets per task; it was not a simultaneous new-task/new-object analogue of strict
dual-cold DTA.

Future frames are real additional within-trajectory sensor observations, so tau increases
within-trajectory supervision density and possibly constraint rank. Its time windows and future
offsets remain correlated. They do not increase the number of independent tasks, objects,
environments, sites, or provenance lineages. DTA pair or quartet expansion is weaker still: it
creates no new sensor observation and preserves the `L-1` centered degrees of freedom of the
original `L` affinity values.

The only transferable abstraction is:

```text
current interaction state + a verifiable intervention
    -> change in interaction state
```

The action in tau is temporally ordered and followed by a synchronized observed response. A
mutation string or medicinal-chemistry comparison becomes an analogous DTA intervention only when
the exact WT-to-single-substitution or directed ligand edit, endpoint, assay/context, inclusion
mechanism, and provenance are frozen. Most public comparisons are selectively observed static
contrasts, not interventions. DTA also has no free "future state": post-edit affinity is the original
label; sequence/SMILES embeddings are generally separable; docking/contact latents reopen the
failed privileged-structure route.

For a proposed frozen teacher `g(t,l)`, the required pair-specific quantity is:

```text
q_m(l)       = g(mut,l) - g(WT,l)
q_centered(l)= q_m(l) - mean_l q_m(l)
T_m(a,b)     = q_centered(a) - q_centered(b)
```

A pair-shaped input alone is insufficient. Before real labels, the teacher's interaction residual
after target-only plus ligand-only projection must be nondegenerate, fold-local, overlap-excludable,
and different under correct versus wrong/shuffled target, mutation, ligand, edit, pose, and
interface assignments. With labels, `T_m(a,b)` must predict the real `C(m;a,b)` under
base-protein/homology/family/provenance aggregation and beat separable, matched-2D, random,
source-only, and within-family wrong-target nulls.

The attachment's bilateral affinity differences, four-cell counterfactual energy, and path
consistency are therefore valid but not new modules. They overlap REWIRE, MISO, DICE/AXIS, WTPAIR,
and R-MAON. When predictions come from one scalar field `s(m,l)`, cycle consistency is an algebraic
identity; adding a cycle loss does not add observations or independent constraints.

The preregistered A0/G0/T0 audit ran in the `drug` environment. It inspected the schema of the
affinity-bearing TRAIN registry and materialized only these preregistered safe metadata columns:
`target`, `conn`, `endpoint`, `scaffold`, `assays`, `docs`, `accession`, `hcluster`, and
`dual_cold_split`. It materialized no affinity column or value. Results were:

| Quantity | Label-free result |
| --- | ---: |
| TRAIN rows | 201,827 |
| Targets / accessions / homology components | 559 / 559 / 517 |
| Ligand parents / scaffolds | 121,401 / 48,234 |
| Assay IDs / document IDs | 23,569 / 9,587 |
| Same-target/endpoint groups with at least 16 ligands and an identical aggregated assay-ID bundle | 3,363 |
| Canonical protein-intervention fields available | 0 of 3 |
| Directed ligand-edit field available | 0 of 1 |
| Admissible pair-specific teachers | 0 |

The 3,363 count is an optimistic metadata-density statistic, not verified same-assay/context
rectangles. Dense observational groups did not identify bilateral interventions. The registry lacks
`base_protein`, `construct_sequence`, `directed_substitution`, and `directed_ligand_edit`. PBCNet2.0
has released weights but failed exact-membership, rights, legal-pose, fold-local, and
overlap-exclusion gates. No
single frozen physical teacher survived prior docking/holo/native/physics controls. LEXOR-MC failed
its MC0 extraction contract and emits no quantitative pair state. Generic separable embeddings and
random/constant teachers were retained only as negative controls.

Formal results:

```text
TAU_G0_BILATERAL_INTERVENTION_TOPOLOGY_NOT_IDENTIFIABLE_STOP
TAU_T0_NO_ADMISSIBLE_PAIR_SPECIFIC_TEACHER_STOP
TAU_CONCEPT_TRANSFERABLE__PAIR_SPECIFIC_TEACHER_ABSENT__NO_TRAIN
```

Required future controls are absolute-state versus delta targets, removal and shuffle of intervention
conditioning, wrong post-edit state, target-only plus ligand-only projection, within-family wrong
target/mutation, ligand/edit shuffle, pose/interface shuffle, pretraining-overlap exclusion, matched
random teacher, uniform versus frozen label-free weights, and weight shuffle. Failure cannot be
rescued by more offsets, views, loss weight, rank, capacity, backbone, epochs, seeds, or teacher
concatenation.

Reopen G0 only with exact recurring bilateral interventions across independent base proteins,
families, assay contexts, and provenance lineages. Reopen T0 only with a frozen version, available
weights or deterministic engine, resolved training lineage, frozen rights, legal fold-local inputs,
pretraining-overlap exclusion, pair-specific output, and executable destruction controls. Synthetic
S0 would reuse the R-MAON direct regular-null carrier and earn engineering credit only. Real I0
requires OpenMut I0; predictive M0 additionally requires OMUT-C0.
This order applies only to a non-structural, independently measured teacher. Any structure-, pose-,
docking-, contact-, or complex-conditioned teacher remains deferred until OMUT-M1 passes and may
enter only through OMUT-R0.

### 2.11 Evaluation targets must remain distinct

| Success type | What it establishes | What it does not establish |
| --- | --- | --- |
| Representation/pretext | An encoder captures sequence, structure, or compatibility regularity | Affinity interaction transfer |
| Measurement model | Better assay-offset, censoring, or uncertainty modeling | Unseen-target or unseen-ligand mean/reordering prediction |
| Synthetic estimator | Correct null calibration, power, and recovery | Existence of a real mechanism |
| Train-only oracle | Interaction geometry within the observed panel | Target-cold, ligand-cold, or cross-source transfer |
| Mechanism control | Correct semantics beat matched destruction | A practically material performance effect |
| Predictive gain | One split beats a baseline | Correct mechanism, source independence, or cross-family transfer |
| Calibration/conformal | Controlled risk or coverage | Better point accuracy |

Primary performance is a component-macro paired contrast against ligand-only and an otherwise
identical shared-global null. Mechanism inference separately requires the true coordinate to beat
all matched corruptions. RMSE cannot substitute for reordering; AUCPR or listwise pair count cannot
hide a weak continuous effect; seed averages cannot substitute for family/source replication.

### 2.12 Rights and reader boundaries

The `drug` environment supports CUDA, synthetic estimators, and the test suite. Current blockers are
rights and auditable data access:

- DAVIS-Complete GitHub code has no detected license; only the separately licensed CC0 Dataverse
  data may be considered;
- PLATINUM permits only rights and schema inspection until usage and redistribution terms are
  frozen; no labels or manually reviewed gold subset may be opened;
- BindingDB and ChEMBL may supply registry metadata only after safe projection and must retain
  sequence, construct, endpoint, assay, DOI, and provenance;
- MdrDB may be inspected for rights/schema and source indexing only; its mixed outcomes cannot be
  bulk labels;
- `davis_complete.tab` has not been downloaded or read;
- if a preregistered label-blind projection cannot be enforced, topology is `unknown` until X0/I0;
  opening an affinity-bearing table and promising not to inspect one column is not label-free.

An engineering pass does not unlock a scientific claim.

### 2.13 Heterogeneous open-data pretraining can model observations without identifying interaction

The new open-evidence proposals correctly separate a deployable latent interaction from its
measurement process:

```text
eta(t,l) = F(target, ligand)
y(t,l,a,s,e) ~ p_e(y | g_e(eta), assay=a, source=s)
```

Endpoint, assay, document, site, and source may enter an observation head but not `F`. Exact `Kd`
and exact `Ki` require separate heads; censoring requires interval likelihoods; inactive, binary,
ordinal, and `IC50` evidence retains its own semantics. This is more defensible than pooling every
record into one pseudo-pK label, but it solves only an observation-model problem.

Three identification problems remain:

1. If target family, assay, source, or document are nearly deterministic functions of one another,
   architecture cannot separate biological interaction from the observation process. Independent
   bridge cells and graph closure are required.
2. If most ligands occur against one target, `F(t,l)` can still collapse to ligand potency.
   Same-ligand cross-target comparisons are useful only when endpoint, campaign, protocol,
   inclusion, family, and provenance are matched and tested negatives are retained.
3. A label-bearing pretraining checkpoint inherits every target, homology, ligand, scaffold,
   chemical-neighbor, assay, document, and provenance exposure in its training graph. Strict
   dual-cold evaluation therefore requires fold-specific pretraining after transitive closure; a
   global checkpoint with a new output head is not cold.

MEDIP is a conditional system hypothesis, not an established model innovation. Its observation
heads overlap standard censored/multi-task likelihoods and prior NARD, de-noising, STRATA, and
R-MAON work; its ranking and mixed-difference terms overlap SAFSA, TCOPA, REWIRE, WTPAIR, and other
closed routes. The potentially distinct package is endpoint separation plus assay/source exclusion,
matched same-ligand selectivity, fold-specific closure, and a downstream direct regular-null
operator. Every element must survive its own ablation.

MEDIP-S0 then falsified the current package using generated coordinates and outcomes across six
fixed seeds and five equal-budget variants. The correct model recovered mixed differences
(`0.998710`) and cross-target order (`0.984848`), endpoint merging worsened exact-mean RMSE by
`1.707017`, and the separable null lost `0.397727` ordering accuracy. However, metadata shuffle
reduced mixed-difference correlation by only `0.005862`, below the frozen `0.10` gate, even though
it worsened calibration RMSE by `1.677035` and reduced ordering from `0.984848` to `0.913826`.
Selectivity-label shuffle changed ordering by `-0.000947`, rather than the required loss of at
least `0.08`. Thus separate observation heads are useful for calibration and the interaction
carrier can recover the simulated truth, but this redundant design did not isolate incremental
interaction value from metadata or selectivity. It did not compare decoupling with metadata
injection, its selectivity pairs were latent-effect filtered, and the preregistration omitted some
DGP constants and source hashes. The formal verdict remains
`MEDIP_S0_ENGINEERING_CALIBRATION_STOP`; no loss-weight, rank, capacity, seed, or duration rescue is
permitted. `reports/active/medip_s0_independent_audit_2026-07-28.md` freezes the corrected claim.

The most distinct coordinate hypothesis is a frozen, directional local mutation delta learned from
ProteinGym, MaveDB, or biophysical mutation tasks. Those outcomes usually describe fitness,
expression, stability, or function, not affinity. They may teach a substitution operation but
cannot receive DTA credit until one coordinate at a time beats identity, family, composition,
pooled ESM, BLOSUM permutation, matched random PSD, and wrong-mutation/target controls on a powered
source-resolved substrate.

Binder2030 corrects one historical source statement: DOI
`10.1016/j.slasd.2026.100299` and PMID `41740721` establish that the paper is real and indexed. The
abstract does not establish release of the full matrix, retained non-hits, inclusion probabilities,
shared-ligand depth, assay metadata, independent provenance, or model-training rights. Binder2030
is therefore an `OMUT-D0` candidate only; the older "not indexed" statement is superseded, while
its scientific inadmissibility remains unresolved.

---

## 3. Target-Specific Reordering Exists Locally, but No Transferable Coordinate Is Identified

### 3.1 KirHub mutation double differences

| Quantity | 5%-95% nonsaturation result |
| --- | ---: |
| Eligible point-mutant constructs | 222 |
| Eligible genes | 34 |
| Eligible kinase families | 22 |
| Kinase groups | 7 |
| Descriptive informative ligand pairs | 75,596 |
| Gene-macro rank-reversal rate | `0.1133` |
| Reversal 95% LCB | `0.0865` |
| True-WT minus wrong-WT Spearman | `+0.5010` |
| Pairing-advantage 95% LCB | `+0.4481` |
| Median gene absolute double difference | `15.33` percentage points |

This rejects the universal claim that target-specific reordering never occurs. It remains a
single-source, kinase-only, aggregate 1-uM residual-activity dataset without public raw-replicate
reliability. Its pair count is not independent `n`, and 34 genes at `SD=0.10` give MDE80 about
`0.048`. The decision is `SIGNAL_PRESENT_CONFIRMATORY_GATE_BLOCKED`: ordinal direction check only,
no affinity predictor.

### 3.2 CROSSDOC gives a small cross-table indication

The strict subset required exact single-document, endpoint-separated measurements, at least five
common ligands, and potency residualization. It retained 13 target-document units and 11 homology
components:

- direct KirHub-to-ChEMBL rank correlation `+0.5337 [+0.3101,+0.7336]`;
- group-residual correlation `+0.4946 [+0.3156,+0.6727]`;
- ligand-permutation null `+0.0163 [-0.0136,+0.0426]`;
- observed minus ligand permutation `+0.4783 [+0.3031,+0.6524]`;
- within-document/group target permutation `p=0.001499`.

This supports directional recurrence in a small subset, but it missed the frozen coverage gate of
30 units and 25 components. RECRO later showed why documents are not automatically provenance
families. CROSSDOC is a credible indication, not independent confirmation.

### 3.3 WTPAIR directly supervised mixed differences and still failed

WTPAIR used 25 strict homology-by-chemical-component folds, 308 homology components, 20,000
mixed-difference examples per fold, and a 256-coefficient bilinear ridge:

| Arm | Component-macro Spearman |
| --- | ---: |
| WTPAIR pooled-ESM operator | `0.0323 [0.0112,0.0538]` |
| Ligand-only | `0.0429 [0.0236,0.0622]` |
| KLIFS-group centroid | `0.0911 [0.0703,0.1119]` |
| Matched 256-parameter cellwise bilinear | `0.0419 [0.0227,0.0611]` |
| Within-group protein shuffle | `0.0099 [-0.0102,0.0303]` |
| Matched random protein | `0.0302 [0.0093,0.0512]` |

WTPAIR minus ligand-only was `-0.0107`, minus group centroid `-0.0588`, minus matched bilinear
`-0.0097`, and minus random protein `+0.0021 [-0.0247,+0.0300]`. This rules out three simple rescue
stories on that substrate: missing pairwise supervision, insufficient bilinear capacity, or a useful
continuous direction already present in pooled ESM.

### 3.4 Joint interpretation

> Local target-specific ligand reordering is real, but no source-independent, multi-family,
> semantically specific substitution coordinate has been shown to transfer it into strict dual-cold
> prediction.

---

## 4. Route Evidence Matrix

Each row distinguishes what survived from why extrapolation stopped. A `STOP` applies to the tested
information conditions, not to all biology.

### 4.1 Measurement, pretraining, and early interaction routes

| Route | Evidence retained | Decisive failure | Status or reopening condition |
| --- | --- | --- | --- |
| Measurement-process de-noising | Real assay offsets improved held-out measurement error `24%`; shuffle only `5%` | Downstream cold-target macro-Spearman `-0.051`; 3/3 seeds worse | `DENOISE_FAIL`; measurement-layer use only |
| Pair-compatibility pretraining | Retrieval AUC `0.635`; broken-pair control `0.512` | Downstream `k0` RMSE `+0.0234` versus random initialization; all three seeds worse | `MULTISEED_FAIL`; a new pretext must first be isomorphic to reordering |
| Residue teacher | KL `2.5366 -> 0.0742`; shuffle `0.2609` | RMSE at `k=0/4/16` worsened `+0.0651/+0.0476/+0.0412` | `MECHANISM_PASS_UTILITY_FAIL_STOP`; frozen teacher output must beat ligand-only and target shuffle |
| Atom-residue field | Executable atom graph, ESM residues, and structural priors | B0/B1/B2 RMSE `1.1745/1.2331/1.2461`; Spearman `0.2439 -> 0.2260`; AF prior only 149/977 and gate near initialization | `SINGLE_SEED_FAIL_STOP`; larger field/pocket is not a rescue |
| GO hierarchy | True hierarchy beat shuffled labels | Lost to random initialization over multiple seeds | `BIOLOGICAL_PRIORS_NOT_YET_LOAD_BEARING` |
| NARD robust Student-t/weighting | Seed 1729 briefly improved RMSE and calibration | Three-seed RMSE difference `+0.0035 +/- 0.081`, mixed signs | `MULTISEED_NULL`; no single-seed credit |
| FIELD/REWIRE | Registered interaction field and four-cell estimand; 471 independent transformations | Mechanism gain did not clear MDE; REWIRE MDE about `0.0645` | Data gate closed; reopen only with a new factorial substrate |
| MISO | Endpoint-separated target/compound/edit OOD and structure-fallback contract | Pretraining mechanism did not transfer to the primary estimand | `PRETRAIN_TRANSFER_FAIL` |
| MISO-OR | Valid cross-fitted additive nuisance and orthogonal residual form | `R2-R1=-0.001881 [-0.004785,+0.001019]`; unknown-target and joint OOD negative | Formal STOP; residual-encoder capacity cannot reopen it |
| DICE/MESA-DICE/AXIS | Directed edit by pocket difference, antisymmetric/bilinear carrier, executable edit interface | RAW-CHEM `+0.00119 [-0.00815,+0.01006]`; 518 transformation groups; `BASE_SHORTCUT_FAIL` | Edit is a coordinate arm only; reopen with continuous MMP-dense multi-source data |
| Public true-MMP scale-up | About 2.16M target-conditioned MMPs and 1.45M transformations engineered | Pair scale did not yield target-disjoint interaction gain; continuous high-SNR data and MMP density were inversely related | Engineering scale, not mechanism evidence |

### 4.2 Structure, posterior, Bayesian, and adaptation routes

| Route | Evidence retained | Decisive failure | Status or reopening condition |
| --- | --- | --- | --- |
| BIOFORGE-R2 | One seed gave `B-A=-0.0206` | Three seeds by three folds aggregate `-0.0044 [-0.0143,+0.0057]`; highly correlated heads without residue/atom axis | Archived; new bound-complex, source-resolved substrate required |
| STRATA-M/MX/EV | Executable posterior field and evidence-native likelihood; useful calibration/abstention | Pocket-prior gain `+0.0002 [-0.0001,+0.0004]`, equal to random/shuffled; censor route harmful or neutral | Retain PSF/CRA reliability ideas only |
| PHAROS/CARP-SAGE/MOSAIC/ORBIT | Explored weak supervision, partial pooling, transport, robust objectives | Repeated support/representation controls found no load-bearing target channel; several single-seed false positives | Closed; renaming or expanding the backbone does not reopen |
| SAFE-BMA/MAPLE/HARP/BORA | Coherent posterior, firewall, AHP, and executable registry engineering | Dependency clustering made protocols underpowered; some registries incomplete; priors lost controls | `*_PROTOCOL_UNDERPOWERED`; engineering only |
| DAVIS-HARP/PCMA/SABER/PME | Unified pipeline and lightweight real-data execution | Support router/posterior/meta arms lacked correct-support specificity or confirmation power | Includes `SABER_PHASE_B_FAIL_STOP`; new independent units required |
| HQ-GBMA Stage D/PARC | Synthetic coordinate sensitivity `+0.0431 [0.0185,0.0696]`; pocket slightly beat pooled ESM | Unstable ratio, 77 units, position did not beat composition; shared-global beat ESM but its advantage over PARC was unresolved | `MODEL_SIDE_ROUTE_NOT_IDENTIFIABLE_STOP`; a real 3Di arm needs its own matched control |
| Hierarchical Bayes/Transformer meta | Some RMSE point gains around `+0.02` | Gains disappeared when protein controls were broken; family prior did not carry load | `PROTEIN_CONDITIONED_PRIOR_NOT_LOAD_BEARING` |
| SCGD/QACO | Tested correct-support adaptation on the small Reinecke panel | Wrong-target support was no worse than correct; protein-free no worse than protein | Both FAIL-STOP; a larger adapter or more support is not a rescue |
| RB-DR-QMAPD oracle | Teacher A increased monotonically with target evidence; secondary `Delta_info +0.0333` | Primary `+0.0154 [-0.0155,+0.0464]`, MDE `0.0452`; total `+0.0262` below `0.05` | Oracle gate failed; no student built; needs more independent-component power |

### 4.3 Physical-structure program

Gate-P2A showed conditional chemistry identification on 82 disjoint native holo systems: full win
`0.7033`, LCB `0.6382`, with bounded drops under chemistry, pair, ligand, and receptor shuffles.
This establishes chemistry conditional on a reference pose, not affinity.

Four later falsifications closed the affinity route:

1. Davis AlphaFold-docking within-target Spearman about `0.005`.
2. Experimental holo minus AlphaFold on eight matched targets `-0.019`.
3. PLINDER native-pose upper bound: joint `0.257 < ligand marginal 0.394 < heavy-atom 0.429`.
4. Physics reduced descriptor baseline `0.381 -> 0.125`, gain LCB `-0.395`.

Formal status: **`PHYSICAL-STRUCTURE PROGRAM CLOSED`** for pose/contact/docking rescue on the current
observational graph. Reopening requires a genuinely new strict-dual-cold factorial substrate where a
frozen structure coordinate first beats ligand descriptors, pose/target shuffles, and matched
capacity.

### 4.4 Registries, panels, and source routes

| Route | Evidence retained | Decisive failure | Status or reopening condition |
| --- | --- | --- | --- |
| Original dual-cold few-shot baseline | Target/ligand registry and baseline engineering | Protocol mismatch and repeated cell reuse during registry expansion; early `+0.093` claim withdrawn | Never cite as performance evidence |
| CFRI Gate Z/BM0/BM1/BM2 | Explicit target, ligand, component, and support firewalls | Correct support approximately matched wrong/permuted support; dependency MDE above threshold | FAIL-STOP; firewall specification retained |
| Metz PA/PB/PC | Real interaction structure: PA2 about `0.663`, PA5 `p=0.000488`, held-cell signal reproducible | One paper, kinase-only, many within-target ties; protein-conditioned subspace lost to shared-global | Historical train-only mechanism substrate |
| Local `panel_davis` | Clean endpoint and complete kinase panel | Median query depth 12; adjusted MDE `0.1596`; WT overlap with DAVIS-Complete | Target-conditioned gate unconsumed but closed; F0 must irreversibly preserve or retire it |
| Reinecke 2024 | Open pKd_app panel; 109 targets, 104 components | 80 scored components, median query depth 5, MDE80 `0.0668` | Development-only; SCGD/QACO/SI0/MIF stopped |
| KirHub SPKOP A1 | Pooled ESM beat ligand/shuffle/random; `ESM-ligand +0.0290` | Below `0.030`, lost to group centroid; one source and kinase-only | Weak protein signal, not a transferable coordinate |
| CAPIT/ASPIRE pocket oracle | Aligned pocket `0.4539`, above group `0.4253`, ESM `0.3983`, shuffle `0.3562`, random `0.3670` | Used training-target held-ligand profiles; pocket-group `+0.0286`, short of threshold by `0.0014` | Ligand-warm oracle only; not deployable |
| Taxonomic Resolution TR-0 | Directly tested own-group dependence | Leaving out own-group centroid retained about 77% of group gain | `TR0_PREMISE_FAIL_STOP`; coarse taxonomy is not the mechanism |
| PFSC | Attempted document-isolated, multi-family, single-cold evaluation | Only 19 scorable components after strict document-set isolation; MDE `0.064` | `PFSC0_FAIL_STOP`; public overlap insufficient |
| AMOB | Found apparent ordinal/cross-environment signal | Missing assay/document IDs; superseded by RECRO provenance audit | `OPEN_DATA_INSUFFICIENT_FOR_AMOB` |
| RECRO R0 to L0 | Raw ligand potency family-disjoint correlation remained `+0.700` | Target-specific residual collapsed after provenance isolation; correct-minus-wrong negative | `RECRO_SIGNAL_EXPLAINED_BY_PROVENANCE` |
| RAMCI/SPD | Multi-family roster and standardized assay clues | Median 14 compounds per gene-assay; about 91.3% censored inactive; no continuous dual-cold rectangles | `STANDARDIZED_PANEL_NOT_DUAL_COLD_CAPABLE` |
| OpenBind | Dense local support | One protein only | Local SAR evidence; no target-cold claim |

### 4.5 Coordinate, operator, and evidence routes

| Route | Evidence retained | Decisive failure | Status or reopening condition |
| --- | --- | --- | --- |
| SAFSA family anchors | Plausible family-selectivity pretraining hypothesis | Insufficient supervision structure and independence after firewall | STOP; family taxonomy cannot be relabeled as mechanism |
| MMP-X | Plausible chemical-transformation anchor | Inadequate label/accession firewall and continuous-endpoint topology | STOP; MMP remains a control arm |
| TCOPA | Target-contrast pretext idea | Same EPA source, incompatible endpoint, no independence | STOP |
| Papyrus F0-P | Large aggregated activity source | No document-replicated cells after resolution; conflicting aggregate provenance | `PAPYRUS_F0_RAW_PROVENANCE_INSUFFICIENT_STOP` |
| PANEL-EVIDENCE/Mamba E0 | Local evidence layer and component engineering passed | Supplies neither supervision nor interaction information | Engineering only |
| K-LBP R1/R3 | Deterministic label-free proxy was nontrivial | Factorized null nonregular; R3 `NO_DECISION` | R3 frozen; replace only with a direct regular-null carrier |
| LEXOR/evidence compiler | Source-bound evidence cards are auditable in principle | Live-source transcription/observability failed; old LEXOR-R/P and unrestricted RAG retired | Reopen only with independent, citable, label-memory-free fixed evidence |
| R-MAON G0 | Direct regular-null synthetic calibration passed | No prospective multi-family topology; real I1/affinity untested | Module-only pass; `RMAON_G0_TOPOLOGY_OR_POWER_STOP` |
| LOCK G0 | Alignment, correct sequence, target assignment, and conservation carry information | Fixed BLOSUM lost to identity/permutation/random PSD | `LOCK_G0_REORDERING_NOT_IDENTIFIED_STOP` |
| Recent teacher/paper audit | PBC pairing, FLOWR dense adaptation, and conformer pooling are useful in their source tasks | Membership, privileged poses, private SAR, wrong estimand, or noncompliant topology | `NO_NEW_IDENTIFIABLE_INTERACTION_SOURCE__NO_TRAIN` |
| Tau A0/G0/T0 | Training-only intervention-conditioned delta supervision can increase within-unit supervision density | Current registry has no bilateral intervention fields; PBCNet, physical, and LEXOR teachers all fail admission | Concept pass only; `TAU_CONCEPT_TRANSFERABLE__PAIR_SPECIFIC_TEACHER_ABSENT__NO_TRAIN` |
| Open evidence / MEDIP / local mutation coordinate | Endpoint-separated observation heads, fold-specific closure, and directional local mutation deltas define falsifiable hypotheses | MEDIP-S0 metadata/selectivity destruction failed; same-ligand topology, source independence, task alignment, power, and real semantic transfer remain untested | MEDIP current form stopped; D0 metadata only; `OPEN_EVIDENCE_RECOVERY_ADMISSIBLE__INTERACTION_PRETRAINING_HYPOTHESES_UNTESTED__NO_REAL_LABEL_TRAIN` |
| OpenMut/Delta2Rank | Data recovery is conditionally feasible; canonical WT-mutation by directed ligand change remains untested | D0 incomplete, F0 `NO-GO`, D1 topology/power unknown; quartet/operator naming is not novel | Only active route; affinity training remains prohibited |

---

## 5. Data Substrate Ledger

Counts use their explicitly stated scopes. Full registry, TRAIN, development, mutation subset, and
strict-component subset counts are not interchangeable. The role column describes historical or
conditional scientific use; the authorization column states what can be done now.

The physical local inventory contains 115 files in 10 top-level source directories, totaling
1,470,606,023 bytes (1.37 GiB). These files form a layered evidence store, not one exchangeable
affinity table. Exact paths, formats, footer-level schemas and counts, license status, consumption
state, permitted roles, and prohibited uses are recorded in
`reports/active/local_dataset_inventory_2026-07-28.md`. That inventory opened no affinity row or
value. File presence is not training authorization.

| Substrate | Frozen scale or signal | Historical or conditional role | Current authorization and decisive limit |
| --- | --- | --- | --- |
| ChEMBL-37 dual-cold TRAIN | 201,827 rows; 559 targets; 121,401 ligand parents; 517 target components; 2,197 scaffold-or-document components; largest ligand component 90,288 (`74.37%`), covering 163,117 rows (`80.82%`) and 527 targets | Historical non-affinity TRAIN topology and engineering audit only | No new D1 credit or affinity use; historical PA2 `0.356 < 0.5`, missing provenance-family, and giant component prevent a compliant split |
| Metz `CHEMBL1201862` pKi | TRAIN 12,574 cells, 112 targets, 101 components, 619 ligands; development median query depth 43; PA2 `0.66292` | Historical train-only interaction/estimator development | No new training; kinase-only, one paper, values rounded near 0.1 pK, median within-target ties about 81.7%; development already used |
| Local `panel_davis` pKd | TRAIN 2,069 cells, 116 targets, 102 components; development 1,360 cells and 26 query ligands; median query depth 12 | Single-use target-conditioned confirmation candidate | Closed pending F0; `consumed=false` only for target-conditioned use, while arm-blind power labels were read; retraining SD `0.2314`, raw MDE `0.0688`, adjusted `0.1596`; overlaps DAVIS-Complete WT |
| DAVIS-Complete v3.0 | Dataverse DOI `10.7910/DVN/RTQGP1`, CC0; 444 FASTA headers; 39 candidate Hamming-one records/pairs, including 37 clean canonical-WT substitutions across about 10 named bases, paper convention 11; 72 ligands | Metadata and, only after F0, development diagnostics | Affinity table unread; optimistic projected MDE80 `0.089/0.085` assumes paired SD `0.10`; GitHub code has no license; never independent confirmation |
| PLINDER TRAIN | 5,804 rows; 2,221 clusters/3,828 accessions; 2,717 ligands; 1,081 components; 217 clusters with at least four ligands, median depth 7 | Historical engineering diagnostics | No new label work; too sparse and lacks an independent provenance/profile design |
| Reinecke pKd_app | 826 cells; 109 targets, 104 components, 171 ligands, 121 scaffolds; 80 scored components; median query depth 5 | Historical development-only substrate | No confirmation role; SD `0.2930`, MDE80 `0.0668` |
| KirHub | Mutation DD: 222 constructs/34 genes/22 families/7 groups; LOCK: 353 genes/303 components/92 ligands | Historical ordinal mutation-mechanism audit | No affinity training; aggregate 1-uM endpoint, no raw replicates, one source, kinase-only |
| CROSSDOC/RECRO | CROSSDOC 13 target-document units/11 components, Spearman `0.4946`; RECRO 8,432 comparisons/984 documents to 463 provenance families | Small directional indication and provenance counterexample | Historical evidence only; `79.75%` exact duplicates, `91.6%` comparisons within one provenance family; family-disjoint residual unresolved |
| Papyrus 05.7++ | 707,461 aggregated rows; 147,434 strict resolved; zero document-replicated cells | Historical source audit | No current label use; semicolon provenance cannot be split; 101,520 rows have document-field conflicts |
| BindingDB July 2026 | About 3.24M records; staff-curated CC BY 3.0, ChEMBL-derived CC BY-SA 3.0; Articles zip 18,114,757 bytes containing a 328,109,536-byte TSV; old exact/native subset had only 38 pKi and 6 pKd targets with article blocks of at least 40 ligands | Conditional OpenMut evidence registry | D0 metadata/hash/schema only; affinity rows unopened; requires a safe blind projection, DOI-level collapse, exact construct/endpoint review, and lineage plan |
| ChEMBL-37 variants | Variant accession, mutation, residue, and sequence fields exist upstream | Conditional bounded API or full-database registry | D0 acquisition/schema planning only; full local layer missing and variant-to-canonical component referential integrity requires review |
| PLATINUM | 687,373-byte CSV; over 1,000 mutation-affinity records from over 180 papers and over 200 ligands | Possible future manually reviewed gold source | Rights and schema audit only; no explicit data license is frozen, so no labels or manual gold indexing |
| MdrDB | Mutation/resistance index with academic-use terms | Possible paper index for direct biochemical records | Rights/schema/source inventory only; commercial terms unresolved and cell response, predicted structure, docking, and simulation are excluded |
| KLIFS 2026-07-22 | 555 human kinases, 13,325 structures, 4,213 ligands, and 11,250 valid interaction-fingerprint complexes in the frozen local snapshot | Pocket alignment, kinase taxonomy, and mechanism infrastructure | Never an affinity source; the prior MNI-0 target-conditioning audit failed, and family/pocket labels cannot count as interaction evidence |
| OpenBind EV-A71 2A | One protein with 44 local files; the `raw` directory is empty and the retained objects are a local SAR/structure benchmark | Local structure and engineering evidence | No current label use; one-protein topology cannot establish target-cold or cross-family transfer |
| Novartis SPD 2023 | Systematic multi-family safety panel with retained inactive/censored measurements | Censoring, source-design, and power context | No current label use; about 91.3% censored inactive and insufficient continuous query depth for the primary DTA estimand |
| ToxCast invitrodb v4.3 projection | Lightweight target-contrast projection from a systematic screening release | Metadata/pretext design context | No current label use; it is not an independent continuous-affinity source and cannot confirm overlapping public chemistry |

### 5.1 Provenance normalization and blocking rules

Normalization does not delete raw evidence:

1. Copies across databases that trace to one DOI, PMID, patent, institution/experimental lineage, or
   original experiment collapse into one provenance lineage.
2. Identical measurement cells repeated under multiple document IDs collapse as database copies.
   Raw replicate IDs inform reliability but never increase biological `n`.
3. Ligand salts, aliases, and repeated connectivity normalize to one parent. Scaffolds and
   high-Tanimoto neighbors need not be deleted, but transitive closure puts them in one chemical
   component.
4. Target aliases collapse only when exact construct and sequence agree; genuine isoforms and
   mutations remain distinct.
5. A Papyrus semicolon-combined aggregate row cannot be split into pseudo-observations.
6. Primary inference uses paired complete-case biological units after joint, transitive closure over
   base protein, homology/family, and provenance. Marginal counts remain descriptive.
7. Pair, quartet, conformer, teacher-pair, fold, seed, and technical-replicate expansion never
   increases independent `n`.
8. Operationally independent sites are provenance and reliability strata. They do not multiply the
   number of biological WT-mutant pairs.

### 5.2 Minimum fields for an OpenMut primary record

```text
source/version/license/record_id
DOI/PMID/patent/document/assay/institution/provenance_lineage
endpoint/relation/value/unit/temperature/pH/replicate
WT_accession/WT_construct/WT_full_sequence
mutant_sequence/directed_mutation/mutation_type
ligand_parent/canonical_SMILES/scaffold/chemical_component/edit_class
evidence_location/extraction_method/reviewer
```

A record without original evidence location, exact construct, endpoint, and assay/context cannot
enter a primary rectangle.

---

## 6. Minimum Candidate Estimator, Only If All Data and Semantic Gates Unlock

No estimator is currently qualified or authorized for affinity training. Only after
`D0 -> F0 -> D1 -> X0 -> I0 -> C0` all pass may the following candidate be fitted:

```text
delta_m(l) = y(mut,l) - y(WT,l)
z_m(l) = delta_m(l) - weighted_mean_l(delta_m(l))

s(m,l)
  = s0(l)
  + [u(m) - mean_train(u)]^T
    Theta
    [phi(l) - mean_train(phi)]

C_hat(m;a,b) = s(m,a) - s(m,b)
```

### 6.1 Why the directly centered operator is the allowed carrier

- Full and null use the identical ligand-only/shared-global `s0(l)`.
- Train-only centering prevents target intercept and generic ligand potency from entering the
  interaction channel.
- `Theta` is optimized directly, so `Theta=0` is the unique regular null and exactly recovers `s0`.
- A scalar potential gives ligand-swap antisymmetry and cycle consistency; extra cycle/Hodge loss is
  redundant.
- The model fits each mutation's centered response rather than expanding ligand pairs into
  pseudo-independent rows.

### 6.2 Forbidden parameterizations

```text
Theta = U V^T
Theta = gamma * a * c^T
```

Any direction unidentified under the null is forbidden in the primary test. A rank-1 sensitivity
arm is allowed only when both directions are frozen entirely without affinity labels and only one
signed scalar is fitted.

Protein coordinates cannot be concatenated. Each coordinate enters the same estimator separately
with identical `s0`, weights, ridge, split, and budget so that coordinate quality is not confounded
with capacity.

### 6.3 Frozen coordinate and destruction battery

| Axis | Required coordinates | Required destructive controls |
| --- | --- | --- |
| Protein/mutation | Position only; aligned identity/one-hot; pocket composition; fixed BLOSUM; physicochemical mutation delta; pooled ESM difference | BLOSUM-label permutation; directed-substitution-label permutation; matched random PSD; sequence shuffle; position shuffle; matched wrong mutation/target; within-family wrong mutation/target |
| Ligand/edit | Ligand-only; whole-molecule Morgan/physicochemical difference; MMP edit if topology supports it; fixed ligand basis | Ligand-edit permutation; quartet recombination; scaffold/chemical-neighbor block; random interaction |
| Nuisance/source | Additive null; shared-global `s0`; document/assay/source-only predictor | Provenance permutation; wrong document/source; endpoint split |
| Deferred structure | Conservation; local PLM; CLOCK; static structure; conformer ensemble | Matched capacity; structure shuffle; pose shuffle; family-held-out map; 2D ligand-only |

### 6.4 C0 success criteria

A coordinate passes only if all conditions hold:

1. The paired component-blocked contrast of true semantics against the strongest nested
   nonsemantic null has LCB95 `> 0`.
2. The same contrast exceeds the separately frozen material-effect threshold and empirical MDE at
   80% power.
3. It beats identity, family/taxonomy, composition, pooled ESM, and matched random geometry.
4. Correct target/mutation beats matched-wrong and within-family-wrong assignments.
5. The ligand coordinate beats ligand-only and reasonable whole-molecule controls.
6. Source-only, assay-only, and document-only predictors cannot reproduce the effect.
7. No worst-family or worst-source collapse occurs.
8. Inclusion, censoring, exact-only, and missingness sensitivities preserve the conclusion.

Failure stops that coordinate. Bandwidth, rank, capacity, threshold, folds, seed count, and training
duration cannot be changed as rescue. Registry work may continue only to locate a materially
different substrate.

---

## 7. Ordered OpenMut Queue

| Stage | Current state | Action | Unlock condition |
| --- | --- | --- | --- |
| `OMUT-D0` | **Complete, pass** | Freeze source/version/rights/URL/size/hash/schema and safe-download plan; no values | Satisfied: `manifests/omut_d0_source_registry.v1.json`, reproduced by `research/omut_d0.py --verify` |
| `OMUT-F0` | **Complete, pass** | Freeze Davis cell-overlap policy and DOI-level provenance firewall | Satisfied: `PRESERVE_CONFIRMATION` recorded in `manifests/omut_f0_davis_policy.v1.json`; enforced by `research/omut_f0.py::davis_access_guard` |
| `OMUT-D1` | **Complete, adequate** | Label-free topology, independent-unit graph, effective rank, projected/optimistic MDE | Satisfied but concentrated: 62 `k=4` components / 17 accessions / 11 connected components / effective rank 25 / MDE80 0.034, dominated by HIV Gag-Pol and ABL1/EGFR kinase resistance |
| `OMUT-X0` | Blocked | Build evidence-bound BindingDB/ChEMBL/PLATINUM/supplement registry | D1 passes; exact endpoint/construct evidence available |
| `OMUT-I0` | Blocked | Exact/censored reordering variance, replicate/site reliability, graph, exact-only selection, empirical MDE | X0 passes; real signal exceeds noise with adequate power |
| `OMUT-C0` | Blocked | One coordinate at a time in one fixed estimator | I0 passes; true semantics beat all nested and destructive controls |
| `OMUT-M0` | Blocked | Direct centered operator, one seed | X0, I0, and C0 pass |
| `OMUT-M1` | Blocked | Multi-seed held-provenance validation | Mechanism and effect reproduce |
| `OMUT-R0` | Blocked | Frozen local PLM or structure/CLOCK, separately | M1 passes under the same estimator and budget |
| `OMUT-T0` | Blocked | Multi-substitution/additivity test | Cross-family/source mutation evidence exists |
| `OMUT-T1` | Blocked | Strict family-cold plus scaffold-cold DTA | T0 and coordinate controls pass |
| `CONFIRM` | Blocked | Newly sealed independent measurement source | Never DAVIS-Complete, a training source, or a database duplicate |

### 7.1 D0 blocker resolution, 2026-07-28

D0 closed each source contract with a recorded determination rather than with acquisition. A
recorded absence of licence is a passing D0 record and a blocking downstream fact. Full evidence:
`reports/active/omut_d0_decision.md`, `reports/active/omut_d0.json`.

1. **PLATINUM rights: unresolved, and now the binding constraint.** The data page shows an Open
   Knowledge badge and states no licence identifier, no terms of use, and no redistribution grant.
   Determination `no_explicit_license`, plan `blocked_rights`. Its 51-column schema is nevertheless
   exactly the shape this program needs: `affin.k_wt`, `affin.k_mt`, `affin.delta_k`,
   `affin.fold_change`, `affin.unit`, with `mutation`, `mut.uniprot`, `mut.is_single_point`,
   `mut.wt_pdb`, `mut.mt_pdb`, `mut.doi`, `mut.pmid`. Availability is not a licence.
2. **ChEMBL variant layer: bounded acquisition, resolved.** There is no `/variant_sequence`
   resource; the layer is a nested seven-field schema on `assay`. 20,150 assays carry a
   `variant_sequence`; 119,801 activities carry an `assay_variant_mutation`, which is 0.49% of
   ChEMBL_37's 24,527,044 activities, so no full dump is required. ChEMBL's own warning that
   variant sequences are not referentially linked to component sequences makes document-level
   engineered-versus-disease review an `OMUT-X0` requirement.
3. **BindingDB: schema frozen from the verified local copy.** 202607 Articles archive SHA-256
   verified, one TSV member of 328,109,536 bytes, 640 columns read header-only, carrying
   `Article DOI`, `PMID`, `Institution`, `Authors`, per-chain UniProt identifiers, and per-chain
   target sequences. The blind projector and DOI-level lineage collapse remain X0 work.
4. **MdrDB: academic use granted, scientific exclusion unchanged.** The download page grants a
   non-exclusive, non-transferable academic licence and requires contact for commercial use; the
   GitHub mirror declares none. MdrDB aggregates GDSC, DepMap, AIMMS, KinaseMD, PLATINUM, TKI, and
   RET, so it re-exports PLATINUM and is not an independent provenance lineage. Usable at most as a
   paper index; cell response, predicted structure, docking, and simulated affinity stay excluded.
5. **Binder2030: unresolved.** Europe PMC reports the article as not open access with no indexed
   supplementary files; the version of record is CC BY-NC-ND, an article licence and not a data or
   model-training licence. Eight fields remain blocking, including whether any WT-to-single-
   substitution construct pair exists in it at all.
6. **ProteinGym (PG_v1.3, MIT) and MaveDB (API 2026.2.7, AGPL-3.0): unresolved.** Both licences
   cover software, not the aggregated assay data. Fitness, expression, stability, function, and
   MAVE scores are not affinity labels under any circumstance.

### 7.2 Davis F0 is an irreversible choice

DAVIS-Complete contains the original Davis WT panel plus modified entries and can never be an
independent confirmation source. Before any `davis_complete.tab` affinity value is read, freeze
cell-level overlap keys and choose exactly one policy:

1. **Preserve local `panel_davis` for target-conditioned confirmation.** Exclude every overlapping
   DAVIS-Complete WT value, which generally prevents the desired WT-mutant contrast.
2. **Use overlapping DAVIS-Complete WT values.** Retire local `panel_davis` from confirmation before
   the read and permanently mark all Davis-derived records development-only.

The manifest's `consumed=false` means the target-conditioned gate is unused. It does not erase the
historical arm-blind power audit. Silent reclassification is forbidden, and any later confirmation
claim needs a new independent source.

### 7.3 D1 must be genuinely label-free

D1 may use only:

1. separately published entity, sequence, construct, ligand, endpoint, assay, document, and
   observation-presence metadata; or
2. a preregistered blind projector that irreversibly drops numeric affinity, relation, and censor
   fields before materialization and logs source hash, projected schema, dropped columns, and row
   count.

If neither is possible, report topology as `unknown` and defer the source to X0/I0. Reading an
affinity-bearing table and ignoring its value column is not label-free.

D1 reports only:

- counts of mutations with at least 4, 8, and 16 shared ligands;
- candidate mutations, canonical base proteins, homology components, broad families, documents,
  institutions, and provenance lineages;
- directed substitution classes recurring across independent proteins;
- mutation-position by family by document contingency;
- recurring ligand-edit classes, scaffold diversity, and chemical components;
- component-graph connectivity and effective rank;
- component-level optimistic/projected MDE.

D1 must not inspect relation/censor fields, replicate covariance, reliability, exact-only selection,
or empirical MDE. Those belong to I0 after X0. If the label-free graph cannot power the
preregistered material effect after clustering, stop before affinity access or modeling.

### 7.4 Tau is a conditional queue, not the active program

| Stage | Executed state | Meaning |
| --- | --- | --- |
| `TAU-A0` | Complete, conceptual pass | Preserve only training-time intervention-conditioned delta supervision; no independent-`n` claim |
| `TAU-G0` | Complete, stop | Current TRAIN metadata cannot identify canonical protein and ligand interventions |
| `TAU-T0` | Complete, stop | No frozen pair-specific teacher passes lineage, rights, input, overlap, and destruction gates |
| `TAU-S0` | Blocked | Synthetic calibration only after G0 and T0 reopen; must reuse the direct regular-null operator |
| `TAU-I0` | Blocked | Non-structural measured-teacher mixed-difference association only after OpenMut I0; structure/pose waits for OMUT-R0 |
| `TAU-M0` | Blocked | Identical inference model with or without auxiliary loss only after G0/T0/S0/I0 and OMUT-C0; structure/pose additionally requires OMUT-M1/R0 |

The stopped tau route does not change the priority `OMUT-D0 -> OMUT-F0 -> OMUT-D1`. A newly
downloadable architecture, more teacher pairs, or additional generated constraints is not a
reopening event.

### 7.5 Open-evidence pretraining is a conditional queue

| Stage | Executed state | Meaning |
| --- | --- | --- |
| `OE-D0` | Active only as part of `OMUT-D0` | Freeze Binder2030, ProteinGym/MaveDB, assay-corpus, rights, version, schema, and safe-reader metadata; no outcomes |
| `MEDIP-S0` | Complete, stop | Recovery, endpoint separation, null, architecture, numerics, and reproducibility passed; metadata and selectivity destruction failed; no rescue or biological credit |
| `OE-G0` | Blocked by `OMUT-D1` | Count same-ligand target contrasts, within-family/campaign matches, recurring substitutions/edits, independent sources, effective rank, and projected MDE |
| `OE-X0/I0` | Blocked | Establish evidence-bound outcomes, reliability, selection, real reordering variation, and empirical power |
| `MUTCOORD-C0` | Blocked | Test one frozen local mutation coordinate at a time against identity, family, composition, pooled ESM, permutation, random PSD, and wrong-mutation/target controls |
| `MEDIP-C0` | Current S0 form stopped | No hyperparameter rescue; a future arm requires a scientifically different preregistered identifiability question plus X0/I0/C0 admission |
| `MEDIP-M0/M1` | Blocked | Small direct regular-null downstream carrier, then multi-seed held-family/provenance validation |

AssayMatch is a train-data selection arm, not proof that assays are identical. It must beat a
simple endpoint/protocol filter and an assay-text shuffle inside inner TRAIN folds. Same-ligand
selectivity is admissible only with retained tested negatives and matched endpoint, campaign,
protocol, family or close homology, inclusion, and provenance.

---

## 8. Prospective Measurement and Real Confirmation

### 8.1 A0 reliability and variance pilot

The minimum prospective design is:

- 12 constructs representing exactly 6 WT-single-mutant base-protein pairs;
- one pair from each of 6 broad protein families;
- at least 16 shared, scaffold-diverse ligands;
- 2 operationally independent sites/provenance lineages;
- separate reagent lots, operators, instruments, raw-data systems, and analysis lineages by site;
- one preregistered Ki or Kd endpoint;
- complete randomized `12 x 16 x 2 = 384` cell inclusion before technical replication;
- retention of inactive, censored, failed, BQL/AQL, and out-of-quantification results;
- ligand roster, parent, scaffold, chemical-neighbor, and plate layout frozen before outcomes;
- blinded processing; panel inclusion probability 1, or frozen nonzero probabilities with IPW if
  incompleteness is unavoidable.

The two sites are operationally independent provenance and reliability strata. They do not turn six
biological WT-mutant pairs into 12 targets. At paired `SD=0.10`, `n=6` gives
`MDE80=0.1143748`, so A0 estimates:

- cross-site reliability;
- mutation-specific reordering variance;
- assay noise and censoring;
- empirical MDE;
- the evidence-based scale required for a later mechanism or predictive study.

A0 is not strict dual-cold performance confirmation.

### 8.2 Conditional later scale

Expansion is permitted only if A0 empirical variance and completeness support it:

- about 88 independent joint-blocked units is only an optimistic mechanism floor for a `0.03`
  effect under `SD=0.10`;
- the historical 70-155 component operator pilot was exploratory R-MAON planning and is authorized
  only after X0, I0, and C0; it cannot establish prediction;
- strict T1 prediction planning is about 423 independent multi-family components;
- strict T1 also requires at least 40 randomized, scaffold-diverse query ligands per target;
- interaction amplitude/PA2 must exceed the frozen floor;
- at least two independent provenance lineages are required;
- final confirmation requires a new sealed source.

After A0 establishes reproducible cross-site reordering and supplies a frozen covariance and
missingness envelope, a D-optimal development tranche may choose complete mini-blocks, for example
two constructs by four ligands. The objective is information about candidate coordinates and
`Theta`, not selection of strong binders. Every candidate block keeps a known nonzero inclusion
probability, adaptive blocks remain development-only, and a randomized independent confirmation
lineage is preserved. Compare D-optimal with uniform random blocks on minimum information
eigenvalue, parameter MDE, family balance, and worst-source balance. Stop if the information matrix
remains rank-deficient, selection collapses to one family/scaffold, or random blocks are equivalent
within the frozen margin.

Final confirmation cannot be DAVIS-Complete, a mirror of local `panel_davis`, or any republication
of a training source.

---

## 9. Established, Not Established, and Prohibited Claims

### 9.1 Established

1. Assay, document, and measurement processes contain real structure; de-noising and provenance
   audits are useful at the measurement layer.
2. Generic ligand potency is a strong baseline that reproduces across provenance.
3. Target-specific ligand reordering exists in some local protein/mutation contexts.
4. Correct aligned pocket, sequence, and target assignment carry information.
5. Pooled ESM carries coarse protein/family signal on some substrates but does not reliably express
   the required reordering.
6. Fixed LOCK is a valid nontrivial kernel, but its BLOSUM semantics are not identified in the
   current DTA graph.
7. K-LBP's factorized null is nonregular, so R3 cannot decide the mechanism.
8. The R-MAON direct operator is calibrated under a synthetic regular null.
9. Public-source registry recovery is conditionally feasible, but compliant topology, power, and
   rights are not yet established.
10. MEDIP-S0 can recover its synthetic interaction and calibrate separated observation processes,
    but its redundant supervision does not establish incremental interaction value from metadata
    or selectivity; the current package is stopped.

### 9.2 Not established

1. Fixed BLOSUM, CLOCK, local PLM, or a structure-conditioned coordinate transfers across families.
2. Any protein coordinate consistently beats identity, family, composition, and matched random
   geometry.
3. Target-specific reordering reproduces in provenance-independent, multi-family continuous
   affinity.
4. Delta2Rank, OMRO, R-MAON, or any deep model works on compliant real affinity.
5. A structure teacher, known pose, conformer ensemble, or UQ method improves strict dual-cold DTA.
6. Any module produces a strict target-cold plus ligand-cold performance gain.
7. MEDIP, AssayMatch, same-ligand selectivity, ProteinGym/MaveDB mutation deltas, or heterogeneous
   open-data pretraining supplies transferable ligand-reordering information.

### 9.3 Prohibited statements

- "Biological information is useless" or "all protein representations fail."
- "LOCK/CLOCK succeeds in dual-cold DTA."
- "CLOCK failed."
- "K-LBP proves the effect is zero."
- "R-MAON is validated on real affinity."
- "KirHub's 75,596 pairs are 75,596 independent samples."
- "CROSSDOC/RECRO independently confirms across laboratories."
- "A public database is downloadable, therefore source-independent and adequately powered."
- "Binder2030 is indexed, therefore its full data are open and training-ready."
- "MEDIP recovered a synthetic interaction, therefore heterogeneous pretraining is validated."
- "8.6M teacher pairs are 8.6M independent labels."
- "A larger Transformer, Mamba, GP, or ensemble can bypass the data gate."

---

## 10. Closed Routes and Reopening Conditions

| Route | Why closed now | Only valid reopening information |
| --- | --- | --- |
| Pooled ESM reparameterization | Substrate-specific failure: in KirHub it beat shuffle/random but not the group centroid; in Metz/PARC it lost shared-global and matched derangement | New multi-family, multi-source factorial panel where it beats all controls in the same estimator |
| Fixed LOCK/BLOSUM | Lost aligned identity, BLOSUM permutation, and matched random PSD | New source-resolved mutation rectangles where true semantics pass component/family inference |
| True CLOCK | Not tested and therefore not closed as a biological hypothesis | Externally frozen structure map, leave-family-out mapping, structure shuffle, and enough mutation landscapes |
| Pose/docking/native-complex rescue | Four affinity framings failed to carry load | Materially new structure information and strict dual-cold factorial evidence |
| Unrestricted support/meta/posterior | Correct support lacked specificity or power | Correct support first beats wrong, cross-target, and permuted support |
| Factorized `gamma*a*c^T` | Direction unidentified under the null | Direct regular-null `Theta`, or directions frozen entirely without affinity labels |
| Larger rank/backbone/epochs/seeds | Adds no information or independent units | Not a reopening condition; requires new data or a newly validated coordinate |
| MMP/edit branding | Historical edit features lost whole-molecule controls and target-cold tests | Independent continuous rectangles with recurring directed ligand edits by mutations |
| Provenance-blind aggregation | Document/database duplication manufactured signal | DOI, patent, institution, and experimental-lineage transitive closure |
| LEXOR-R/P or unrestricted RAG/LLM prediction | Source observability and memory contamination cannot be controlled | Frozen, citable, label-memory-free evidence channel with destructive controls |
| External teacher with privileged pose/membership | Inputs unavailable at deployment or training membership unauditable | Exact downloadable manifest, legal label-blind inputs, and overlap-excluded evaluation |
| Tau-style Pair-JEPA/BCEL | No canonical bilateral intervention graph and no admissible pair-specific teacher; algebra overlaps existing routes | Exact recurring WT-mutation and ligand-edit interventions plus a frozen fold-local teacher that beats separable/wrong/random controls |
| MEDIP-S0 current package | Metadata and selectivity corruptions did not remove interaction recovery despite successful calibration | Not a weight/rank/capacity rescue; only a scientifically different preregistered identifiability question after a powered X0/I0 substrate |

---

## 11. Compact Decision Timeline

| Date/stage | Immutable decision |
| --- | --- |
| 2026-07-15 to 2026-07-20, FORT/FORGE | Established cold-target task, sealed-test rules, and base ligand/protein pipeline; multiple biological modules showed successful pretexts without downstream transfer |
| 2026-07-20 to 2026-07-21, measurement/pretraining | De-noising, pair compatibility, NARD, residue/atom field, and GO hierarchy completed matched controls; no multi-seed cold-target gain |
| 2026-07-16 to 2026-07-21, FIELD/REWIRE/MISO/DICE | Exhausted the main mixed-difference, orthogonal-residual, directed-edit, and interaction-field algebra; data and transfer gates failed |
| 2026-07-21 to 2026-07-23, PHAROS/STRATA/HARP/BORA/DAVIS | Separated posterior, Bayesian, few-shot support, and reliability claims; engineering/calibration survived, accuracy or protocol power did not |
| 2026-07-24, FIRE/BridgeFIRE | Conditional pose chemistry was identifiable; four affinity falsifications closed the physical-structure program |
| 2026-07-25, dual-cold registry audit | Withdrew the early `+0.093` few-shot claim after registry/evaluation mismatch; independent units and MDE became mandatory |
| 2026-07-25 to 2026-07-26, dense kinase panels | Metz established train-only interaction; Davis/Reinecke exposed inadequate query depth and cross-arm variance |
| 2026-07-26, KirHub/OpenBind | Mutation DD established local reordering; SPKOP pooled ESM lost the group control; OpenBind remained single-protein |
| 2026-07-26, WTPAIR/CAPIT/CROSSDOC | Pooled-ESM mixed-difference predictor failed; aligned pocket oracle nearly passed; CROSSDOC was strong but only 11 components |
| 2026-07-26, TR/PFSC/AMOB/RECRO | Taxonomy premise failed; public overlap was sparse; RECRO attributed nominal cross-document signal to provenance duplication |
| 2026-07-27, anchor/PARC/model audit | SAFSA/MMP-X/TCOPA/Papyrus stopped; PARC exposed ratio and position-versus-composition problems; model-side rescue closed |
| 2026-07-27 to 2026-07-28, K-LBP | R1 label-free proxy passed; R3 returned `NO_DECISION` because the factorized null was nonregular |
| 2026-07-28, R-MAON | Direct regular-null synthetic calibration passed; absent prospective real topology stopped the program |
| 2026-07-28, LOCK/CLOCK | Fixed LOCK identified alignment/target information but not BLOSUM semantics; final category 3 |
| 2026-07-28, recent-paper audit | PBC/FLOWR/conformer/structure/UQ supplied no compliant new interaction source |
| 2026-07-28, tau transfer A0/G0/T0 | Conceptual supervision-density transfer passed; label-free topology and all pair-specific teacher candidates stopped; no synthetic or affinity stage authorized |
| 2026-07-28, OpenMut cross-review | Public data recovery conditionally feasible; D0 active, F0 `NO-GO`, all later stages blocked, and no affinity training authorized |
| 2026-07-28, open-evidence cross-review | MEDIP, local mutation coordinates, AssayMatch, Binder2030, and D-optimal blocks registered only as gated hypotheses; no real-outcome training |
| 2026-07-28, MEDIP-S0 | Synthetic recovery and endpoint calibration passed, but metadata/selectivity destruction failed; current package stopped without rescue |
| 2026-07-28, `OMUT-D1` | Full BindingDB stream (93,712 rows) plus a bounded 5,000-row ChEMBL variant sample. 62 candidate WT-vs-mutant components at k=4 shared ligands, 17 accessions, 98 shared-ligand-scoped documents, 11 connected components, effective rank 25, projected MDE80 0.034 (0.03 threshold). Cleared the informal 25-component adequacy bar but concentrated in HIV Gag-Pol (P04585) and ABL1/EGFR kinase resistance, same pattern as KirHub and the promiscuous-ligand block; mutation tokens are unverified regex candidates. Unlocks `OMUT-X0` |
| 2026-07-28, `OMUT-X0` | Full no-outcome reconstruction stopped the public-data route. BindingDB produced 37 exact Ki/Kd k>=4 construct components but no source-native assay locator. The 119,801-row ChEMBL census produced 225 k>=4 WT/mutant candidates and recovered all 2,622 requested assays; 227/241 mutant sides were sequence-exact but 0/241 WT sides exposed an exact assay-level construct. Primary topology: 0 components. Verdict `OMUT_X0_EVIDENCE_INADEQUATE_STOP`; I0/C0/M0 and real-outcome training remain blocked |
| 2026-07-28, `OMUT-X1` | Source-native ChEMBL descriptions recovered four strict non-BRAF components: A3EZI9 R155K Ki (6 ligands), JAK2 V617F Kd (4), EGFR L858R Kd (13), and LRRK2 G2019S Ki (13). All execution/firewall gates passed, but only three broad families were represented against the frozen 25-component/six-family requirement. Verdict `OMUT_X1_DESCRIPTION_REGISTRY_INADEQUATE_STOP`; 25 other candidates had only 1-3 exact ligands and motivate a label-free source-accessibility audit |
| 2026-07-28, `OMUT-X2` | Eighty-seven documents carried near-exact pairs; Europe PMC matched 83, exposed nine OA full-text records, missed two, and dispositioned two patents. The optimistic accessible topology reached 16 components, six families, and largest-accession share 0.375 but missed the 25-component floor. All-doc theoretical upper bound is 110 components / 17 families / 0.145 largest share, motivating licensed OA version discovery rather than threshold relaxation |
| 2026-07-28, `OMUT-X3` | OpenAlex formal execution stopped on HTTP 429 with anonymous daily budget zero. No result file or registered X3 verdict was produced; no body or outcome was read |
| 2026-07-28, `OMUT-X3C` | Exact-DOI Crossref metadata dispositioned all 87 documents (85 matched, two no DOI). Two version-matched CC-licensed XML links plus nine EPMC records produced an optimistic 32-component / seven-category / 0.50 largest-accession topology and verdict `OMUT_X3C_LICENSED_LINK_RECOVERY_FEASIBLE`. Sixteen components are one-accession/one-document `P61073` candidates, so this passes only the source upper bound and unlocks X4 construct recovery, not training |
| 2026-07-28, `OMUT-X4` | Corrected source transport found nine valid EPMC main texts and two Elsevier metadata-only XML envelopes. Two safe construct fragments were projected, zero candidate/documents met the exact common-construct forms, and the topology reverted to the four X1 components / three categories. Verdict `OMUT_X4_CONSTRUCT_REGISTRY_INADEQUATE_STOP`; official EPMC supplementary archives remain unconsumed and are the only X5 continuation |
| 2026-07-28, `OMUT-X5` | Nine official EPMC supplementary ZIPs yielded nine parsed PDFs while CSV tables and images remained unread. Poppler review accepted `P08581 D1228V` on the shared `His6-cMet(1038-1348)` construct and rejected an EGFR L858R hit as a multi-substitution TMLR construct. Final strict topology: five components / three categories; verdict `OMUT_X5_SUPPLEMENT_CONSTRUCT_REGISTRY_INADEQUATE_STOP` |
| 2026-07-28, `OMUT-X6` | Label-free scan of 145 WT assays used by frozen near pairs found zero labeled catalog/product tokens, zero HTTPS product URLs, and zero actionable supplier candidates. Verdict `OMUT_X6_REAGENT_VERIFICATION_INSUFFICIENT_STOP`; no supplier page or outcome was read |
| 2026-07-28, `OMUT-F0` | Davis role frozen as `PRESERVE_CONFIRMATION` by recorded human decision: `panel_davis` stays sealed single-use confirmation, unconsumed; the DAVIS-Complete WT-mutant contrast is foreclosed; 0 Parquet rows read. The overlap key, DOI collapse rule, and re-export edges became importable, tested objects (`davis_access_guard`, `provenance_unit`). Unlocks `OMUT-D1` on non-Davis sources only |
| 2026-07-28, `OMUT-D0` | Source freeze complete under a no-value firewall: 16 sources, 25 reads, 0 rows materialized, 0 violations. ChEMBL variant acquisition bounded at 20,150 assays / 119,801 activities (0.49% of ChEMBL_37); PLATINUM carries the exact WT/mutant paired schema but no licence, so rights became the binding constraint; Binder2030, ProteinGym, and MaveDB stay unresolved. Unlocks `OMUT-F0` only; category 3 unchanged |

---

## 12. Authoritative Artifacts

### 12.1 Current decision and OpenMut

- `task.md`: authority for the current queue and allowed action.
- `reports/active/omut_d0_preregistration.md`, `omut_d0.json`, and `omut_d0_decision.md`: the
  `OMUT-D0` source freeze, its no-value firewall ledger, and the gates. Anchor
  `substantive_registry_sha256` `1953cff4c1d7301c51d1ef934c0c5c913f7c154022ad515c53e416bbce8f82f9`;
  `registry_sha256` covers live transport sizes and is deliberately not a stable anchor.
- `manifests/omut_d0_source_registry.v1.json`: the frozen 16-source registry alone, reproducible
  with `research/omut_d0.py --verify`.
- `reports/active/omut_f0_preregistration.md`, `omut_f0.json`, and `omut_f0_decision.md`: the
  recorded one-way Davis role decision, the sealed-panel state check, the frozen overlap key, and
  the DOI-level collapse rule.
- `manifests/omut_f0_davis_policy.v1.json`: the immutable policy record. Reversal requires a new
  recorded human decision, not a code change; `research/omut_f0.py::davis_access_guard` and
  `provenance_unit` are the enforceable objects later stages must import.
- `reports/active/omut_d1_preregistration.md`, `omut_d1.json`, and `omut_d1_decision.md`: the
  label-free topology measurement on BindingDB (full stream) and a bounded ChEMBL sample; the
  62-component `k=4` result, its HIV/kinase concentration, and the unresolved simplifications
  (accession-only grouping, unverified mutation tokens) that `OMUT-X0` must close.
- `reports/active/omut_x0_preregistration.md`, `omut_x0.json`, and `omut_x0_decision.md`: the
  full evidence-bound reconstruction, explicit disposition of every D1 candidate and additional
  ChEMBL accession, complete WT/assay/document queries, zero primary components, and the stopped
  public-data route.
- `reports/active/omut_x1_preregistration.md`, `omut_x1.json`, and `omut_x1_decision.md`: exact
  source-native description pairing, four recovered non-discovery components, and the remaining
  source-evidence coverage gap before any outcome stage.
- `reports/active/omut_x2_preregistration.md`, `omut_x2.json`, and `omut_x2_decision.md`: the
  near-exact document topology, exact Europe PMC availability dispositions, nine OA full-text
  sources, and the 16-component optimistic recovery bound.
- `reports/active/openmut_delta2rank_feasibility_crossreview_2026-07-28.md`: OpenMut/Delta2Rank
  mathematical, source, power, and adversarial cross-review. Superseded on stage status by the
  `OMUT-D0` decision.
- `reports/active/open_evidence_pretraining_crossreview_2026-07-28.md`: MEDIP, local mutation
  coordinate, AssayMatch, Binder2030, D-optimal, fold-closure, and dilemma-by-dilemma cross-review.
- `reports/active/local_dataset_inventory_2026-07-28.md`: directory-complete local source, schema,
  license, role, consumption, and prohibited-use inventory.
- `reports/active/medip_s0_preregistration.md`, `medip_s0.json`, and
  `medip_s0_decision.md`: synthetic-only observation/pretraining module falsification and stop.
- `reports/active/medip_s0_independent_audit_2026-07-28.md`: independent reproduction, claim
  correction, preregistration limitations, and valid reopening conditions.
- `reports/active/recent_interaction_paper_audit_2026-07-28.md`: paper-by-paper teacher,
  conformer, adaptation, UQ, and structure audit.
- `reports/active/recent_interaction_paper_gate_2026-07-28.json`: machine-readable label-free
  paper gate.
- `reports/active/tau_transfer_preregistration_2026-07-28.md`: tau A0/G0/T0 definitions and
  frozen stop rules.
- `reports/active/tau_transfer_decision_2026-07-28.md`: cross-reviewed transfer decision.
- `reports/active/tau_transfer_crossreview_2026-07-28.md`: independent paper/teacher/adversarial
  cross-review; its grouped T0/T1/T2 names are mapped to the canonical live queue.
- `reports/active/tau_feasibility_a0_t0_2026-07-28.json`: executed label-free topology and teacher
  admission results.
- `manifests/tau_teacher_admissibility.v1.json`: one-at-a-time teacher contracts and negative
  controls.

### 12.1a DCST R14-R17 decision packets

- `reports/active/dcst_r14_r17_complete_record_2026-07-29.md`: unified summary of
  R14 absolute transport, SISMT/R15, DTIOD/R16, and KLIFS/R17, including frozen
  estimands, thresholds, controls, numerical results, protected-asset boundary,
  and reopening conditions.
- `reports/active/dcst_r14_transport_support_preregistration_2026-07-28.md`,
  `dcst_r14_transport_support_seed1729.json`, and
  `dcst_r14_transport_support_decision_2026-07-29.md`.
- `reports/active/dcst_r15_sismt_preregistration_2026-07-29.md`,
  `dcst_r15_sismt_seed1729.json`, and `dcst_r15_sismt_decision_2026-07-29.md`.
- `reports/active/dcst_r16_dtiod_preregistration_2026-07-29.md`,
  `dcst_r16_dtiod_t1_seed1729.json`,
  `dcst_r16_dtiod_t1_target_blocks_seed1729.csv`, and
  `dcst_r16_dtiod_t1_decision_2026-07-29.md`.
- `reports/active/dcst_r17_klifs_bridge_preregistration_2026-07-29.md`,
  `dcst_r17_klifs_bridge_seed1729.json`,
  `dcst_r17_target_support_correction_2026-07-29.md`, and
  `dcst_r17_klifs_bridge_decision_2026-07-29.md`.

### 12.2 Protein coordinates and estimators

- `reports/active/lock_clock_g0_decision.md`
- `reports/active/lock_clock_g0_final_audit_2026-07-28.md`
- `reports/active/lock_clock_g0_label_free.json`
- `reports/active/lock_clock_g0.json`
- `reports/active/protein_conditioned_signal_audit_2026-07-27.md`
- `reports/active/parc_m0_decision.md`
- `reports/active/klbp_r3_decision.md`
- `reports/active/post_r3_multiagent_deep_route_selection_2026-07-28.md`
- `reports/active/rmaon_g0_decision.md`
- `reports/active/rmaon_g0_final_audit_2026-07-28.md`

### 12.3 Reordering, provenance, and power

- `reports/active/kirhub_dd_decision.md`
- `reports/active/kirhub_wtpair_c1_decision.md`
- `reports/active/kirhub_pocket_oracle_c2_decision.md`
- `reports/active/crossdoc_reordering_c3_decision.md`
- `reports/active/recro_l0_decision.md`
- `reports/active/pfsc_gate_decision.md`
- `reports/active/rb_dr_qmapd_power_audit.md`
- `reports/active/post_failure_measurement_design_exploration_2026-07-27.md`

### 12.4 Source and firewall records

- `dataset/public/chembl_37/processed/dualcold/manifest.json`
- `dataset/public/chembl_37/processed/panel_metz/manifest.json`
- `dataset/public/chembl_37/processed/panel_davis/manifest.json`
- `dataset/public/plinder_2024_06_v2/processed/dualcold/manifest.json`
- `manifests/open_sources.json`
- `reports/active/open_data_only_amendment.md`
- `reports/active/panel_davis_registration.md`
- `reports/active/rmaon_g0_data_compliance_audit_2026-07-28.md`

Old append-only history referenced intermediate reports that no longer exist. Existing frozen reports
take precedence over this summary. The canonical physical-structure termination numbers and some
archived residue/atom architecture results are retained here and in Git history rather than linked
to missing files.

---

## 13. Firewall and Reproducibility State

### 13.1 Permanent firewall ledger

```text
new_affinity_labels_read_this_round = false
davis_target_conditioned_confirmation_consumed = false
davis_arm_blind_power_labels_historically_read = true
sealed_test_consumed = false
historical_chembl_confirmation_labels_read = true
```

The final line records the five rows displayed during an earlier CROSSDOC schema inspection. They
entered neither training nor statistics, but that partition remains permanently isolated.

### 13.2 Last recorded verification

Using `D:\anaconda\envs\drug\python.exe`:

- current repository-owned suite: `497 passed, 1 pre-existing unregistered
  slow-mark warning in 119.02s`;
- `pytest.ini` fixes collection to `tests/` and `lexor/tests/`, so downloaded third-party code
  under `tmp/` is not mistaken for a FORT test suite;
- `history.md` and `task.md`: zero Han characters, zero non-ASCII characters, and balanced code
  fences;
- all current Tau/OpenMut/open-evidence/MEDIP/LEXOR authoritative paths checked here exist;
- `git diff --check` and trailing-whitespace checks pass for the edited tracked ledgers;
- the X0 frozen-result regression, forbidden-key traversal, gate/verdict
  consistency check, complete 119,801-row census check, and zero-primary-topology
  assertion pass;
- `tools/raw_data_guard.py verify` passes both registered Plinder raw-file
  size and SHA-256 checks;
- tau result/manifest hash binding passes, with zero admitted teacher, no synthetic or affinity
  authorization, and all current-run development/confirmation/sealed flags false; these do not
  overwrite the historical ChEMBL confirmation flag;
- MEDIP-S0 repeated formal suites are byte-equivalent after runtime removal, with canonical
  SHA-256 `7b0e65e54d9fec47352d8dd97b52aea4274ceb78f27c9131fdeae103510a08e5`;
- `davis_complete.tab`: absent locally and unread;
- no affinity predictor started, and no new development, confirmation, or sealed labels were read.

### 13.3 Formal decision

```text
OPENMUT_DATA_RECOVERY_CONDITIONALLY_FEASIBLE__DELTA2RANK_NO_TRAIN
OMUT_X0_EVIDENCE_INADEQUATE_STOP
TAU_CONCEPT_TRANSFERABLE__PAIR_SPECIFIC_TEACHER_ABSENT__NO_TRAIN
OPEN_EVIDENCE_RECOVERY_ADMISSIBLE__INTERACTION_PRETRAINING_HYPOTHESES_UNTESTED__NO_REAL_LABEL_TRAIN
MEDIP_S0_ENGINEERING_CALIBRATION_STOP
```

Final category remains:

> **3: current data do not identify a transferable substitution-geometry mechanism. New prospective
> measurements or a newly recovered, source-resolved, adequately powered public substrate are
> required.**

---

## 14. UBSE-A1-v2 Stage-1 Static Preregistration, 2026-07-30

The Stage-1 extractor/FG/topology contract was frozen without reading a
coordinate body, CCD body, contact/event value, affinity, confirmation, or
sealed outcome.

Authoritative new artifacts:

- `manifests/ubse_a1v2_stage1_extractor_contract.v1.json`, SHA-256
  `e2e66f1143c93ce886af86a976abe4dcb28f5677ab332451998b9d301e638c92`;
- `research/ubse_a1v2_stage1_contract.py`;
- `tests/test_ubse_a1v2_stage1_contract.py`;
- `reports/active/ubse_a1v2_stage1_extractor_topology_preregistration_2026-07-30.md`.

The frozen primary extractor is PLIP 3.0.1 at commit
`2f4911d307490479ac023b22d6faa8f59b577ca8`, source-archive SHA-256
`1d9f3ecdfb02b84f957e83c4cf71cfae5d6a618a2cf77305139cca12b7a48816`.
ProLIF 2.2.1 at commit
`7b7055288a037e60eab38f24a7ea979bb5c91f30`, wheel SHA-256
`f7239aceb11949cb4a4b4d2c32f1f7c9b46792cc65eedc69bc44cb3c492bd806`,
is a non-voting implementation audit. Union, consensus, voting, and backfill
labels are forbidden.

The contract now freezes:

- seven ordered channels and explicit native-record to ligand-role to FG maps;
- role-scoped multi-hot FG ownership and RDKit automorphism collapse;
- formal-charge versus ionizable-heuristic evidence separation;
- model 1 of the deposited asymmetric unit, exact auth locators, whole-instance
  altloc selection, occupancy, HOH, metal, covalent, and missing-atom rules;
- one nullable `0/1/null` row per legal cell per extractor, with sparse
  positives forbidden from defining negatives;
- candidate tetrads from the event-blind mask separately from informative
  post-extraction `1,1,0,0` or `0,0,1,1` patterns;
- A1-R's 153 units, 124 frozen dependency components, assignment SHA-256
  `46bbff8c6be38ffce699b03faabd86871261ad660f2d667f987ea6680a14e83c`,
  and worst-case Kish effective size `98.7721519`;
- fixed assignment under attrition with surviving sizes/Kish recomputed;
- hash-first, no-redirect, write-once coordinate and wwPDB CCD semantics;
- all 15 execution authorizations as false.

The directional-power correction is binding. For `0.60` versus `0.55`,
one-sided Bonferroni `0.05/7`, power is `0.054711` at `n=98` and `0.418954`
at `n=512`; 1,065 independent Bernoulli units are required for 80% power.
The `200 targets` and `4 x 50` requirements are coverage only. Directional
accuracy is auxiliary. The primary component-held-out cycle
log-score/deviance procedure is frozen, but its numeric material threshold
and power curve remain pending A1-S-TRAIN-only calibration before development
or A1-C release.

All 21 declared Python distributions were installed from local hash-frozen
archives in `D:\anaconda\envs\drug`. Version, local artifact SHA-256, PEP 610
installed-archive identity, critical source hashes, imports, safe inputs, and
SMARTS compile checks pass. The focused suite reports `76 passed`; the full
repository suite reports `808 passed` in 141.36 seconds with only the
pre-existing slow-marker warning and an MDAnalysis deprecation warning;
`pip check` reports no broken requirements. CUDA remains available as Torch
`2.6.0+cu124` on one NVIDIA GeForce RTX 4060 Laptop GPU, but Stage 1 is CPU
work.

An independent validator audit found one pre-authorization P1 issue:
`checkerboard_direction` did not require the event-blind evaluability mask.
It was corrected before this record; non-evaluable cells and null values now
raise, and dedicated tests pass. Installed archive provenance is bound, and
every hashed installed file is now checked against its distribution `RECORD`
before the runtime can be ready.

The supplied modeling suggestions were narrowed as follows:

- accept no-P0A fixed-margin tetrads as the primary Stage-2 coupling route;
- defer pose-conditioned geometry to a separately gated Stage-3 incremental
  arm with license, membership, deployment, interaction-recovery, and
  destruction checks;
- retain hierarchy only as diagnostics and retain PLIP as the sole teacher;
- forbid checkerboard-driven reselection of the already frozen roles;
- use the name Fixed-Margin Tetrad Interaction Distillation, not
  counterfactual distillation;
- keep Stage-4A as `B0 + theta^T z_int` and the flexible kernel deferred to
  conditional Stage-4B.

The formal state remains:

```text
NO_DECISION_A1V2_STAGE1_EXTRACTOR_TOPOLOGY_PENDING
coordinate_get = false
coordinate_parse = false
ccd_snapshot_get = false
event_extraction = false
affinity_read = false
student_training = false
confirmation_read = false
sealed_access = false
```

Without additional user authorization, only read-only review and preparation
of a separate coordinate-acquisition authorization packet are allowed.

---

## 15. UBSE-A1-v2 Task-Aligned Action Review, 2026-07-30

The latest model-strengthening recommendations were reconciled against
`task.md` without reopening the frozen Stage-1 contract. The binding
interpretation separates three claims:

1. deposited typed events may add real ligand-conditioned 3D information;
2. fixed-margin tetrads may make that information identifiable beyond row and
   column marginals;
3. a separately gated pose ensemble may add deployment-available
   pair-conditioned geometry.

The accepted route remains **Fixed-Margin Tetrad Interaction Distillation**.
No-P0A is primary, P0A is incremental only, PLIP is the sole teacher, ProLIF
is non-voting, and Level-0/Level-1 aggregates remain diagnostics. Tetrads are
Stage-2 coupling evidence and do not replace the 224-dimensional Stage-4A
readout `theta`. Result-driven A1 role reselection, extractor consensus/union,
the word counterfactual, a single docking pose as truth, and model-capacity
rescue remain rejected.

The ordered recommendation is recorded in
`reports/active/ubse_a1v2_task_aligned_feasible_actions_2026-07-30.md`.
The only current state-changing preparation is an explicit hash-first
coordinate/CCD authorization packet with status `PREPARED_NOT_AUTHORIZED`.
No downloader, coordinate GET, CCD GET, parse, event extraction, model
training, affinity read, confirmation, or sealed access was authorized or
executed.

A runtime audit found that NumPy's installed `RECORD` contains a hashed
`__pycache__` entry that changes when imported. The validator was narrowed to
except only a standard current-interpreter cache whose same-name `.py` source
is `RECORD`-listed and hash-valid. Source-less caches, non-current-tag cache
files, legacy `.pyc`/`.pyo`, ordinary source, and all other hashed payloads
remain checked. The one explicit runtime exception is
`numpy/distutils/__pycache__/conv_template.cpython-311.pyc`, bound to the
hash-valid companion `numpy/distutils/conv_template.py`.

Post-correction verification in `D:\anaconda\envs\drug`:

- frozen contract SHA-256:
  `e2e66f1143c93ce886af86a976abe4dcb28f5677ab332451998b9d301e638c92`;
- `contract_valid=true`, `runtime_ready=true`, and zero runtime failures;
- all 15 execution authorizations remain false;
- focused suite: `77 passed, 1 warning`;
- repository suite: `809 passed, 2 warnings in 126.61 seconds`;
- `pip check`: no broken requirements.

The formal scientific decision remains:

```text
NO_DECISION_A1V2_STAGE1_EXTRACTOR_TOPOLOGY_PENDING
```

---

## 16. UBSE-A1-v2 Mask Closure and Prepared Acquisition Packet, 2026-07-30

A final static audit identified two ambiguities that could not be carried
into A1-S or Stage 3:

1. the deposited-complex 6.5 A mask is valid for A1-R reliability but would
   reveal the true ligand neighborhood if passed to a deployable student;
2. "log-score or deviance" and lexical-basis averaging did not uniquely
   specify the nonlinear primary statistic.

The binding correction is
`reports/active/ubse_a1v2_mask_statistic_closure_amendment_2026-07-30.md`,
SHA-256
`42e0684792595067d703f9b01ef347246ec153caed134ddff32f6d54293b7b76`.
It freezes:

- `M_rel` for deposited A1-R reliability only;
- an all-sequence, 2D-FG, event-role-compatible `M_deploy` that cannot use
  PDB, true pocket, pose, P0A, event, or affinity information;
- complete A1-S/A1-C `M_obs` enumeration without a 6.5 A prefilter and with
  strict `0/1/null` semantics;
- conditional tetrad probability `sigmoid(y * Omega)` and primary deviance
  gain against the exact rank-one null `Omega=0`;
- all elementary informative tetrads for the target-normalized likelihood,
  while the lexical independent basis remains a rank certificate only;
- bitwise mask, margin, dustbin, tetrad-set, and weight preservation for
  destruction controls.

Synthetic tests establish zero deviance at the additive/rank-one null,
positive gain only for the matching checkerboard direction, and rejection of
non-evaluable, non-finite, or invalid-direction inputs.

The machine-readable A1-R/CCD request is
`reports/active/ubse_a1v2_coordinate_acquisition_authorization.v1.json`,
SHA-256
`0b650526ba64160355e1da92d18219239bffbc7f44800b804f7744c26d37e7a2`.
It binds the parent contract, closure amendment, A1-R locator manifest,
selected URLs, and strict HEAD records. Its scope is 421 unique A1-R
coordinate bodies, 459 coordinate instances, and one CCD snapshot. The
packet status is `PREPARED_NOT_AUTHORIZED`; no runner exists, network GET is
forbidden, and all 15 effective authorizations are false.

The installed-file cache exception was further narrowed to an exact
current-interpreter `__pycache__` filename grammar. Substring cache-tag
matches, wrong tags, source-less entries, non-cache `.pyc`, and `.pyo` remain
failures.

Final verification in `D:\anaconda\envs\drug`:

- contract CLI: `contract_valid=true`, `runtime_ready=true`, zero failures;
- acquisition packet:
  `PREPARED_NOT_AUTHORIZED`, exact SHA-256 as above;
- runtime exceptions: exactly one bound NumPy cache with a hash-valid source;
- all 15 execution authorizations: false;
- focused suite: `81 passed, 1 warning`;
- repository suite: `813 passed, 2 warnings in 225.24 seconds`;
- `pip check`: no broken requirements.

No coordinate body, CCD body, event, affinity, confirmation, or sealed
outcome was read. No downloader or acquisition runner was created. The
scientific decision remains
`NO_DECISION_A1V2_STAGE1_EXTRACTOR_TOPOLOGY_PENDING`.

---

## 17. UBSE-A1-v2 Falsification Round: Phase A, 2026-07-30

The user instruction at SHA-256
`e06b32c7f6f8d99252ffb160e57bdbc92458e81dd089d4aee4a543ab92083983`
switched A1-v2 from static review to the ordered `C1 -> C2 -> C3 -> C4`
execution chain. The immutable Stage-1, mask/statistic, and prepared
acquisition parents retained their original hashes.

The Phase A child authorization has SHA-256
`8a9e4cd4427504fbe966136c3cb2a4b637c5c05c6f6ee7189801d1097e1f9ba0`.
Only coordinate GET and CCD snapshot GET were true. The 22-test acquisition
fault matrix and 121-test related Stage-1 suite passed before execution.

The first aggregate run exposed a Windows startup race while resolving a
not-yet-created output parent. Other workers had already completed 419
write-once bodies. Recovery found exactly two missing identities (`1at1`,
`1b0o`), 419 matching body/record pairs, no orphan and no `.part`. Re-execution
under the same unchanged child verified those 419 pairs and acquired the two
missing bodies plus CCD.

Phase A therefore passes with 421 unique coordinate bodies, 459 instances,
one CCD snapshot, 313,452,618 compressed bytes, and 1,360,881,606
decompressed bytes. Coordinate/CCD ledger hashes are
`61454d03f5bc22faa85a7389cadccfa051f2522d3185a4998f81ade7aeaf764b`
and `75076d1649e3458bb9b0666895b7dd6859acffc902b8272f439e1489968ce6c8`.
No coordinate or CCD was parsed and no event or affinity was read. The next
action is a separate hash-bound Phase B child for full A1-R C1 reliability.

Phase B was issued at SHA-256
`44b44369545e85181a2369c6b7ea44162638a150b133d948f99ad57327191b4f`
with only coordinate parsing and event extraction true. Its engineering smoke
exposed entity-sequence non-closure and a non-scientific OpenBabel InChIKey
plugin issue. A global 459-instance structure preflight was therefore run
before spending the full extractor budget.

Only 143/459 instances had an entity sequence exactly equal to the frozen
target sequence. Only 32/153 target units had all three exact-sequence
depositions, against the frozen minimum of 128 complete units. This is already
an upper bound before ligand resolution, metals, covalent links, FG mapping,
PLIP/ProLIF completion, or repeatability. The ledger at SHA-256
`352fcd538b5dd23b3eb267e42244cd43a24c43db34a40dab5c1a499ee9338f33`
therefore returns:

```text
STOP_C1_A1R_STRUCTURE_COVERAGE_BELOW_FROZEN_MINIMUM
```

No threshold, membership, entity-sequence rule, or model was relaxed after
seeing this result. The full PLIP/ProLIF run, C1 repeatability statistic,
mandatory controls, C2 coupling, C3 student, and C4 affinity were not executed
because the earlier frozen construction gate failed. The single next action
under this authorization is to end the A1-v2 route without model rescue.

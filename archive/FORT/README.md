# FORT

FORT is a target-level long-tailed drug-target affinity research program. Its
core task is abundant-to-scarce transfer: measure whether data-rich head
targets improve prediction, ranking, uncertainty, and experiment selection for
scarce recipient targets at the same support budget as no-transfer controls. The
former strict unseen-target k=5 task is retained as a secondary stress test,
not the active architecture objective.

The current authorized executable work is metadata-only audit. Gate D0 returned
`DATA_NOT_READY`; affinity and neural architecture training are disabled:

```powershell
D:\anaconda\envs\drug\python.exe main.py topology-audit --split train
D:\anaconda\envs\drug\python.exe main.py natural-tail-audit
D:\anaconda\envs\drug\python.exe main.py audit
D:\anaconda\envs\drug\python.exe main.py interaction-audit
```

`main.py` coordinates all runnable stages. Model code belongs in `model`,
protocol and preprocessing code belongs in `scripts`, experimental work enters
through `research`, reports remain in `reports`, and raw/processed data remain
in `dataset`.

The natural-tail audit found only 40 pKi recipients in its most optimistic
parent/document/assay-closed upper bound, below the frozen minimum of 50. The
strict time/source/scaffold/provenance roster has zero recipients, and all
candidates join into one source/homology dependency component. pKi, pKd,
pKd_app, and pIC50/IC50 are not pooled. The corrected A2S source-support control
is descriptive, the cross-fitted router is NO-GO, and no model training is
active.

`main.py a2s-baseline` fits the pooled source-only ridge and evaluates
recipient calibration, source routing, and an abstention gate on target-side
single-cold pKi/pKd episodes. It is a falsification baseline, not the final
architecture.

`main.py train` is fail-closed on the Gate D0 `DATA_NOT_READY` decision. The
retained `scripts/train.py` is a legacy strict k=5 baseline and is not an active
HTL-DTA runner.

`interaction-audit` is a TRAIN-only, model-free 2x2 exact-ligand audit. It
reports endpoint-separated rectangles, replicate-noise propagation, and
homology-pair bootstrap units. The next action is provenance-rich pKi source
acquisition and an outcome-blind D0 rerun, not another model experiment.

`topology-audit` is metadata-only and writes
`dataset/processed/htl_target_topology.v1.json`. It reads only the requested
`train` or `development` split and does not read affinity labels or confirmation
rows.

`natural-tail-audit` is metadata-only and writes
`dataset/processed/a2s_natural_tail_d0.v1.json`. Its binding decision is in
`reports/active/a2s_natural_tail_d0_decision_2026-07-31.md`.

`a2s-trace-q1`, `a2s-trace` and `a2s-trace-headroom` are the source-only A2S-TRACE
stage. They open the `fit` and `probe` roles of the balanced v2 lock only; the
`locked` source role and the A2S recipient roster are never requested. They need
the CUDA `drug` environment:

```powershell
D:\anaconda\envs\drug\python.exe main.py a2s-trace-q1
D:\anaconda\envs\drug\python.exe main.py a2s-trace --max-epochs 12
D:\anaconda\envs\drug\python.exe main.py a2s-trace-headroom
```

`a2s-trace-q1` measures *where* correctly assigned support labels carry
transferable ranking information, by varying only the support policy and the
support-query nearest-Tanimoto stratum. Its decision is
`INFORMATION_ADMITTED_IN_A_LOCAL_RELATION_STRATUM`
(`reports/active/A2S_TRACE_Q1_STRATUM_DECISION_2026-08-01.md`): the gain lives
entirely at nearest Tanimoto >= 0.55 and is null below 0.35 in every policy.

`a2s-trace` is the mechanism stage. TRACE meta-learns a label-free per-pair
transport reliability with zero target-specific parameters, exactly nesting
fixed Tanimoto KRR. Its decision is
`POSITIVELY_CONTROLLED_NULL_LEARNED_TRANSPORT_NOT_ADMITTED`
(`reports/active/A2S_TRACE_Q2_MECHANISM_DECISION_2026-08-01.md`): the learned
part adds -0.0001 CI [-0.0006, +0.0005] over a bar that includes one global
transport scale, while the same learner recovers +0.016 to +0.026 CI of an
injected pair-reliability signal. The one surviving positive is that global
scale itself (+0.009 to +0.010 CI over unscaled KRR); it is a baseline
parameter, not a mechanism, and every future comparison must grant it.

`a2s-mode-gates` is the pre-implementation measurement suite for the A2S-MODE
meta-adaptation route (`reports/active/A2S_MODE_MECHANISM_PROPOSAL_2026-08-02.md`),
which replaces support-similarity transport with a small discrete target state:

```powershell
D:\anaconda\envs\drug\python.exe main.py a2s-mode-gates
```

Its decision
(`reports/active/A2S_MODE_GATES_A0_A4_DECISION_2026-08-02.md`) is
`PREMISE_CONFIRMED; K_SHOT_INFERENCE_NOT_YET_IDENTIFIABLE_WITH_AN_UNSHAPED_DICTIONARY`.
Gate A0 finds per-target ranking headroom of +0.052 to +0.085 CI in **every**
relation stratum, including the ones where all transport operators measure zero;
Gate A1 shows a small discrete mode set carries part of it and is not a better
global ligand model; Gates A2/A4 show an unshaped k-means dictionary is not
separable from k<=5 noisy residuals even when the world is exactly the model.

`a2s-mode-generalization` then answers whether that object is worth building:

```powershell
D:\anaconda\envs\drug\python.exe main.py a2s-mode-generalization
```

Verdict `GENERALIZABLE_BUT_NOT_FEW_SHOT_REACHABLE`
(`reports/active/A2S_MODE_GENERALIZATION_DECISION_2026-08-02.md`). The per-target
head survives a within-target Murcko-scaffold-disjoint split (+0.052 CI
[+0.029, +0.075]), so the object is real. But source-target heads have a nearly
flat spectrum, and projecting a target's own head onto the top-2 source
directions retains **-6%** of its gain, so there is no small shared structure to
transfer; the protein embedding predicts nothing zero-shot; and the measured
label learning curve has its knee at **k ~ 10**, with the best k<=5 cell at
-0.003. This refutes both the shared-mode dictionary and the low-rank-code
family on this substrate, and reframes the open question as whether
meta-learning can move that curve left from k~10 to k~5.

`a2s-rip-r0` tests whether a *certified subset* of ranking interventions beats
applying the whole few-shot head:

```powershell
D:\anaconda\envs\drug\python.exe main.py a2s-rip-r0
```

Decision `SELECTION_CEILING_REAL; NOT_REACHABLE_FROM_AN_OBSERVABLE_MARGIN`
(`reports/active/A2S_RIP_GATE_R0_DECISION_2026-08-02.md`). A hindsight-selected
40 % subset beats the wholesale head by **+0.075 CI** at k=5 -- more than the
entire fully-supervised head is worth -- but the best observable statistic
reaches only AUC 0.555, and the implementable rule is indistinguishable from
random selection and from a magnitude rescale. Preregistered triggers P4 and P5
both fired, so A2S-RIP is retracted; the certification layer (cross-task
threshold transfer) passed and is carried forward.

`a2s-transfer-object` is Gate T0, the measurement that follows nine falsified
mechanisms. It trains nothing and asks the two questions all nine assumed:

```powershell
D:\anaconda\envs\drug\python.exe main.py a2s-transfer-object
```

Decision (**revision 2**, `NO_TRANSFERABLE_CHEMICAL_HEADROOM_OBJECT_IS_MEASUREMENT_CONTEXT`,
`reports/active/A2S_TRANSFER_OBJECT_GATE_T0_DECISION_2026-08-02.md`). All four
gates fail. **This gate is EXPLORATORY only** -- it evaluates on `probe`, whose
outcome was consumed once by PIRS and may not drive model selection.

The one robust, deterministic, basis-independent result: a chemistry-free
**document-mean oracle** scores +0.061 CI, **beating** the full 26-dimensional
per-target chemical head (+0.052). The accompanying provenance audit shows the
within-target scaffold-disjoint split is not leakage-free -- **91.1 %** of query
rows reuse a support-side document, **88.8 %** reuse a support-side assay, and all
52 targets share documents across the split. Once same-document pairs are scored,
the per-target chemical remainder (+0.029 [+0.005, +0.056]) does **not** clear its
0.005 admission threshold, and cross-target transfer of a fitted head is
**negative** (-0.018 [-0.044, +0.005]).

Revision 1 of this gate reported a positive transfer result. **It is retracted.**
Its basis was irreproducible (`torch.svd_lowrank` takes no generator; the verdict
flipped on replay), it omitted the same-document control on the transfer arm, and
its information-theoretic "closure" depended on an arbitrary definition. The basis
is now a deterministic eigendecomposition -- which affects **every** gate built on
it, A0-A4, G1-G4, R0 and HOTSPOT included.

Every A2S gate must now report the **same-document** contrast, and the
document-mean oracle is a mandatory control alongside the magnitude-matched and
random-selection controls.

`a2s-nea-preconditions` is the programme's terminal measurement:

```powershell
D:\anaconda\envs\drug\python.exe main.py a2s-nea-preconditions
```

Decision `NO_CHEMICAL_ADAPTATION_OBJECT_SURVIVES_SEPARATION_STOP_PROGRAMME`
(`reports/active/A2S_NEA_PRECONDITIONS_D0_N0_N1_DECISION_2026-08-02.md`). It runs
on source `fit` only, leaving `probe` and `locked` untouched. The same 93 targets
are evaluated under two splits differing only in what they separate. Under
scaffold-only separation the per-target chemical head scores +0.0610 and a
chemistry-free **document-mean oracle scores more** (+0.0671). Under simultaneous
scaffold + document + assay separation the head's advantage collapses **93 %** to
+0.0044 [-0.0161, +0.0242], while the document oracle measures **exactly zero** --
a structural self-validation, since on a document-disjoint split it can only
predict a constant. Base concordance is essentially unchanged between regimes, so
the separated task is not harder; the head's advantage evaporates.

Gate N0 quantifies why: **per-context offsets explain 68 % of all residual
variance**, with an offset SD of ~1.7 pKi -- larger than the target-specific
effects being sought -- and a material scale term, so the acting nuisance group is
the full affine group. Gate N1 shows coverage was never the obstacle (85.9 % of
passive k=5 draws contain a usable within-document contrast). The mechanism was
buildable; there was nothing for it to learn.

The honest bound is **not** "the effect is zero": the same-document point estimate
is +0.0123 and positive in both regimes, and resolving it would need **~445
homology components against the 92 available** -- roughly 4.8x this corpus.
Reopening requires more independent components (Papyrus 05.7, BindingDB), not a new
architecture.

`reports/active/A2S_META_ADAPTATION_MECHANISM_DESIGN_2026-08-02.md` (revision 2)
records the A2S-FBA design and **withdraws it as a Stage 1 proposal**: its premise
was the retracted transfer result, its entropy and harmlessness terms admitted
trivial collapse, and a method-specific novelty analysis places it inside existing
modular / mixture-of-experts meta-learning. Gate F1 does **not** run. The
replacement first action is Gate D0 -- rebuild the split with simultaneous target,
scaffold, document and assay separation, then re-measure what headroom survives.

The predecessor branch `research/a2s-conformational-free-energy-state-20260802`
stopped at its first real-data gate
(`reports/active/A2S_CFES_C0B_DECISION_2026-08-02.md`,
`CFES_C0B_SEMANTICS_NOT_ADMITTED_STOP_CFES`): a ligand-only predictor of
protein-ligand contact profiles beats the explicit ligand-by-pocket model by
0.092, and the tiny surviving cross term is matched by a parameter-matched
non-multiplicative residual and by frozen random projections.

An earlier failure opened branch `research/a2s-hotspot-sparse-20260802`
(`reports/active/A2S_HOTSPOT_BRANCH_CHARTER_2026-08-02.md`). Its theoretical
basis is the binding hot-spot principle plus Free-Wilson additivity: a target's
response head should be **sparse on a target-specific support**, which
retrodicts the flat head spectrum, the failure of pooled protein embeddings, and
the measured k~10 knee (`s log(d/s)` with the measured `s~8`, `d=26`). The
in-basis test already passed -- 8 coordinates retain 63 % of the head's gain
where 8 principal directions retain 16 %.

Read `task.md` for the active contract and `history.md` for the compact evidence
ledger. The complete A2S-CFRA design and falsification matrix is in
`reports/active/a2s_dta_master_design_2026-07-31.md`. External literature is
reference material only; no external model code or data split is imported into
the active implementation.

# A2S PIRS representation gate decision

Date: 2026-08-02  
Branch: `research/a2s-interaction-state-20260802`  
Decision: **`INTERACTION_STATE_REPRESENTATION_NOT_ADMITTED`**

Artifacts:

- `reports/active/a2s_interaction_state_gate_2026-08-02.json`
- `reports/active/a2s_interaction_state_gate_records_2026-08-02.parquet`
- `reports/active/a2s_interaction_state_gate_weights_2026-08-02.pt`
- `research/a2s_interaction_state_gate.py`
- `tests/test_a2s_interaction_state_gate.py`

Only source `fit` and source `probe` labels were opened. Source `locked` and
recipient labels were not requested. This was the single registered PIRS probe
run; the result is not eligible for architecture, threshold, seed, or capacity
tuning.

## Direct answer

**Is the mechanism generalizable? No.** A two-coordinate interaction response
state learned on fit components did not transfer to 52 scaffold-disjoint probe
targets in 50 independent components. The negative result includes a successful
synthetic positive control, sufficient component count, and exact dense-rotation
invariance.

**Does it retain significant potential or substantial value? No on this
substrate.** The full-support oracle over the learned coordinates failed before
the k-shot state estimator was considered. A learned support-to-state operator
cannot repair a representation whose oracle lower confidence bound is below
zero. R1 is therefore prohibited.

## Primary results

The primary endpoint is component-bootstrap CI gain. At k=5:

| Cell or contrast | Mean | 95% interval | Decision |
|---|---:|---:|---|
| learned-channel oracle minus base, pooled | +0.00460 | [-0.00137, +0.01023] | fail |
| learned-channel oracle minus base, Tanimoto <0.35 | -0.00077 | [-0.00997, +0.00767] | fail |
| EB state minus base, pooled | +0.00025 | [-0.00273, +0.00329] | fail |
| EB state minus base, Tanimoto <0.35 | -0.00284 | [-0.00707, +0.00082] | fail |
| correct minus label derangement, pooled | +0.00206 | [-0.00024, +0.00452] | fail |
| correct minus norm-matched wrong-target residual | +0.00179 | [-0.00042, +0.00403] | fail |
| correct minus protein transplant | -0.00089 | [-0.00371, +0.00175] | fail |
| segment minus ligand-only | -0.00157 | [-0.00443, +0.00126] | fail |
| segment minus pooled protein | +0.00072 | [-0.00273, +0.00405] | fail |
| segment minus frozen random coordinates | -0.00013 | [-0.00305, +0.00299] | fail |

The segment arm's absolute pooled CI was 0.55326 versus 0.55278 for the frozen
base. NDCG@10 changed by -0.00162, Spearman by +0.00192, and the RMSE contrast
favoured the base by 0.00417. Only 48% of components had positive mean CI gain.

The k=3 to k=5 mean was monotone, but neither point had a positive lower bound.
The apparent local-chemistry k=5 gain was +0.00393 CI with a positive lower
bound of +0.00052; this is below the registered 0.005 material threshold and is
exactly the regime already covered by TRACE/KRR. It cannot rescue the mechanism.

## Harness validity

The synthetic two-state control passed:

| Budget | Correct gain | 95% lower bound | Correct minus wrong | 95% lower bound |
|---:|---:|---:|---:|---:|
| 1 | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| 3 | +0.1556 | +0.1247 | +0.1609 | +0.1221 |
| 5 | +0.3069 | +0.2823 | +0.2933 | +0.2625 |

The k=1 rank edit was exactly zero. Nine focused tests passed, including
support-null behaviour, rotation invariance, protein-segment permutation
invariance, the source-role firewall, and nonzero gradients through the
analytic state solve into both the interaction output and segment path. The run
used CUDA on the RTX 4060 and completed in 1,183 seconds.

Training loss decreased while held-fit-component validation generally flattened
or deteriorated. This is compatible with task-specific representation overfit,
but it is not the binding reason for the stop: the separately held probe oracle
itself has no admitted headroom.

## Scientific interpretation

PIRS tested the last unresolved distinction after BIR/MODE/TCRS: whether a
query-dependent `phi(p,x)` learned from coarse local sequence segments could
turn sparse recipient measurements into one or two transferable ranking
coefficients. It did not.

The result adds four constraints to the sequential evidence chain:

1. protein-conditioned interaction coordinates are not sufficient merely
   because they are nonlinear and query-dependent;
2. coarse ESM segment attention does not identify a transferable binding
   environment from scalar affinity tasks;
3. the failure is upstream of support inference because full-support oracle use
   of the coordinates also fails; and
4. the only positive cell is chemically local and therefore does not escape the
   TRACE transport boundary.

HyperPCM, DrugBAN, and PSICHIC already establish protein-conditioned parameter
generation and sequence-derived interaction representations. PIRS could not
have claimed those ingredients as innovations. Its prospective contribution
was a measured, budget-matched support state plus a learned operator that beat
analytic inference. Since R0 failed, that contribution was never reached.

## Successor branch

No larger PIRS encoder, extra epoch, or relaxed gate is authorized. The next
branch changes the biological information object:

> Binding affinity is a population-weighted free-energy difference between
> physically localized bound and unbound protein-ligand conformational states,
> not a static compatibility between a ligand vector and coarse sequence
> segments.

The prospective state is a one- or two-logit **conformational free-energy
state**. Source structural data teach state-specific bound-minus-free response
channels; k=3 can identify at most one centered population contrast and k=5 at
most two. Query predictions are changed through bounded state-specific energy
surfaces, not support-label interpolation.

This principle is supported biologically by conformational selection and
population shift (`10.1038/nchembio.232`) and computationally by quantitative
protein reorganization thermodynamics (`10.1021/acs.jcim.4c01612`). PLINDER
provides leakage-aware apo/holo and predicted-state resources
(`10.1101/2024.07.17.603955`), while IPBind makes single-state bound-minus-free
encoding a mandatory non-novel baseline (`10.1109/OJEMB.2026.3667030`).

The successor must first pass label-blind structure/accession coverage and an
external apo/holo semantic gate, then a fit-component-only affinity
representation gate. Probe will not be reused for method selection. Source
`locked` and recipient labels remain sealed until a complete learned-operator
gate is frozen.

## Promotion status

There is no major breakthrough. Nothing is copied into `model/` or `script/`.
The gate code and artifacts remain research evidence only.

Artifact hashes:

- JSON content: `58e18568afdf6478347697a8671c70cbc130f4f9fd35ea23cd906cd24b93e2c8`
- records: `4def869f6a2190101fdfecf9fa35b9aa1335f6eb1ea39be82c2d323d3166766f`
- weights: `69bb76ad7586b7cada2a5a76de615207e8d2bf618a5dc1d04e93291858b0d4b8`


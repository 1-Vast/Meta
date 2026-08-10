# Phase 2B S5D evidence consolidation

## Terminal result

```text
D1 ligand-steering collapse ............... NOT CONFIRMED — the registered mechanism is FALSIFIED
D2 conditional estimand E1, E2, E3 ........ ALL FAIL
terminal verdict .......................... LIGAND_DIRECTION_COLLAPSE_NOT_CONFIRMED
S4R verdict ............................... UNCHANGED
heldout-B .................................. NOT CREATED AND NOT READ
affinity value reads ....................... 0
trainable parameters introduced ............ 0
```

S5D trained nothing and reused the frozen S4R checkpoints byte-for-byte. It did
not reopen the ligand representation route, which the S4R stopping rule closed.

## 1. The hypothesis this stage registered, and its falsification

S5D was written to test one mechanism for S4R's R3 failure. The registered
claim was that the estimator collapses every ligand difference onto
approximately one residue direction per protein, so that a foreign ligand pair
produces almost the same field and therefore almost the same AP.

**That claim is wrong, and the preregistration required saying so plainly
rather than reaching for a second explanation.**

On the 131 heldout-A constructs with at least three eligible pairs, the top
principal energy fraction `rho` of the mean-centred unit fields is:

| quantity | median |
|---|---:|
| `rho_dg`, unit ligand differences — the data-side upper bound | 0.4550 |
| `rho_graph`, candidate residue fields | 0.4793 |
| `rho_base`, baseline41 residue fields | 0.5758 |
| excess `rho_graph - rho_dg` | 0.0138 |
| true-versus-foreign field cosine, over 46,817 pairs | 0.4487 |

The registered rule required a median `rho_graph >= 0.80` and an excess over
the data-side bound of at least `0.10`. Observed: `0.4793` and `0.0138`. The
estimator adds essentially **no** collapse of its own beyond the collinearity
already present in the ligand differences, and only `10.7%` of constructs sit
above `rho_graph = 0.90`.

The median true-versus-foreign field cosine of `0.4487` makes the point
directly: swapping in a foreign ligand pair changes the residue field
substantially — the two fields share under half their direction. The estimator
really is steering on the ligand.

The `rho_base > rho_graph` ordering is a second, independent confirmation that
the S4R representation change did what the audit predicted: the graph-aware
statistic produces **more** diverse residue fields than the mean-pooled 41-D
baseline, not fewer.

## 2. The conditional estimand

D2 restricted every comparison to the residues that actually changed. Inside
that set both classes are pocket residues, so pocket membership cancels exactly
and non-parametrically — no gauge, no projection, no tuning. The estimator is
`ap_symdiff_conditional`, already implemented and registered in
`p2b_residue_residual.pair_metrics` and already aggregated by the parent Phase
2B runner; S3R and S4R computed it per pair and never aggregated it.

40,157 of the 46,818 primary pairs are eligible, across 107 closure components.
Median changed residues per pair `7.0`, median gain fraction `0.5000`, so the
conditional problem is balanced and its chance AP sits near `0.64`.

| arm | component-macro AP_cond |
|---|---:|
| candidate, graph-aware | 0.655030 |
| foreign ligand pair | 0.655470 |
| conditional chance | 0.643744 |
| baseline41, mean-pooled | 0.638830 |
| trained permuted-label learner | 0.628586 |

| Gate | delta | LCB95 | required | result |
|---|---:|---:|---:|:---:|
| E1 candidate - conditional chance | +0.011285 | -0.007749 | +0.05 | FAIL |
| E2 candidate - foreign pair | -0.000440 | -0.021814 | +0.03 | FAIL |
| E3 candidate - permuted learner | +0.026444 | -0.002977 | +0.03 | FAIL |

| non-gating contrast | delta | LCB95 |
|---|---:|---:|
| E4 candidate - baseline41 | +0.016199 | -0.005947 |
| E5 baseline41 - chance | -0.004914 | -0.023844 |

Every Gate fails, and E1 fails on its lower bound as well as its margin: the
candidate's conditional advantage over chance is not distinguishable from zero
at the closure-component level. `baseline41` is *below* conditional chance. The
foreign arm ties the candidate to within `0.0004`.

## 3. What the two diagnostics say together

They are more informative jointly than either alone, and the joint reading is
sharper than the S4R result was.

1. The estimator **does** convert ligand differences into distinct residue
   fields. True and foreign ligand pairs give fields with median cosine
   `0.4487`, so roughly half the field direction is ligand-determined.
2. That ligand-determined half carries **no** information about the labels.
   Under `AP_bidir` the foreign arm costs `+0.000644`; under the conditional
   estimand it costs `-0.000440`. Two different metrics, the same answer.
3. Under an estimand that removes pocket membership exactly, the whole
   above-chance effect disappears: `+0.011285 [LCB -0.007749]` against a chance
   level of `0.643744`.

So the failure is not that ligand information is lost on the way in, and not
that the estimator ignores it. Ligand information arrives, it moves the residue
field by a large angle, and the direction it moves in is unrelated to which
residues actually gained or lost contact. Point 3 is consistent with the S4R
`AP_bidir` gain having come from sign-agnostic pocket-like structure that the
two-dimensional `span{1, b^P}` gauge does not remove — but S5D did not register
a test of that, D1 falsified the mechanism it *did* register, and the
consistency is therefore offered as an observation, not as an established
result.

## 4. Numerical limitation, disclosed

`top_principal_energy_fraction` normalizes to unit vectors before centring.
Rescaled copies of one direction differ by about `1e-16` after normalization,
so a totally collapsed construct floors near `0.87` rather than at exactly
`1.0`; only bit-identical fields return `1.0`. The floor biases `rho`
**downward** in the degenerate limit, so it can understate collapse and cannot
manufacture it. It cannot have produced this verdict: total collapse still
clears the registered `0.80` threshold, and the observed median is `0.4793`,
nowhere near the degenerate regime. Both properties are asserted in
`tests/test_s7_l2b_phase2b_s5d.py`.

## 5. Governance

- Trains nothing. Zero trainable parameters introduced. No new representation,
  no new arm beyond the four registered, no capacity change.
- Heldout-B was not created and not read. R6 not opened. Affinity value reads
  `0`. DAVIS, KIBA, recipient, ChEMBL and BindingDB reads `0`.
- One label view opened: `heldoutA_residue_masks.json.gz`.
- **Heldout-A has now been consumed three times, by S3R, S4R and S5D.** Every
  number in this document is development evidence. None of it confirms
  anything, and the panel is correspondingly weaker as evidence each time.
  This is the third look, and the registered stopping rule forbids a fourth
  estimand variant on it.
- No threshold, seed, margin, arm or eligibility rule was changed after a
  statistic was read. One run, as registered.
- The S4R terminal verdict `REAL_RESIDUE_DIRECTION_STILL_NOT_IDENTIFIED` is
  unchanged and was not reinterpreted.

## 6. Remaining boundary

Nothing here authorizes a new stage. The representation route stays closed by
S4R, and the conditional-estimand route is now closed too: it was the
mathematically clean way to remove the pocket confound, and it found nothing.

What has changed is the shape of the remaining question. Before S5D one could
argue that ligand information was either being destroyed upstream or diluted by
the metric. Both are now measured and neither holds. The ligand signal reaches
the residue field at full strength and points somewhere biologically wrong. The
natural reading is that the missing ingredient is **correspondence** — which
ligand substructure sits against which residue — and that a pose-free
sequence-plus-2D estimand has no channel to supply it. Testing that would be a
separately governed information stage about geometry, with its own
preregistration, and it is explicitly not authorized by this stage.

No biological statistic is admitted to `z`. Affinity, selectivity, few-shot
sectioning, heldout-B, R6, DAVIS/KIBA/recipient labels, CSMO, Band and the
frozen operator `A(F,z)=K(B(z)F(z))` remain untouched.

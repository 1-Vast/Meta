# Preregistration — S7/L2B Phase 2B

## One ligand-conditioned residue residual head, supervised by the same-protein ligand differential

Stage identifier: `S7_L2B_PHASE2B_LIGAND_CONDITIONED_RESIDUE_RESIDUAL`

Written: 2026-08-10. Repository commit `623602e`.

**Authorizing evidence.** Phase 2A terminal verdict
`LIGAND_CONDITIONED_RESIDUE_SIGNAL_WITHOUT_EDGE_COUPLING`
(`report/s7_l2b_r0r/PHASE2A_VERDICT.json`), whose mandated next action is
exactly one item: *preregister one ligand-conditioned residue residual head*.
Phase 2A registration SHA-256
`4e01401d0468e3241bd05cde93b2a517919209d770f8f75ab471e42228f1b12e`.

**This document authorizes nothing by itself.** No Phase 2B run has occurred, no
Phase 2B code exists, and execution requires explicit authorization. It is
written now so that the design is fixed before any Phase 2B number is seen.

---

## 1. Why this head and not a pair-coupling head

Three Phase 2A measurements determine the design, and they point away from the
architecture one would naively reach for:

1. The labels **are** ligand-conditioned at the residue level: `ΔJ = +0.258`
   [LCB +0.234] over a measured replicate noise floor, with chemistry
   association `ρ = +0.322` [LCB +0.299].
2. Edge-level coupling is **not identifiable** — in the teacher (median
   `z = +0.413` against a degree-preserving null, threshold 2.0) or in B5
   (coupling `+0.0060` over its null, below the registered 0.01 margin).
3. The **additive** problem is nowhere near solved: B5 reaches `0.0698` of a
   well-posed label-fitted additive ceiling of `0.389`, i.e. 17.9%.

A pair-coupling head would optimise a term worth `0.011` while leaving `0.32` of
additive AP unclaimed, against labels whose edge structure is statistically
indistinguishable from its own margins. The residue-differential head targets
the one signal that Phase 2A proved is present.

## 2. The model — a residual only

The B5 protein-pocket prior is **frozen**:

```text
b_r(P) := the residue-marginal component alpha_r of the weighted additive
          projection of the sealed B5 logits, as computed in Phase 2A
```

`b_r(P)` is label-free, already materialised, already hashed, and is not
retrained, refit or rescaled. The only trainable object is a residual:

```text
logit p_r(P, L) = b_r(P) + delta_r(P, L)
```

`delta_r` is **one small low-rank bilinear interaction** between the existing
frozen residue-local state and the existing ligand-atom state:

```text
delta_r(P, L) = sum_{k=1..K} (U h_r)_k * (V g(L))_k ,   K <= 8
```

where `h_r` is the frozen ESM2-650M residue state already extracted and hashed
(`ada2413765...`), and `g(L)` is a permutation-invariant pooling of the existing
41-dimensional atom features. `U` and `V` are the only new parameters.

**Forbidden, and not present in the design:** any additional PLM, any
cross-attention stack, any typed-interaction head, any geometry branch, any
knowledge graph, any parallel SSL module, any affinity head, any PU loss, any
increase in the frozen ESM2 size, and any change to the atom branch.

### 2.1 Mandatory projections

Before use, `delta` is projected away from the directions it must not be allowed
to re-learn, within each complex:

```text
delta <- delta - Proj_span{ 1 , b(P) , c(L) * 1 } delta
```

removing (i) the constant direction, (ii) the frozen generic pocket prior, and
(iii) ligand-only / global interaction-mass directions. The achieved residual
orthogonality `||X^T delta|| / (1 + ||delta||)` must be `<= 1e-8`; violation is a
fail-closed contract error, not a reported number.

## 3. The objective — the differential, by construction

Primary supervision is the same-protein ligand differential, which removes the
generic pocket marginal by construction rather than by penalty:

```text
true       Delta y_r(P; L_a, L_b) = y_r(P, L_a) - y_r(P, L_b)
predicted  Delta delta_r          = delta_r(P, L_a) - delta_r(P, L_b)
```

Loss: binary cross-entropy on the **symmetric-difference residues** only, i.e.
residues `r` with `Delta y_r != 0`, labelled `+1` on `R_a \ R_b` and `0` on
`R_b \ R_a`. Residues in neither or both contribute nothing, because they carry
no differential information and would simply reintroduce the marginal.

Pairs are drawn only from **within one exact construct** (`seq_key`), so residue
indices are identical by construction, and only from training closure
components. All loss weights are frozen before any test scoring.

## 4. Split, inference unit, and data

- Split: the existing frozen protein-closure split, unchanged. Training pairs
  come only from training components. Scoring is on held-out A.
- Inference unit: the protein closure component. Atom–residue rows and residue
  pairs are never inference units.
- Held-out A carries 47,016 within-construct scaffold-distinct pairs and 254
  replicate pairs across 27 components with both pair types.
- No ligand graph shared with training enters held-out A (the existing
  disjointness filter is unchanged).

## 5. Metric

**Differential AUPRC.** For each held-out ordered pair `(L_a, L_b)` of one
construct, rank the symmetric-difference residues by `Delta delta_r` and compute
average precision against the label `+1 on R_a \ R_b`. Chance level is
`|R_a \ R_b| / |R_a symmetric-difference R_b|` and is reported per pair.
Summarised by component-macro; intervals by paired component bootstrap, 10,000
resamples, seed `20260901`.

Tie handling uses the exact tie-aware expectation registered in Phase 2A
amendments 01/02.

## 6. Registered ceiling — so a modest number is read correctly

The replicate floor of `J = 0.636` means the differential target carries
irreducible annotation noise. The **replicate oracle** is therefore computed and
registered as the ceiling: predict `Delta y_r(P; L_a, L_b)` using the mask of a
*replicate* of `L_a` in place of `L_a` itself. Phase 2B results are reported
both in absolute terms and as a fraction of that ceiling. A differential AUPRC
below the replicate oracle is not evidence of model failure and must not be
reported as such.

## 7. Gates — frozen now

All are component-macro with one-sided 95% lower bounds.

| id | contrast | margin |
|---|---|---:|
| `D1` | differential AUPRC − per-pair chance level | ≥ 0.05 |
| `D2` | differential AUPRC − wrong-ligand control (`delta` evaluated with a foreign ligand substituted for `L_a`) | ≥ 0.03 |
| `D3` | differential AUPRC − capacity-matched random head (identical parameter count, frozen random `U`, `V`) | ≥ 0.05 |
| `D4` | differential AUPRC − chemistry-shuffled control (ligand identities permuted within construct) | ≥ 0.03 |
| `D5` | non-inferiority: full pair AP of `b + delta` must not fall below sealed B5 by more than 0.005 | — |

Every gate requires LCB95 `> 0` in addition to its margin. No gate may be
lowered. Failure of `D1` terminates the stage; the sequence-plus-2D residue
route is then closed and no larger model may be substituted for it.

## 8. Module-participation audit — mandatory, fail-closed

A module present in the architecture but without a causal training contribution
must be rejected. Before any gate is read, the following are computed and
reported for `U` and `V`:

- gradient norms per parameter block across training;
- parameter movement from initialisation (relative Frobenius);
- activation variance of `U h_r` and `V g(L)` on held-out data;
- branch ablation: zeroing `V g(L)` must collapse differential AUPRC to chance;
- input shuffling: permuting ligands must collapse it to the `D4` control;
- gradient blocking: detaching `h_r` must change the outcome measurably;
- capacity-matched random control (`D3`).

If `delta` has no causal contribution, the stage is reported as such and the
head is rejected regardless of any metric.

## 9. What Phase 2B does not touch

Real ChEMBL/BindingDB affinity training; DAVIS, KIBA and recipient labels;
independent confirmation scoring; few-shot section adaptation and any `k`-shot
claim; admission of any statistic into production `z`; CSMO, Band, mesh,
positive ridge and `A(F,z) = K(B(z)F(z))`; P2–P4. All remain frozen. All Phase
2B code stays under `research/`.

A Phase 2B PASS would establish a **structural** ligand-conditioned residue
statistic and nothing more. It would not establish affinity semantics, and it
would not authorize a source-affinity Gate, which continues to require
closure-component OOF `correct − ligand ≥ 0.03` and `correct − wrong protein
≥ 0.03`, both with 95% lower bounds above zero, followed by a sealed transfer
Gate.

## 10. Registered expectations

1. The differential AUPRC is expected to be **modest in absolute terms**,
   because roughly 56% of the observed alternative-ligand mask difference sits
   at the replicate noise level. The fraction-of-ceiling figure is the one to
   read.
2. `D2` and `D4` are expected to be the discriminating gates, not `D1`.
3. A **large** differential AUPRC should first be suspected of construct leakage
   — the same PDB entry or a near-duplicate ligand appearing on both sides of a
   pair — and the pair census in `PHASE2A_CONSTRUCT_GROUPS.json` is the first
   artifact to inspect.

# PFSC-0 pan-family single-cold-start gate — decision

Verdict: **`PFSC0_FAIL_STOP`** / `PAN_FAMILY_SINGLE_COLD_UNDERPOWERED_UNDER_DOCUMENT_ISOLATION`.

Preregistration `reports/active/pfsc_preregistration.md`; result `reports/active/pfsc_gate.json`;
runner `research/pfsc_gate.py`. Train partition only; no development/confirmation/sealed label read.
Scope: **single cold-start** (novel target, known ligands) — a declared relaxation of the ligand-cold
axis, not a dual-cold claim.

## Result (pKi promiscuous block: 771 ligands × 187 targets × 172 homology clusters)

Under strict leave-homology-cluster-out + cross-family + **full doc-set-disjoint** neighbours, only
**19 independent homology-component units** were scorable (MDE80 0.064 ≫ the 0.03 threshold).

| arm | component-macro ρ (n=19) |
|---|---|
| ligand_only (B0) | −0.1017 [−0.3136,+0.1131] |
| esm_kernel | −0.0861 [−0.2936,+0.1160] |
| protein_shuffle | −0.0163 [−0.2111,+0.1788] |
| random_protein | −0.0171 [−0.2346,+0.1888] |
| far_kernel | −0.0937 [−0.3034,+0.1205] |

| contrast | mean | 95% | gate |
|---|---:|---|---|
| G1 esm − ligand_only | +0.0156 | [−0.1865,+0.2110] | FAIL |
| G2 esm − protein_shuffle | −0.0698 | [−0.3172,+0.1666] | FAIL |
| G2 esm − random_protein | −0.0690 | [−0.3230,+0.1843] | FAIL |
| G0 esm − far_kernel | +0.0076 | [−0.2467,+0.2449] | FAIL |

All five gates fail. Two independent reasons: (1) **underpowered** — 19 units, CIs ≈ ±0.3; (2) **no
positive point signal** — esm_kernel sits below its own shuffle/random controls, so there is not even a
directional hint of protein-proximity transfer on the resolvable units.

## Why: the document confound is the binding constraint

The block *looked* like a powered pan-family panel because its promiscuous ligands are measured across
many families. But that density is largely **within-document** (kinase/selectivity profiling papers
reuse a compound across many targets in one document). Strict, mandatory cross-document isolation — the
only way to exclude the same-protocol shortcut the whole program firewalls against — collapses the
usable substrate:

* first-doc approximation (understates the confound): 55 components;
* correct full doc-set disjointness: 28 components;
* correct isolation + leave-cluster-out CV: **19 components**.

**Preregistration disclosure:** PFSC-0 expected ~83 components from a feasibility audit that
approximated cross-document with a first-document check. The correct strict test is underpowered. No
threshold or isolation rule was relaxed to compensate; the stricter, correct condition was kept and the
route stopped.

## Decision

PFSC-1 (predictor + calibration) is **not authorised**: the identifying gates failed and the substrate
cannot resolve the effect under valid document isolation. No rescue is admissible — relaxing to
within-document neighbours would be the exact assay/protocol shortcut the contract prohibits.

This is the third independent confirmation, now at the single-cold-start resolution, that the binding
constraint is **document-independent factorial overlap in open continuous-affinity data**, not model
architecture or estimand:

* dual-cold cross-source (CROSSDOC): 11–13 units vs 25–30 required;
* dual-cold entity-disjoint 2×2 (BM2-PIRR P0): 160 blocks, sign-MDE80 0.111 > 0.10;
* single-cold-start pan-family (PFSC-0): 19 units, MDE80 0.064.

Program verdict remains **② `SIGNAL_PRESENT_EVIDENCE_INSUFFICIENT`**: a small real coarse protein
signal exists when document isolation is relaxed or the source is single (KirHub within-source +0.029;
dense-panel PB real−random +0.051), but under the correct document-isolated evaluation no open
substrate — dual- or single-cold — is powered to confirm it. The remaining route is a **prospective
measurement design**: a factorial panel where the *same* compound library is deliberately measured
across many multi-family targets under *independent* documents/assays, sized for ≥25–30
document-independent components. No architecture, estimand reframe, larger encoder, extra seed or
threshold change substitutes for that measurement.

`sealed_test_consumed=false`; `confirmation_labels_read=true` (pre-existing; PFSC read train only).

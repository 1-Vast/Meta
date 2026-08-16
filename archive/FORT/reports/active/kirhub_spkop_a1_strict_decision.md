# A1 SPKOP strict-firewall decision

Verdict: `A1_PROTEIN_NECESSITY_FAIL_STOP_A2`.

This is the authoritative A1 result. It supersedes the earlier KLIFS-family/Murcko-only run, which
did not explicitly close cross-scaffold high-Tanimoto edges or use the project-wide full-sequence
homology components.

The corrected split contains 324 full-sequence homology components and 79 chemical components.
Chemical components transitively join canonical parent identity, equal Bemis–Murcko scaffold or
Morgan Tanimoto >=0.50; the maximum similarity between distinct components is 0.478873. Five-fold
assignment is performed only after connected-component closure. Across all 25 fold pairs, the
minimum support is 59 target profiles, 51 homology components and a median of eight evaluable query
ligands per target. The estimator, seed, neighbours, arms, metric and gates are unchanged.

Component-macro Spearman:

| Arm | Spearman | LCB95 | UCB95 |
| --- | ---: | ---: | ---: |
| ligand-only | 0.0429 | 0.0236 | 0.0622 |
| true frozen ESM | 0.0719 | 0.0515 | 0.0920 |
| protein shuffle | 0.0261 | 0.0058 | 0.0470 |
| random protein | 0.0200 | 0.0002 | 0.0396 |
| KLIFS-group centroid | 0.0830 | 0.0628 | 0.1028 |

Across 308 evaluable homology-component units, true ESM minus ligand-only is +0.0290
[+0.0083,+0.0497]. It beats protein shuffle by +0.0458 [+0.0214,+0.0692] and random protein by
+0.0519 [+0.0283,+0.0755]. However, the point gain is below the preregistered minimum +0.030 and
true ESM does not beat the KLIFS-group centroid: -0.0110 [-0.0318,+0.0099]. The MDE80 envelope is
+0.0160 at paired SD 0.10, so the failure is not attributable to inadequate component count.

Interpretation: frozen ESM carries a real protein-dependent signal under strict dual cold, but the
effect is just below the minimum substantive threshold and is not demonstrably richer than coarse
kinase-group taxonomy. This supports final category ② (signal present, evidence insufficient), not
authorization of A2. No mechanism revision, extra seed or learned model is allowed.

The single KIRHub release cannot isolate assay, document or publication components. Therefore this
result is a within-source mechanism probe only. That common-source condition makes the task easier;
it cannot be used as positive evidence for cross-assay DTA generalization.

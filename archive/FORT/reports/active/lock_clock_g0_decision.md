# LOCK/CLOCK G0 decision

Accepted A1 artifact:
`reports/active/lock_clock_g0.json`, SHA-256
`555522b0956834b807972b002b7de79b01dba100a9680a7d4d449e6cb7f022e0`.

The accepted rerun is scientifically identical to the invalidated pre-A1 run after removing source
hash fields. A1 changed only typed artifact comparison.

## Result

G0-R used 353 activity-eligible genes, 303 frozen target components, 92 shared ligands, and the same
fixed top-eight squared-similarity estimator for every coordinate. Common ligand masks and non-finite
Spearman cells were removed atomically across all arms.

Primary group-mode component contrasts:

| contrast | mean | 95% bootstrap interval | n | MDE80 |
| --- | ---: | --- | ---: | ---: |
| LOCK - group centroid | `+0.02997` | `[+0.01223,+0.04799]` | 301 | `0.02558` |
| LOCK - aligned identity | `-0.01370` | `[-0.02304,-0.00459]` | 301 | `0.01301` |
| LOCK - composition | `+0.03565` | `[+0.02081,+0.05087]` | 301 | `0.02150` |
| LOCK - pooled ESM-2 | `+0.04548` | `[+0.02826,+0.06271]` | 301 | `0.02464` |
| LOCK - position shuffle | `+0.08327` | `[+0.06315,+0.10300]` | 301 | `0.02854` |
| LOCK - sequence shuffle | `+0.09184` | `[+0.07168,+0.11210]` | 301 | `0.02889` |
| LOCK - BLOSUM permutation | `-0.00161` | `[-0.01138,+0.00840]` | 301 | `0.01430` |
| LOCK - random PSD | `-0.00444` | `[-0.01553,+0.00684]` | 301 | `0.01590` |
| LOCK - matched wrong target | `+0.09403` | `[+0.07372,+0.11534]` | 301 | `0.02951` |

The primary gain missed the frozen `0.030` substantive bar by `0.0000257`. More importantly, fixed
LOCK was worse than exact aligned identity and did not beat either the BLOSUM-label permutation or
the correlation-matched random PSD control. These failures exclude a substitution-specific mechanism
claim even though LOCK beat composition, pooled ESM-2, position/sequence shuffles, and wrong targets.

Family-restricted component gates had 215 paired components and MDE80 `0.00921`; LOCK beat the uniform
family centroid by `+0.00767 [+0.00135,+0.01411]` and within-family wrong target by
`+0.01766 [+0.00933,+0.02643]`. This cannot rescue the group mechanism gate. At the family-macro level,
LOCK minus family centroid was `+0.00426 [-0.00024,+0.00888]` over only 44 families.

`conservation_LOCK` beat fixed LOCK by `+0.00664 [+0.00125,+0.01219]` and composition by
`+0.04229 [+0.02703,+0.05756]`. It is not CLOCK, has no structure embedding, is not lower-dimensional,
and was not allowed to rescue fixed LOCK. Its group-mode mean (`0.45374`) remained below aligned
identity (`0.46079`).

## Decision

Verdict: `LOCK_G0_REORDERING_NOT_IDENTIFIED_STOP`.

Final category: **3 - current data cannot identify the substitution-geometry mechanism; prospective
measurements are required.**

No predictive model, affinity regressor, GP, R-MAON, K-LBP, Transformer, or deep encoder was trained.
Strict dual-cold performance and cross-family/provenance transfer remain untested.


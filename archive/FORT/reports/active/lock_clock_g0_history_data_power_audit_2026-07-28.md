# LOCK/CLOCK G0 independent history, data, and power audit

**Date:** 2026-07-28
**Environment:** `drug`
**Scope:** read-only audit of label-free metadata and already-emitted result artifacts.
**Firewall:** KirHub Table S4 was not opened by this audit. No development, confirmation, Davis,
ChEMBL-confirmation, or sealed label was read.

## 1. Independent conclusion

G0-L is an admissible and genuinely new mechanism audit. Earlier work tested pooled ESM-2, exact
aligned pocket identity, pocket composition, target shuffles, random protein coordinates, and several
physical or learned-geometry routes. It did **not** test the fixed BLOSUM substitution geometry of LOCK
against position shuffle, amino-acid-label permutation, and a correlation-matched random PSD token
kernel.

The label-free substrate is large enough to run G0-L:

- 372 genes are in the human 85-position KLIFS pocket / frozen ESM-2 / strict-component intersection;
- they span 319 strict full-sequence homology components, 93 KLIFS families, and 9 KLIFS groups;
- 371 pocket strings are unique;
- all residues use exactly the registered 20 canonical amino acids plus the gap token.

The family-mode G0-R gate is plausible but cannot be certified before G0-R. Under the exact existing
greedy five-fold allocator applied to the 372-gene label-free universe, 270 query genes in 233 strict
components have at least two same-family training candidates. This is an upper bound before activity
eligibility, held-ligand missingness, and complete-case contrast pairing. It is not permission to replace
the preregistered requirements of at least 70 actually paired components and measured MDE at most 0.03.

The claim boundary remains narrow. Even a complete G0 pass would establish a single-source,
ligand-warm, within-kinome reordering mechanism. It would not establish strict dual-cold prediction,
cross-family transfer, provenance-independent replication, absolute affinity prediction, or true CLOCK.

## 2. Bound artifacts and audit inputs

| artifact | SHA-256 |
| --- | --- |
| `reports/active/lock_clock_g0_preregistration.md` | `0f5932ac130ecb70f127878f23d8b347794f959c1883a3924b3cc7a2a583c762` |
| `manifests/lock_clock_g0_sources.json` | `68d92261d0330ee57854bf09dadb28519ed6e15ae2b1d79dab114070ab58316b` |
| `dataset/public/kirhub_2026/processed/strict_components.json` | `0122780145a5504a61d78a709829d65d093f5636ccb4bca117d5e39cc2c4a7d4` |
| `dataset/public/kirhub_2026/processed/target_esm2.npz` | `062e4accc076527014adc5b4865788eed31a17ff43d1feff2b9e883b5c0cff79` |
| `dataset/public/klifs_2026_07_22/raw/kinase_information.json.gz` | `9a64e3e9027aeb729812fcdb244b07c65e162b7c68d3e2785b939e36c05eb4d9` |
| `dataset/public/kirhub_2026/processed/h0_component_registry.json` | `fc3bd42fb7618c548ad7b661fa7907cfacd7afa70a768782c99ea83edf1825ed` |

The implementation amendment binds the first two hashes and supplies the missing normalization,
composition-projection, random-PSD, derangement, and measured-MDE definitions. This audit does not
change them.

Existing result JSON and decision reports were used only as historical evidence. Raw activity cells
were not re-read or reconstructed.

## 3. Label-free census

The intersection was formed by gene name using the same pocket rule as the prior ASPIRE loader:
retain the first human KLIFS record with a pocket string of length 85, then intersect with the ESM key
set and `strict_components.target_axis.component_by_gene`.

| quantity | count |
| --- | ---: |
| human KLIFS genes with an 85-position pocket | 521 |
| frozen ESM-2 genes | 377 |
| strict-component genes | 377 |
| three-way G0-L intersection | **372** |
| strict components represented | **319** |
| KLIFS families represented | **93** |
| KLIFS groups represented | **9** |
| unique pocket strings | **371** |
| within-family unordered gene pairs | **1,300** |

The five ESM/strict genes without a valid 85-position pocket are `EEF2K`, `PDK2`, `PDK3`, `PDK4`,
and `TRPM7`. The only duplicated pocket string belongs to `PAK1` and `PAK3`.

Group sizes are:

| group | genes |
| --- | ---: |
| AGC | 50 |
| Atypical | 1 |
| CAMK | 65 |
| CK1 | 11 |
| CMGC | 49 |
| Other | 47 |
| STE | 41 |
| TK | 78 |
| TKL | 30 |

The pocket alphabet is exactly
`-ACDEFGHIKLMNPQRSTVWY`. There are 80 gap occurrences among 31,620 aligned positions and no unknown
token requiring an after-result encoding rule.

### 3.1 Candidate-set feasibility

For this audit, strict components were assigned with the existing deterministic allocator:
components are ordered by decreasing represented gene count, ties by component ID, and placed in the
least-loaded of five folds, ties by fold number. This gives label-free fold loads of
`75 / 75 / 74 / 74 / 74` genes.

| candidate restriction | queries with at least 2 training candidates | strict components represented | candidate count median | range |
| --- | ---: | ---: | ---: | ---: |
| same KLIFS group | 371 / 372 | 318 / 319 | 40 | 0-64 |
| same KLIFS family | 270 / 372 | 233 / 319 | 3 | 0-25 |

The one group-mode exception is `MTOR`, the only Atypical query. In family mode, 50 genes have zero
training candidates and 52 have one under this fold allocation. These queries must not be silently
promoted to a group or global candidate set.

Per-fold family-mode upper bounds are:

| held target fold | query genes | strict components | genes with at least 2 candidates | represented components |
| --- | ---: | ---: | ---: | ---: |
| 0 | 75 | 64 | 53 | 47 |
| 1 | 75 | 64 | 50 | 43 |
| 2 | 74 | 63 | 52 | 45 |
| 3 | 74 | 64 | 61 | 51 |
| 4 | 74 | 64 | 54 | 47 |

These counts are label-free feasibility bounds, not G0-R sample sizes. Existing ASPIRE retained 353
activity-eligible genes with valid pockets and produced 302 evaluable strict components in group mode.
Family restriction and complete-case pairing will reduce that number further.

### 3.2 Derangement feasibility

On the 372-gene label-free universe:

- within-group sequence derangement has one unavoidable gene fixed point, the singleton Atypical block;
- within-family wrong-target derangement has 18 unavoidable gene fixed points from 18 singleton
  families;
- every non-singleton group or family can be cyclically deranged with zero fixed points;
- the matched `(group, eligible-ligand quartile)` fixed-point count is intentionally unknowable before
  G0-L passes because the quartile depends on Table S4 non-saturation counts.

The runner must recount and report fixed points on the final G0-R eligible universe. The label-free
counts above cannot be copied into the result after activity filtering.

## 4. Power audit

The amendment defines measured MDE on complete-case paired component deltas:

```text
MDE80 = (z_0.975 + z_0.80) * sd(delta_component) / sqrt(n_components)
      = 2.8016 * sd(delta_component) / sqrt(n_components)
```

The earlier KirHub reports used an envelope with an assumed paired SD of 0.10. At that planning value:

| paired components | MDE80 at SD 0.10 |
| --- | ---: |
| 70 | 0.0335 |
| 92 | 0.0292 |
| 233 | 0.0184 |
| 302 | 0.0161 |
| 308 | 0.0160 |

Therefore the two family gates are independent. Exactly 70 paired components pass the count gate but
fail the 0.03 MDE gate when paired SD is 0.10. At 70 components, empirical paired SD must be at most
0.0896; at the label-free upper bound of 233 components it may be at most 0.1635.

The execution must:

1. average the five ligand-fold Spearman values within target;
2. average target values within strict full-sequence component;
3. form each contrast on components observed in both arms;
4. compute both the bootstrap interval and empirical-SD MDE on that same paired vector.

Fold cells, targets, and ligands are not independent power units. An arm-specific total component count
cannot be used for a contrast when its paired complete-case intersection is smaller.

## 5. Historical identification matrix

| question | strongest existing evidence | audit reading |
| --- | --- | --- |
| Does pooled ESM add target-conditioned signal? | A1 strict: true ESM minus ligand-only `+0.0290 [+0.0083,+0.0497]`, but true ESM minus group centroid `-0.0110 [-0.0318,+0.0099]` over 308 components. | A small single-source signal exists, but pooled ESM is not richer than coarse kinase coherence. |
| Is the signal specifically KLIFS-group taxonomy? | TR-0: group minus ligand-only `+0.0400 [+0.0181,+0.0616]`; group minus leave-own-group centroid only `+0.0091 [-0.0076,+0.0256]`. | Own-group identity is not load-bearing. A generic kinase centroid retains most of the gain. |
| Does an aligned pocket help a shared-panel oracle? | ASPIRE: aligned identity `0.4539`, group centroid `0.4253`, pooled ESM `0.3983`, pocket shuffle `0.3562`, random same-group target `0.3670`; pocket minus group `+0.0286 [+0.0116,+0.0459]`. | Aligned pocket information exists, but the registered substantive threshold was missed by 0.0014. |
| Is aligned position stronger than pocket composition? | PARC diagnostic: positional pocket minus composition median `+0.0079`, positivity `0.53`, Wilcoxon `p=0.857`; only 77 evaluable components and the inherited containment mean is invalid. | Suggestive evidence against a positional increment, not a decisive refutation. G0's direct geometry controls remain warranted. |
| Do generic random-coordinate controls explain prior gains? | A1 true ESM beat protein shuffle and random protein. ASPIRE aligned identity beat within-group pocket-target shuffle and random same-group targets. | Random coordinates are partially excluded, but no prior arm matched BLOSUM token correlation or permuted amino-acid labels. |
| Did learned protein geometry survive corruption elsewhere? | SCGD: true minus protein shuffle `-0.0002 [-0.0023,+0.0016]`; true minus random protein `+0.0049 [-0.0035,+0.0182]`; protein-free trained arm was better. CFRI true minus random target coordinate `-0.0260 [-0.1070,+0.0447]`. | Learned geometry on pooled ESM repeatedly failed specificity. This does not test fixed LOCK. |
| Did structure-mechanism supervision establish both axes? | MIF-NK: joint minus target shuffle `+0.0256 [+0.0130,+0.0384]`, below 0.03; joint minus ligand-only `+0.0052 [-0.0104,+0.0211]`; conformation heads failed. | The high contact AP was largely a ligand/contact-frequency prior, not identified target conditioning. |
| Is there provenance-independent cross-family support? | PFSC scored only 19 document-isolated components, MDE80 0.064; ESM minus random protein was `-0.0690 [-0.3230,+0.1843]`. | No. KirHub remains single-source and kinase-only. |
| Can ordinary physical structure substitute for the missing geometry? | Deployment docking `rho` about 0.005; holo versus AlphaFold delta `-0.019`; native-pose joint `0.257` below ligand marginal `0.394`; adding physics changed descriptor base `0.381 -> 0.125`. | The closed pose/contact route mostly captured ligand size/composition and does not answer the substitution-kernel question. |

### 5.1 What is actually new in G0

ASPIRE's similarity was the fraction of exactly matching aligned positions. It did not assign partial
credit to conservative substitutions. PARC changed the coordinate in a low-rank containment model and
did not compare:

- fixed LOCK versus exact aligned identity;
- fixed LOCK versus BLOSUM-aware composition;
- fixed LOCK versus within-target position shuffle;
- fixed LOCK versus amino-acid-label-permuted BLOSUM;
- fixed LOCK versus a token-correlation-matched random PSD kernel.

Consequently the fixed LOCK geometry is not a reparameterization of a previously failed test. It is a
new, bounded mechanism hypothesis. Conversely, `conservation_LOCK` is not CLOCK: no learned positional
structure embedding or mutation-landscape likelihood is present.

## 6. Execution invariants and pre-result risks

1. **Stage order is binding.** G0-L may read only KLIFS, ESM, strict components, and label-free
   taxonomy. Table S4 and the matched eligible-ligand-quartile derangement may be constructed only after
   the complete fixed-LOCK G0-L gate passes.
2. **The fold assignment needs one deterministic universe.** `strict_components.json` stores component
   IDs but no persistent fold numbers. Earlier KirHub runners recomputed balanced folds after their own
   eligibility filtering. The G0 implementation must record the exact gene/component universe and the
   resulting component-to-fold map before calculating any G0-R arm statistic; all arms and controls must
   reuse it.
3. **Family failure is not backfilled.** Queries with fewer than two same-family training candidates
   are absent from family mode. They cannot use group candidates, and their absence must propagate
   through complete-case pairing.
4. **Fixed points are result metadata.** Singleton-block fixed points are retained and counted, not
   dropped after seeing performance. The final count must be recomputed after G0-R eligibility filtering.
5. **Measured MDE uses observed paired dispersion.** The historical SD-0.10 envelope is context only.
   It cannot replace the amendment's empirical paired-delta SD.
6. **No dimension rescue.** Failure of 16 dimensions to retain 0.80 centered energy does not fail the
   main G0-L execution gate and does not authorize increasing rank.
7. **No predictive wording.** Training-target activity at query ligands is used by the oracle. Any pass
   is `TRAIN_ONLY`, within-kinome, and single-source.

## 7. Audit recommendation

Proceed with G0-L under the frozen preregistration and implementation amendment. It has a valid
label-free substrate, a distinct untested mechanism, and controls that directly separate substitution
geometry from taxonomy, composition, exact identity, and random token geometry.

Proceed to G0-R only if the fixed-LOCK G0-L gate passes exactly as written. The label-free census makes
the family count gate feasible in principle, but neither the final paired-component count nor the
measured MDE can be declared in advance. Failure of either must produce
`LOCK_G0_REORDERING_NOT_IDENTIFIED_STOP`; it cannot be repaired by changing top-k, folds, kernel scale,
rank, controls, thresholds, or candidate restrictions.

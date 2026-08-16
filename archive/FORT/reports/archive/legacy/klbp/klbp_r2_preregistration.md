# K-LBP v2 R2 — LLM Mechanism Card compilation and closed-book contamination audit (preregistration)

**Frozen:** 2026-07-27, **before any API call was made.**
**Program:** `task.md` Part 9 (K-LBP v2). **Stage:** R2. **Gating:** for card eligibility only. An
eligible card enters R1-part-2 (identical thresholds) and then, if it survives, R4. R2 authorizes no
affinity claim and no training.
**Authorization:** the user explicitly authorized API budget for the experimental-verification round on
2026-07-27 ("允许调用env文件的API"). This stage is the only stage in the ladder that spends it.
**Design rationale:** `reports/active/model_blueprint_reconstruction_2026-07-27.md` §8; task.md §9.6.

---

## 1. The two questions

> **Q1 (contamination):** With the source document withheld, can the compiler recover actual Metz
> panel train-cell affinity values from memory? (C2 closed-book probe — Busch et al. 2026 blinding
> methodology: progressive masking + cell-level power.)

> **Q2 (name memory):** Does the named card differ from the sequence-only de-identified card in a way
> that indicates target-name memory rather than sequence-derived content? (C1.)

A closed-book recovery **kills the named-card route** for these targets. A failed probe does **not**
prove absence of contamination (registered honestly).

## 2. API configuration (frozen)

```text
endpoint        = env MOONSHOT_BASE_URL (https://api.moonshot.cn/v1), OpenAI-compatible chat/completions
model           = env KIMI_MODEL (kimi-k2.6)          # version string recorded verbatim per call
temperature     = 0
top_p           = 1
response_format = JSON object where supported; else strict JSON system instruction + code-side parse
seed            = 1729 where the API exposes it; recorded otherwise
k_runs          = 1 (single run; k-run agreement is LEXOR's EBE-C concern, not R2's)
budget cap      = <= 700 calls and <= 2.0M total tokens (input+output); exceeding the cap STOPS the
                  stage rather than expanding scope (task.md §6.6.7 discipline)
key handling    = read from env at runtime; never written to any file, log, JSON, or report
recording       = for every call: prompt_hash (SHA-256 of exact prompt), model string, timestamp,
                  input/output token counts, raw response text (verbatim), parse outcome
```

## 3. Card compilation (C1 family)

Targets: the R1 eligible set (≈ 111 Metz train targets). Two compilations per target:

| card arm | prompt content | role |
| --- | --- | --- |
| `card_named` | target **name + UniProt accession + full sequence** + the frozen schema (task.md §9.4) + field definitions + "abstain when unknown; every non-missing field requires a resolvable source id (PMID/UniProt/PDB); no affinity values, no compound/chemotype names, no family labels, no source counts" | the Task B candidate coordinate |
| `card_deidentified` | **full sequence only** — name, accession, organism, and any literature withheld + identical schema and rules | the Task A candidate coordinate AND the C1 de-identified control |

Code-side validation of every card: JSON parses; fields within the frozen enum/ordinal ranges; no
field contains a number interpretable as an affinity value (regex guard: pKi/pKd/Ki/Kd/IC50/Kdapp,
nM/µM magnitudes, log values in [4, 12] adjacent to such tokens); no compound or chemotype names
(guard list); every non-missing field's source id resolves in format (PMID digits / UniProt accession
pattern / PDB id pattern); violations → card **rejected** (not repaired); missing fields → explicit
missing indicator. Compile date, model string, and prompt hash stored per card; a model or prompt
change voids every card. The accepted/rejected/missing counts are reported per arm.

## 4. Closed-book affinity probe (C2)

**Cell sample (frozen before any call):** from the Metz TRAIN registry, stratified by target exposure
and ligand, seed 1729: **100 real cells** (target, ligand, true pKi — the true value is recorded in
the sealed runner input and never shown to the model). **50 fake cells**: real ligands from the same
registry re-paired with targets that have **no** co-measurement anywhere in the panel (verified
against the registry; these pairs are chemically realistic but have no published Metz value).

**Prompt arms (each cell, temperature 0):**

| arm | content |
| --- | --- |
| P1 named | "Report the pKi of compound {SMILES} against {target name} ({accession}) from the Metz et al. 2011 kinase panel. If you do not know, answer NA." |
| P2 de-identified | identical, with target full **sequence** replacing name/accession |
| P3 progressive mask (real cells only, P1 ladder) | name → name with residues masked → accession only → sequence only; records the masking level at which any recovery disappears (Busch progressive-masking method) |

**Recorded per call:** produced value (float or NA), refusal, raw text.

**Statistics (frozen):**
1. **Recovery rate** `rr`: fraction of answered cells with `|pred − true| ≤ 0.30` log (the
   Landrum–Riniker maximal-curation noise floor; 0.45 sensitivity also reported).
2. **Constant-predictor control:** recovery rate of predicting the train-cell global mean for every
   real cell. The model's `rr(real)` must exceed this to count as recovery.
3. **Real-vs-fake contrast:** `rr(real) − rr(fake)` with a two-proportion z-interval; and a binomial
   test of `rr(real) − rr(fake)` against 0.
4. **Rank recovery:** Spearman(pred, true) on answered real cells vs 10,000 label permutations
   (permutation p; p < 0.01 = recovery).
5. **Named−deidentified gap:** `rr(P1) − rr(P2)` with interval; and per-field named-vs-deidentified
   card disagreement rate (C1).
6. **Power statement:** at n=100 real cells, the design detects `rr(real) − rr(fake) ≥ 0.18` at 80%
   power (one-sided α=0.05) — computed and reported.

## 5. Frozen gates

| gate | requirement |
| --- | --- |
| **G-R2-1 (closed-book)** | `rr(real)` exceeds the constant-predictor control **and** `rr(real) − rr(fake)` LCB95 > 0, **or** permutation p < 0.01 → **FAIL = contamination found** |
| **G-R2-2 (name memory)** | named−deidentified `rr` gap LCB95 > 0 → named arm ineligible (name-memory channel); de-identified arm may remain eligible |
| **G-R2-3 (schema discipline)** | card rejection rate ≤ 30% per arm; zero affinity-numeric guard hits accepted (any accepted card containing an affinity numeric → arm ineligible) |
| **G-R2-4 (budget)** | within cap (exceeding = stop, no result) |

```text
G-R2-1 fails on real cells            -> R2_CLOSED_BOOK_RECOVERY_STOP      (named route dead for Metz;
                                                                            de-identified route assessed separately)
G-R2-1 passes, G-R2-2 fails           -> R2_NAME_MEMORY_STOP               (named arm ineligible)
G-R2-1/2/3 pass                       -> R2_CARDS_ELIGIBLE                 (cards enter R1-part-2)
G-R2-4 exceeded                       -> R2_BUDGET_EXCEEDED_NO_DECISION
```

## 6. Declared expected outcome (stated before any call)

Metz 2011 is a famous public panel, likely present in pretraining corpora. The expected outcome is
**partial closed-book recovery on named arms** (P1 rr materially above the constant predictor on
well-known kinases and famous ligands — staurosporine-class), with the de-identified arm near the
constant predictor — i.e., the named-card route is expected to be **dead for Metz** and the
sequence-only card to be the surviving Task A coordinate. The named-vs-deidentified card disagreement
is expected to concentrate in functional-state fields (A) rather than site-chemistry fields (B).

## 7. Prohibited rescues

No threshold, budget, or prompt change after seeing a result. No second batch at a different
temperature. No re-sampling of cells after seeing outcomes. A G-R2-1 failure may not be argued away
by card-side controls. No development/confirmation/sealed label may be probed (real cells are TRAIN
only). The API key must never appear in any artifact. A failed probe is not evidence of cleanliness;
it is recorded verbatim.

## 8. Artifacts

```text
research/klbp_r2_cards_and_probe.py           runner, deterministic seed 1729, budget-capped
reports/active/klbp_r2_cards.jsonl.gz         compiled cards + validation record (SHA-256 in decision)
reports/active/klbp_r2.json                   probe result, parses with allow_nan=False
reports/active/klbp_r2_decision.md            verdict + contamination table + what was NOT shown
tests/test_klbp_r2.py                         schema guards, fake-cell non-co-measurement proof,
                                              prompt-hash stability, budget-cap enforcement (mock client)
```

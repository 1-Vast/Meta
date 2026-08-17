# Superseded (2026-08-16)

This cycle's four narrative documents — `RESEARCH_SYNTHESIS.md`,
`DATAFLOW_AUDIT.md`, `CANDIDATE_COMPARISON.md`, `PREREGISTRATION_DRAFT.md` —
were removed after consolidation. Two of their four load-bearing conclusions
did not survive re-testing, and keeping the originals alongside their
corrections invites citing the withdrawn version.

Read instead: `../a2_exact_probe/FINAL_DECISION.md`.

| this cycle claimed | outcome |
|---|---|
| ordering is interaction-borne (+0.186) | ✅ holds, reconfirmed across five donor strata |
| the interaction branch's ordering is protein-inert | ✅ holds, strengthened |
| the architecture can express what training removes (E3) | ❌ **withdrawn** — one random init; ten move in uncorrelated directions |
| the collapse is in the readout (E4) | ❌ **relocated** to the fusion/pooling inside `ContactGrammar` |
| A2's premise is falsified by the endpoint scalar | ❌ over-reach; superseded by the exact operator in Stage R, which closes A2 on its own gates |
| CPC recommended as the next stage | superseded by `../a2_readiness_v2/PREREGISTRATION_V2.md` Stage P |

Retained here because it is still valid and not duplicated elsewhere:

* `LITERATURE_LEDGER.md` — the primary-literature review;
* `tests/test_centering_excludes_the_level.py` — 11 probes verifying that a
  centered protein counterfactual gives `protein_head` identically zero
  gradient. Stage P depends on this property;
* `branch_ordering_probe.py`, `attention_locus_probe.py`, `_arms.py` and their
  two `.json` results — the measurements the corrections were made against.
  Their `"seal"` strings predate the 2026-08-16 terminology repair and describe
  logical exclusion after parsing, not a physical label seal.

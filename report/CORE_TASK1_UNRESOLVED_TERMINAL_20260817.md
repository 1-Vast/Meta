# Core Task 1 terminal verdict — UNRESOLVED on current local assets

Date: 2026-08-17. Machine-readable:
`CORE_TASK1_UNRESOLVED_TERMINAL_20260817.json`.

## Verdict

**UNRESOLVED on current local assets.**

This is a valid terminal verdict under the Core Task 1 definitions: the data
do not support a SOLVED claim, and the missing positive control plus severe
censoring forbid a FALSIFIED-AS-TESTED claim. It is explicitly **not**
evidence of biological absence.

## Why UNRESOLVED

1. **W0-P failed and cannot be repaired locally.**
   - Local surrogate: 6 near-identical BindingDB sequence pairs, 32 rows.
   - Leave-one-pair-out, 3 seeds, low-capacity bilinear model:
     correct positions sign accuracy 0.240; random positions 0.156;
     BLOSUM-approximate unrelated positions 0.125; global pooled 0.760
     (recorded as under-powered and unexplained, not a pass).
   - Local search cannot enlarge the panel beyond 7 pairs / 34 rows.
   - The required standard panel (resistance/gatekeeper/ortholog mutations
     with matched ligand panels) is absent locally; acquisition contract is
     frozen in
     `tools/research/stageW0b_core1_audit/W0P_ACQUISITION_SPEC.md`.
2. **Label censoring is severe** on the available single-platform panels:
   Davis 71.2%, Metz 60.4%, Klaeger 93.5% at detection floors.
   After censored exclusion, broad all-pairs layers still have support, but
   strict MMP layers remain small (Davis 7 classes / EIU 7; Metz 480 / EIU 9;
   Klaeger 37 / EIU 9).
3. **Cross-platform residual transfer** is positive only between Metz and
   Klaeger (+0.642 [0.482, 0.777]). Pairs involving Davis are negative; per
   the frozen interpretation this closes direct cross-platform residual
   sharing with Davis, not single-platform signal.
4. **Stage W W1 was never trained** and no W1 training metric was read; its
   local operator and data banks are retained as paused artifacts.

## Authorized reopening path

1. Acquire a standard W0-P panel per `W0P_ACQUISITION_SPEC.md` with full
   provenance (source, version, license, SHA-256, label semantics, censoring
   fields).
2. Re-census single-platform support after censored exclusion or
   interval-censoring.
3. Freeze a new W0-P/W1 preregistration from the enlarged panel's effective
   sample size and execute the full control matrix.

## Governance

No threshold was moved. Stage W / W0b / W0P preregistrations are frozen and
unchanged. `meta_test` was never evaluated. `model/` and production `scripts/`
were not modified by this investigation.

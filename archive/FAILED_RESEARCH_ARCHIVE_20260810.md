# Failed research archive index

This cleanup removes 259 superseded tracked experiment implementations,
route-specific tests and verbose reports from the active tree. Their state
immediately before deletion is fully
recoverable from Git commit:

```text
c05d3f95fe59f1f0b1e1cc34163ba473f16ea008
```

Removed active-tree families:

- `research/s7_l2b_r0r/` and `report/s7_l2b_r0r/`;
- `research/ssl_b2_structural_observability/` and its report tree;
- `research/correspondence_router/` and its report tree;
- superseded S7/SSL/correspondence regression tests;
- the failed ChEMBL X1A/X1A-R implementation and reports;
- superseded cycle-feasibility code and reports;
- terminal `crossed_panel_deployability` and `ssl_gauge_fixed` remnants.

The durable findings are consolidated in `history.md` and
`report/EXPERIMENTAL_EVIDENCE_LEDGER.md`. Current BindingDB quotient code,
T-BASIS feature generation, machine evidence, production `model/` and passed
data utilities remain active.

Untracked user artifacts were not added to this archive or to the cleanup
commit. Large third-party data and caches remain governed by
`DATA_AVAILABILITY.md`.

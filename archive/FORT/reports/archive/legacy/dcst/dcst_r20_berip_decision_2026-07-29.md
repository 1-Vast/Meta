# DCST-R20 balanced ERIP decision

Date: 2026-07-29  
Decision: `STOP_BERIP_BALANCED_MANIFEST_INADEQUATE`

## Result

Only the cross-homology and firewall gates passed. Endpoint consistency and
document separation exposed that the apparent R19 rectangle count did not
translate into a broad, endpoint-balanced Stage-1 substrate.

| Stage | Rectangles | Target pairs | Targets | Homology components | pKd |
| --- | ---: | ---: | ---: | ---: | ---: |
| Valid before caps | 324,444 | 233 | 145 | 136 | 340 |
| Target-pair cap | 11,871 | 233 | 145 | 136 | 256 |
| Final balanced | 11,871 | 233 | 145 | 136 | 256 |

The selected set was strongly cross-homology (`99.5367%`; 223 target-pair
blocks), but:

- valid topology missed the frozen 500 target-pair and 150-target floors;
- balanced scale was below 25,000;
- the largest target/homology-pair block still contributed `2.1565%`, above
  the 1% target-pair limit;
- pKd contributed only `2.1565%` (256 rectangles), below both endpoint gates.

No label-free manifest was written because the gate failed.

## Firewall

Only TRAIN metadata/reliability summaries were loaded. Numeric affinity was
not requested or loaded. Development, confirmation, and sealed rows remained
untouched. The CPU-only audit completed in `3.859 s`.

## Consequence

The R19/R20 exact-rectangle family is stopped. Its raw combinatorial count was
pseudo-replication rather than sufficient independent endpoint-balanced
information. A future same-domain Stage-1 experiment must use a different
information unit and carry explicit target-dependence controls; it may not
relax the R20 caps or mix pKi/pKd contrasts post hoc.

Authoritative machine result:
`reports/active/dcst_r20_berip_seed1729.json`.

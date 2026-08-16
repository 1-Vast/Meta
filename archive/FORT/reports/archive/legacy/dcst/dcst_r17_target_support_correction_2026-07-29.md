# DCST-R17 target-support correction

Date: 2026-07-29  
Status: frozen before formal R17 execution

## Contradiction found

The initial R17 gate required both:

1. exclusion of every source pocket with development/confirmation 4-mer
   containment at least `0.40`; and
2. at least 20% development coverage at the same `0.40` threshold.

Those conditions are logically incompatible whenever the protected target has
a KLIFS pocket. This is a design error, not an observed result. The formal
audit has not run, so the gate is corrected before any result exists.

## Correct target-support object

Leakage remains defined by contiguous 4-mer containment at least `0.40`.
Support is instead measured on the separately aligned KLIFS 85-position
functional-site coordinate:

```text
aligned_identity(p,q) =
    mean_r [pocket_p[r] == pocket_q[r]], r=1..85.
```

A development pocket is distantly supported when:

- its maximum allowed-source aligned identity is at least `0.25`; and
- every allowed source remains below the `0.40` 4-mer firewall.

This permits shared catalytic-site chemistry without admitting a homologous
or exact target. Report both identity and containment distributions.

## Corrected gate 4

At least 20 ChEMBL development targets must have an available KLIFS pocket,
and at least 20% of those targets must have an allowed-source maximum aligned
pocket identity of at least `0.25`.

All other R17 thresholds and prohibitions remain unchanged.


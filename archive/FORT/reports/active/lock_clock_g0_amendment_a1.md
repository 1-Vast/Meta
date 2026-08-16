# LOCK/CLOCK G0 amendment A1: typed artifact comparison

**Frozen after invalidating all pre-A1 outputs and before the accepted rerun.**

Bound implementation amendment SHA-256:
`9930f1021ca09ba1e9e5e0335adc26d58c25e1a85097849ed6fbf49709947fa9`.

The recomputation comparator classifies Python and NumPy scalar values into four disjoint kinds:
boolean, integer, floating-point, and other. Kinds must match. Booleans and integers compare exactly;
only float-to-float comparisons may use absolute and relative tolerance `1e-12`. Lists must retain
length and order, and dictionaries must retain the exact key set.

Regression tests must reject:

- boolean `True` versus integer `1`;
- integer `372` versus float `372.0000000001`;
- integer permutation entry `0` versus float `1e-12`.

This is integrity hardening only. All scientific definitions and stop rules remain unchanged.

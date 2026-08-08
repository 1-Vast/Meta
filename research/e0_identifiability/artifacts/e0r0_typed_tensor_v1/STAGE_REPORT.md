# E0R0 Typed Tensor Identifiability

Decision: `TYPED_TENSOR_LEARNED_HEADS_FAIL`.

| Arm | Parameters | Correct CI | Deranged CI | Correct-Ligand | Correct-Deranged | Gate |
|---|---:|---:|---:|---:|---:|---|
| Analytic tensor | 240 | 1.00000 | 0.51250 | +0.51546 | +0.48750 | True |
| Learned full tensor | 240 | 0.69013 | 0.63026 | +0.20559 | +0.05987 | False |
| Learned CP rank 6 | 114 | 0.66776 | 0.60066 | +0.18322 | +0.06711 | False |
| Frozen generic MAP | comparison | 0.68553 | 0.64934 | +0.20099 | +0.03618 | False |

The analytic 240D residual error is `4.51e-07`.
The numerical CP-rank-6 tensor/residual errors are
`2.34e-13` and
`4.51e-07`.

Both learned heads were still improving at the frozen endpoint. Full-tensor
loss changed from `0.97938` to `0.80043`; CP-rank-6 loss changed from `0.98218`
to `0.66852`. Their epoch-60 gradients remained nonzero, and their learned
tensor cosine against the analytic teacher was only `0.21025/0.22827`.

Therefore the typed tensor has an exact realization guarantee, but the frozen
E0 optimization protocol did not identify that realization within its
registered budget. This run does not authorize an epoch extension, optimizer
change or post-hoc model selection.

No affinity labels or DAVIS data were read. PLIP/typed-interaction training,
production integration, CSMO/Band changes and downstream authorization remain
outside this stage.

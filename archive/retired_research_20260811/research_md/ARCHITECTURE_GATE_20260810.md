# Architecture gate for the proposed Theory-Projected Q-PMA

Date: 2026-08-10

## Decision

```text
ARCHITECTURE_GATE_FAIL
MODEL_AND_SCRIPTS_REFACTOR_NOT_AUTHORIZED
REAL_EXPERIMENT_NOT_AUTHORIZED
```

The frozen linear Meta-Section remains the only admitted core. Its task state
is shared by every query and is restricted to the support row space:

```text
a_t = M_S^T (M_S M_S^T + lambda I)^-1 r_S, lambda > 0.
Delta_q = a_t^T m_q.
```

## Failure of the proposed objective

For the prompt's objective

```text
min_w ||M_S^T w - P_S m_q||^2
      + eta ||w - w_tilde_q||^2 + lambda ||w||^2,
```

`M_S P_S = M_S`, so the solution is

```text
w* = [M_S M_S^T + (eta + lambda) I]^-1
     [M_S m_q + eta w_tilde_q].
```

When `eta=0`, this is the existing ridge query leverage and adds no new
information. When `eta>0`, a component of `w_tilde_q` in
`null(M_S^T)` can affect `w*^T r_S`. Thus a query with zero row-space coverage
can receive a nonzero correction from support residual directions outside the
declared family. A separate `w*` per query also breaks the single shared task
coefficient contract.

## Status of a repaired attention arm

A safer construction first freezes the unique `a_t`, then forms a
feature-only, order-equivariant attention coordinate inside the row space:

```text
v_q = P_S m_q
m_att = ||v_q|| / (||M_S^T alpha_q|| + eps) * M_S^T alpha_q
g_q = coverage_q * sigmoid(h_q)
m_section = (1 - g_q) v_q + g_q m_att
Delta_q = a_t^T m_section
```

It must enforce `coverage_q=0 => g_q=0 => Delta_q=0`, and `alpha_q` may read
support/query features but no labels or residuals. This removes the null-space
failure, but changes the family from `m(P,L)` to `m_section(q; X_S)`. It is
therefore a future replacement hypothesis, not an implementation detail of the
frozen core. It requires a new preregistration and fresh evaluation targets.

## Other blocking contracts

- Numeric Ki availability leaves 29 evaluation targets at k=5, below the
  frozen minimum of 30.
- The legacy 288D bank covers only 12 of those evaluation tasks at k=5.
- The proposed 4--5 dimensional state has no registered map into the current
  bounded `[0,1]^28` CSMO interface.
- `project_state.json` keeps both training authorization and confirmation-label
  access closed.

No file is moved into `model/` and no runner is moved into `scripts/` until the
numeric, feature, support, biology and transfer Gates pass.

## Benchmark ownership

The failed Q-PMA does not imply that the project should invent an isolated data
ecosystem. The main benchmark may reference AdaMBind's protein-task and
protein-40 split and CARA's assay-aware numeric cleaning, while the existing
dependency-closed split becomes a strict confirmation lane. These are reference
conventions, not copied architecture. MetaSieve retains ownership of the Ki
estimand, biological representation, support-identifiable adaptation and law
output. CARA's assay-level aggregation and AdaMBind's protein-level task must be
reconciled explicitly before a new corpus is materialized.

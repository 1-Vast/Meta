# Conditional preregistration — X1B crossed interaction existence

Stage: `E-AFF-X1B_CROSSED_INTERACTION_EXISTENCE`

Status: frozen before any direct DD value is computed. Execution is conditional
on `X1A_R_DEPENDENCE_PRECONDITION_PASSED` for the endpoint.

For every exact-assay-aligned, label-blind capped rectangle:

```text
DD = y(P1,La) - y(P1,Lb) - y(P2,La) + y(P2,Lb)
D  = DD / 2
```

Estimate endpoint measurement variance from pooled exact-assay within-cell SSE
and degrees of freedom. Use its one-sided 95% chi-square upper bound
`sigma_rep,U^2`, and define conservatively:

```text
v_D,U = sigma_rep,U^2 / 4 * sum_c(1 / n_c)
Z     = D^2 - v_D,U
tau^2 = weighted mean(Z)
```

The factor of two places tau on the cell-interaction scale used by X0-B.

Inference units are dependency components. The point estimate weights component
means by their frozen capped rectangle counts. Use a one-sided 95% small-G
cluster-t lower bound. For Kd also report all 4096 Rademacher sign patterns and
leave-one-component-out estimates. Rectangles are never treated as IID.

Gates, separately by endpoint:

```text
B1 assay/noise/hash contract valid
B2 X1A-R dependence PASS, n_eff >= 245, largest weight <= 0.25
B3 LCB95(tau^2) > 0
B4 sqrt(max(tau^2,0)) / sigma_rep,U >= 0.5
```

Exactly one endpoint verdict:

```text
X1B_DATA_OR_NOISE_CONTRACT_INVALID
X1B_INTERACTION_NOT_DETECTED
X1B_INTERACTION_PRESENT_BELOW_DESIGN_MARGIN
X1B_REAL_CROSSED_INTERACTION_IDENTIFIED
```

Mean DD is descriptive only and cannot pass a Gate. PASS identifies a
panel-compatible observed crossed interaction beyond repeat noise, not physical
free energy or causal mechanism. PASS authorizes only an X2 preregistration;
X2 is not trained in this programme run.

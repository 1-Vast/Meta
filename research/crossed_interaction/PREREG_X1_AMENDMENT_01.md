# Amendment 01 — X1

Parent: `PREREG_X1_ICC_AND_DATA_CONTRACT.md`, SHA-256
`67e3c651de8d3f932934d76fb955f4554ded668551200ca439253d0548549bbc`,
committed `008c82a`.

Written 2026-08-10, immediately after the first X1A execution exposed that the
parent's ICC estimator is **degenerate by construction**.

## 1. The defect

Parent section 5 specifies that additive target and ligand effects be removed by
a two-way least-squares fit "fitted within panel so that panel-level offsets
cannot leak into `r`", and then decomposes `r` into cluster, panel, cell and
replicate components with

```text
rho = var(u_cluster) / (var(u) + var(v) + var(w) + var(e))
```

A per-panel intercept forces the least-squares residuals of every panel to sum
to zero. Every panel mean of `r` is therefore exactly zero, every cluster mean
of panel means is exactly zero, and `var(u_cluster) == 0` for **any** dataset.

This is not an empirical finding. It is an algebraic identity, demonstrated on
synthetic panels with deliberately injected cluster offsets of 10, 20 and 30
log units, all of which returned residual means of order `1e-15`.

The first X1A execution duly returned `rho = 0.0000` and
`ucb95 = 0.0000` for both endpoints and would have reported
`X1_ICC_PRECONDITION_PASSED`. **That result is void and is not used.** A Gate
that cannot fail is not a Gate, and reporting a pass obtained from one would
have been the most damaging possible outcome of this stage.

## 2. Correction

Parent section 5's first step is replaced:

```text
OLD  fit mu + a[target] + b[ligand] within each panel
NEW  fit mu + a[target] + b[ligand] once per endpoint, globally across that
     endpoint's cells, by sparse least squares
```

Removing additive target and ligand effects is retained, because that is
exactly what `DD` cancels. Fitting them globally rather than per panel leaves
panel-level and cluster-level structure in the residual, where it can be
estimated. Every remaining element of section 5 — the nested
cluster/panel/cell/replicate decomposition, the Henderson-III moment
estimators, truncation of negative components at zero, identification of
`var(e)` from exact-assay replicates, and the 10,000-draw cluster bootstrap at
seed `20260903` — is unchanged.

## 3. What is explicitly NOT changed

No threshold, margin, Gate, seed, endpoint definition, cluster definition,
firewall rule, verdict name or stopping rule is altered:

```text
G1 Ki  UCB95(rho) < 0.0915      unchanged
G2 Kd  UCB95(rho) < 0.0164      unchanged
G3 largest capped cluster share <= 0.25   unchanged
G4 recomputed effective n >= 245          unchanged
```

Ki and Kd remain completely separate. The X0-B design and its statistical unit
are still not rebuilt or replaced.

## 4. Reporting requirement

The corrected run must report **both** quantities side by side: the parent's
within-panel `rho`, labelled as structurally zero, and the corrected global-fit
`rho`. The void first execution is described in the evidence report rather than
deleted.

## 5. Honesty note

This amendment was written after affinity values had been read, which the
parent's own firewall permitted at that point. It changes an estimator that was
provably incapable of producing a non-zero value, in the direction that makes
the Gate able to fail. It does not move a threshold in response to an
unfavourable number, and the direction of the correction is against the
stage's own convenience: the void result was a pass.

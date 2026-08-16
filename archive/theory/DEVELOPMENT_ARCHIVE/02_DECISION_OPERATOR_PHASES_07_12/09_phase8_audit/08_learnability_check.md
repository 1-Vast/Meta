# Learnability Check

## Verdict

`INCOMPLETE_AND_PARTLY_INCORRECT`

## Required separation

### A. Existence of the decision object

Established under the compactness/lower-semicontinuity conditions used for the argmin theorem, or via an explicitly declared approximate action set.

### B. Identification of the decision object

Established set-wise from \(J_Q(O)\), the declared loss/criterion, population ambiguity class, and tie-break. Single-valued identification fails at undeclared ties by design.

### C. Estimation from finite historical tasks

DR-L3 gives an explicit forced/compatible frequency construction for pairwise signs and listwise orders. This is substantive, but the advertised theorem is not valid under its stated assumptions.

## Defects

1. **Exchangeability is not IID.** DR-L2 says joint exchangeability buys Hoeffding/DKW rates and that conditional exchangeability gives the same rates within a fiber. Mere exchangeability does not imply those concentration bounds. For example, let every task variable equal one shared Bernoulli \(Z\). The sequence is exchangeable, but its empirical frequency is always 0 or 1 and does not concentrate around the marginal probability \(1/2\). IID, conditional IID, or an explicit dependence/concentration assumption is required.
2. **Simultaneous listwise coverage is missing.** The displayed \(\eta_n=\sqrt{\log(4/\delta)/(2n)}\) supports a fixed order's two bounds. Intersecting intervals for all \(m!\) orders does not retain probability \(1-\delta\) without a simultaneous argument, such as allocating \(\delta/m!\) or proving a multinomial uniform bound.
3. **The “if and only if” is overbroad.** Per-task exact/set identification is sufficient for this distribution-free forced/compatible construction, but it is not necessary across all statistical models; a declared measurement-error or likelihood model can identify an aggregate law without identifying every historical latent value. The theorem must restrict its claimed necessity to the frozen distribution-free information model.

With conditional IID, corrected simultaneous confidence, a precisely scoped necessity statement, and a properly indexed pushforward target, the finite-history construction can become a valid conditional learnability theorem. As written, Phase 8 has not proved the claimed complete learnability result.

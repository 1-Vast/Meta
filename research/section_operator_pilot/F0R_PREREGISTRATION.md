# F0R preregistration: support design and section reliability

Frozen after F0 returned `F0_CEILING_NOT_VIABLE` and before F0R scoring.
PKIS2 and Anastassiadis2011 remain consumed development panels; KCGS numeric
outcomes remain unread.

## F0 failure localized

At `k=5`, correct adaptation beat wrong-target support and the protein controls
on PKIS2 but lost to support-free prediction. At `k=20`, it beat support-free
and location-only. This pattern identifies estimation variance/conditioning,
not absence of a section signal. F0 also failed to exclude query ligands sharing
a generic Murcko scaffold with support ligands; F0 is retained as a pilot but
cannot support a leakage-safe claim.

## One permitted revision

F0R may change only:

1. support construction: compare scaffold-disjoint random support with a
   label-blind greedy D-optimal support over `[1, phi_1, ..., phi_r]`;
2. select both location and interaction positive-ridge penalties on source-only
   dual-cold episodes;
3. fit a source-episodic reliability gate that mixes support-free and adapted
   predictions. Gate inputs are restricted to label-blind condition/coverage
   statistics plus support residual dispersion and fit residual. It cannot read
   target/ligand IDs or any query label.

For every support molecule, all query molecules with the same generic Murcko
scaffold are excluded. Source/transfer strict scaffold-cold exclusions remain.

## Theory relation

Support design addresses the frozen theory's result that identifiability is a
joint property of family and design. Reliability mixing is implemented as a
convex mixture of prior and adapted coefficient/law outputs, so simplex and
band validity are preserved. It is not an added scalar deployment head.

## Gate

The primary setting remains `k=5`; `k=20` is diagnostic. D-optimal F0R must:

- improve over support-free, location-only, wrong support and label permutation
  on PKIS2 with target-bootstrap 95% lower bounds above zero;
- have positive point estimates against the same arms on Anastassiadis2011;
- beat scaffold-disjoint random support on PKIS2;
- retain a positive correct-minus-deranged-protein point estimate;
- show higher median support minimum singular value than random support.

If this ceiling fails, no neural section operator is trained on this basis.
The only remaining live paths are a different biological pair basis or an
explicit increase in support budget.


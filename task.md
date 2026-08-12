# Active QPSMP Work Contract

## Question

Can one shared, query-loss-trained QPSMP neural meta-potential use correctly bound support labels to
improve protein-specific Cold Target affinity prediction beyond identical-budget additive, level,
ligand-only, analytic-ridge, wrong-protein, shuffled-protein, and design-nuisance controls?

## Allowed Inputs

- source-trained or legally frozen protein residue representations;
- ligand molecular graphs encoded by the shared ligand encoder;
- declared measurement context;
- `k={1,2,3,5}` support observations from the unseen recipient target.

Target IDs may index cached tensors but cannot enter the model. Query labels, target memory, and
recipient-specific trainable parameters are prohibited.

## Primary Module

`QPSMPBioModel` is the primary learned path. Query loss must deliver finite gradients to the protein
encoder, ligand encoder, localizer, crossed scalar head, section basis, and neural support adapter.
The analytic centered ridge is a comparator/diagnostic and cannot be the sole primary arm.

## Frozen Invariants

- support/query rows are disjoint;
- pair inclusion and orientation are outcome-independent;
- support ordering does not change output;
- zero centered support evidence implies zero SAR update;
- `k=1` reports level utility only and cannot pass SAR adaptation;
- delta and rectangle outputs are differences of the retained scalar endpoint path;
- additive, cross-zero-shot, level, and SAR-adaptation channels are reported separately;
- foreign support changes only the transient SAR state;
- validation episodes are fixed before checkpoint selection;
- component/dependency weighting and consumed-development status are explicit.

## Next Gate

The governed raw-biological-input trainer now runs. The next experiment must use a preregistered
budget across `k={1,2,3,5}`, report absolute utility and all matched controls, and reserve untouched
components for confirmation. The existing 20-step CPU result is only a trainability/implementation smoke.

## PASS

Proceed beyond development only if preregistered component-level lower bounds show useful full-scalar
gain and correct-protein crossed/SAR specificity, with no support-binding, target-main, document,
panel, or query leakage.

## STOP

Fail closed if the learned arm only beats zero, only separates an artificially destructive wrong
arm, loses to level/additive/ligand-only controls, uses unstable validation selection, or depends on
overlapping development units. Such failure closes the frozen recipe, not the entire function class.

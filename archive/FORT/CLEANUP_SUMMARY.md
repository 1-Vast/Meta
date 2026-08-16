# FORT Cleanup Summary

Date: 2026-08-05

## Decision

The FORT model-development task is stopped and the working tree has been
cleaned because the available datasets cannot support the required model claim.
No model training or further architecture work is authorized from the removed
materials.

## Consolidated Findings

- The natural-tail data gate is `DATA_NOT_READY`: its strict admitted
  recipient roster is empty, and the most permissive pKi upper bound contains
  40 recipients, below the frozen minimum of 50.
- The formal ChEMBL 37 pKi corpus supports source-side work but does not supply
  a sufficiently valid prospective natural-tail evaluation population.
- The strict unseen-target few-shot route has no established
  protein-conditioned, support-dependent ranking gain. The recorded pilot and
  control results did not exceed appropriate calibration or protein-free
  baselines.
- pKi and pKd are distinct endpoints and cannot be pooled to compensate for
  missing valid units.
- New data acquisition or a rebuilt, outcome-blind provenance-valid roster is
  required before this task can be restarted. Reusing the deleted data for
  further model experimentation would not resolve the data gate.

## Cleanup Scope

All current project content has been removed: datasets, generated artifacts,
source code, configurations, manifests, scripts, tests, temporary files, and
the original report collection. Git metadata is retained so prior repository
history remains auditable. This summary is the sole retained working-tree file.

# LEXOR L0 Decision

## Verdict

`LEXOR_L0_CORPUS_FRAME_INSUFFICIENT_STOP`

This verdict applies only to the frozen local inventory stated in the
preregistration. It does not claim that all open-access literature has been
searched. It authorizes no LLM/API call and no model training.

## Corpus Frame

| measure | value |
| --- | ---: |
| documents in inventory | 766 |
| provenance families | 762 |
| eligible family environments | 0 |
| required family environments | 30 |
| declared components | 766 |

## Gates

| gate | pass |
| --- | --- |
| >=30 eligible provenance families | False |
| >=40 declared scaffold-diverse queries per selected family | False |
| verified acquisition list | False |
| high-noise MDE80 design proxy <=0.03 | False |

## Power

Status: `not_evaluable`.

The proxy is deliberately not an empirical model-power result. L0 consumed no
measurement values, prediction values, or model outputs.


## Why This Local Frame Failed

The frozen inventory contains no eligible provenance-family environment. The
registered blockers are:

* license is absent, unrecognized, or unverified: 702 declared component(s)
* no explicit scaffold-diverse query-ligand count: 766 declared component(s)

This is a local-frame stop, not a claim that all open-access literature has
been searched. Expanding discovery requires a new inventory and a new L0
preregistration; it cannot reopen this frozen run.


## Firewall State

* raw measurement files read: `False`
* external network called: `False`
* LLM API called: `False`
* model trained: `False`
* sealed test consumed: `False`

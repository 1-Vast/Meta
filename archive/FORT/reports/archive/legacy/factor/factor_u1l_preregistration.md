# FACTOR-U U1-L preregistration: local strict-corpus representation proof

Date: 2026-07-26  
Route origin: user-supplied strict unlabeled-corpus expansion; no agent candidate slot consumed.  
Role: local mechanism proof before full-corpus training on the larger machine. F1-C remains locked.

## Authorization and hypothesis

U0-B passed with 138,805 globally firewalled molecules and 60,021 scaffolds. U1-L tests one claim:
whether more chemically diverse, completely evaluation-disjoint unlabeled graphs can repair the
F0-C1 participation-rank, atom-decoy-calibration and scaffold-OOD defects without changing the
encoder family or using activity labels.

## Frozen local proof corpus

Rebuild the exact U0-B retained union. For each nonempty Murcko scaffold (or each acyclic
connectivity as its own key), choose its lowest SHA-256 molecule. Rank scaffold keys by
`SHA256("factor-u1l-scaffold-v1:" + key)` and retain the first 50,000 keys. Thus the local proof has
exactly one molecule per scaffold and is not selected by any evaluation statistic.

Split these scaffolds 90/10 by
`SHA256("factor-u1l-validation-v1:" + scaffold_key)`. Train on the 90% block; use the 10% block only
for structural checkpoint selection. If malformed graphs reduce the corpus below 48,000 total or
4,000 validation molecules, stop as an engineering/data failure.

This 50,000-molecule balanced subset is a local proof, not the final pretraining corpus. A passed
result must be reproduced on all 138,805 molecules and with uncapped pair carriers on the larger
machine before F1-C can be considered.

## Frozen encoder

Use the exact F0-C1 architecture and objectives:

- four residual GINE layers, width 128, edge width 16, LayerNorm, dropout 0.10;
- the same six categorical atom inputs and three edge inputs;
- the same four deterministic connected 15% masks;
- masked element/charge/degree/aromaticity/hybridization, directed bond and BRICS-attachment heads;
- loss = mean atom CE + 0.5 bond CE + 0.25 attachment CE;
- AdamW 1e-3, weight decay 1e-5, batch 64, seed 1729.

Train the single globally disjoint encoder for 12 epochs. This gives more than seven times the graph
exposures of an F0-C1 fold while remaining a local proof. Checkpoint solely by minimum public
unlabeled validation structural loss. No F0-C1 checkpoint is reused.

## Frozen evaluation

After checkpoint freeze, encode the original three evaluation sources once. For each held source,
rebuild the exact F0-C1 strict outer connectivity/scaffold deletion and inner 80/20 scaffold split.
The atlas consists only of that fold's inner-train evaluation molecules; the public pretraining
corpus never becomes an atlas shortcut.

Use the exact F0-C1 atom/pair/motif carriers, 96-pair audit cap, robust scaling, per-role bandwidth,
role/chemistry/environment decoys, equal level weights, rarity weights, source balance and 10,000
source-stratified bootstrap draws.

## Gates

All must pass:

1. public-unlabeled validation masked element margin >=0.15, bond margin >=0.10 and attachment
   margin >=0.10; all losses/gradients finite;
2. in every evaluation fold, 5-NN atom-role macro-F1 margin >=0.15;
3. in every fold, atom participation rank >=16 and ratio >=0.10;
4. every validation decoy family changed fraction >=0.95;
5. every level's calibration false coverage <=0.05;
6. every inner scaffold-OOD median >=0.85 and q10 >=0.60;
7. source-balanced external median >=0.90 and q10 >=0.70;
8. true-minus-decoy grouped LCB95 >0, source weight <=0.40, inherited primitive graph and
   MDE80 pass;
9. all corpus firewalls remain zero; current-run confirmation labels unread; sealed test unconsumed.

Pass: `FACTOR_U1L_PASS_AUTHORIZE_FULL_CORPUS_REPLICATION`.  
If structural/rank/decoy/calibration gates fail:
`FACTOR_U1L_REPRESENTATION_UNIDENTIFIED_STOP`.  
If all representation gates pass but coverage fails:
`FACTOR_U1L_REAL_CHEMICAL_SUPPORT_FAIL`.

No after-result change to subset size, epochs, architecture, seed, bandwidth or thresholds is
allowed.

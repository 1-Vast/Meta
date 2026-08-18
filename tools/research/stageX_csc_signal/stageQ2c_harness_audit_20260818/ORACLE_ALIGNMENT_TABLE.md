# Oracle alignment table (Q2c-0 item 1)

Evidence sources: `stageX0c_measurement_qualification_20260818/Q2_PLANTED.json`
(frozen artifact) and the diag21 run log (`diag21.log`, retained as evidence,
protocol reconstructed from the committed diag21 run script's command line and
the recorded output).

| Dimension | diag21 P@U probe (cited in X0c REPORT) | oracle_protein arm (in Q2_PLANTED.json) |
|---|---|---|
| model class | Q2Model(97, n_prot) - MLP encoders + rank-4 interaction head | Q2Model(97, n_prot) - same class |
| protein input | prot_feats @ lat['U'] (float32) | prot_feats @ lat['U'] (float32) |
| planted truth | generate(1.0, 4, 'dense', False, seed=0) for ALL three runs | generate(1.0, 4, 'dense', False, seed=s) per seed s |
| train/eval graph | same load_duongly_graph splits | same |
| tau / rank / noise | 1.0 / 4 / noise_sd=1.0 | 1.0 / 4 / noise_sd=1.0 |
| endpoint | logit + interval censoring (same censored_loss) | same |
| training | single-phase AdamW, 6000 steps, final-state model evaluated | train_arm_with_restarts: 8 restarts x 6000 steps, init seeds 0-7, BEST-VAL-LOSS checkpoint selected |
| seeds | optimizer seeds 0,1,2 - three runs share ONE truth (seed 0) | generate seeds 0,1,2 - three independent truths |
| projection | none (raw head output) | none for interaction_head metrics (raw head output) |
| observed dz | 0.730 / 0.677 / 0.755 | 0.674 / 0.607 / 0.664 |
| observed Spearman | 0.448 / 0.433 / 0.478 | 0.375 / 0.331 / 0.390 |

Conclusion: the two protocols differ in truth seeds, restart/selection rule and
checkpoint choice. The X0c REPORT sentence citing dz 0.68-0.76 for the oracle
arm is WITHDRAWN; the only admissible oracle evidence is the in-artifact
oracle_protein arm: dz 0.607-0.674, Spearman 0.331-0.390 - below the frozen
0.70 dead-zone threshold in all three seeds. The correct statement is: a
frozen oracle probe shows recoverable information in part of the setting, but
the in-artifact oracle arm has not stably reached the full gate.

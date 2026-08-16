# Paper evidence audit (2026-07-26)

Scope: original papers named in the two user-supplied reading lists. The lists' proposed innovations
are deliberately excluded. This note records only the papers' assumptions, demonstrated results,
limitations, and relevance to the current strict dual-cold affinity problem.

## Direct evidence for the present DTA work

| Paper | What the paper actually demonstrates | Boundary for the present work |
|---|---|---|
| AdaMBind (Nature Communications, 2026) | MAML with a learned easy-to-hard task scheduler. Two BiLSTMs use query loss and support/query gradient similarity to score candidate tasks; the scheduler is optimized on validation tasks. The base predictor is a two-layer molecular GNN plus a three-layer protein 1D CNN. | This is evidence for adaptive episode allocation, not evidence that its encoder or scheduler solves strict scaffold-plus-homology dual cold affinity regression. Its scheduler also uses query labels during meta-training and is an auxiliary training policy, not a new core affinity mechanism. |
| Resolving data bias / PDBbind CleanSplit (Nature Machine Intelligence, 2025) | A filter combining protein TM-score, ligand Tanimoto similarity, pocket-aligned ligand RMSD, and affinity similarity removes train-test overlap and training redundancy. Search baselines and published scoring models deteriorate after filtering; protein removal becomes strongly harmful under CleanSplit. | Strong evidence that random or nominally disjoint splits can reward memorization. The filter itself uses test affinity proximity when deleting training complexes, so it is not a deployable label-blind split rule for our test protocol. The transferable result is the need for multiple label-blind structural/chemical firewalls and matched protein ablations. |
| Learning characteristics of GNN affinity predictors (Nature Machine Intelligence, 2023) | Edge-level explanations across six GNNs show strong ligand memorization, little protein memorization, and increasing use of interaction edges for some models at high affinities. The authors do not conclude that the models comprehensively learn physical interactions. | Supports keeping ligand-only, shuffled-protein, and random-protein controls. It does not validate a sequence-only model or demonstrate generalization under our dual-cold protocol. |
| Co-folding physics stress tests (Nature Communications, 2025) | AlphaFold3, RosettaFold All-Atom, Chai-1, and Boltz-1 often preserve ligand poses after physically disruptive protein or ligand edits and can output high-energy structures with high confidence. | Confidence, attention, or a plausible pose is not evidence that a model learned binding physics. Mechanistic claims require counterfactual perturbations with an expected physical direction. |
| PoseBench (Nature Machine Intelligence, 2025) | Evaluates blind docking without a supplied pocket on recent bioassemblies, including pocket identification, pose accuracy, chemical validity, residue specificity, cofactors, and multiligand clashes. | A useful model-evaluation pattern, but a different task from sequence/few-shot affinity regression. It motivates multi-axis validity checks rather than a single headline metric. |
| ColdstartCPI (Nature Communications, 2025) | Mol2Vec and ProtTrans token features are cross-conditioned with Transformers for binary CPI prediction. Its "blind start" randomly partitions compound and protein identities, then reports AUC/AUPR. | Not comparable to continuous affinity regression. Identity-disjoint random partitions do not enforce scaffold disjointness, target homology disjointness, document firewalls, or temporal independence. Performance and feature visualizations do not by themselves establish induced-fit physics. |
| EquiScore (Nature Machine Intelligence, 2024) | A heterogeneous equivariant complex graph integrates physical priors and is trained with augmented, redundancy-reduced structure data; it performs strongly on external virtual-screening sets. | Evidence for interaction-structured priors when 3D complexes exist. It neither supplies those structures for our sequence-only setting nor establishes few-shot target adaptation. |

## Meta-learning, task selection, and attribution

| Paper | What the paper actually demonstrates | Boundary for the present work |
|---|---|---|
| Deep Kernel Transfer (NeurIPS, 2020) | Meta-learns a neural kernel by maximizing per-task GP marginal likelihood over support and query observations. The GP analytically replaces the inner loop and provides predictive uncertainty. Regression evidence is mainly sinusoidal functions and head-pose trajectories. | Direct precedent for the current Bayesian-kernel family, so a learned kernel plus exact posterior is not a new innovation by itself. Its gains assume related tasks share a transferable covariance. |
| Hierarchical Bayesian Model for Few-Shot Meta Learning (ICLR, 2024) | Uses global and episode-specific weight variables with a Normal-Inverse-Wishart variational model. A diagonal quadratic approximation estimated with SGLD gives closed-form local/global coupling and reduces MAML memory. It improves accuracy/calibration on image classification and synthetic/vision regression. | Strong prior art for global/local Bayesian shrinkage. Its i.i.d. task-distribution, diagonal quadratic, and Gaussian-noise assumptions are not evidence that protein-conditioned covariance is identifiable in sparse DTA panels. |
| In-Context Learning Is Provably Bayesian Inference (2025 preprint / ICLR 2026 submission) | Under squared loss, decomposes risk exactly into model-dependent Bayes Gap and irreducible Posterior Variance. Bounds are derived for a uniform-attention Transformer under boundedness and Hölder-smooth Bayes predictors; task-mixture uncertainty can decay quickly with context. | A useful diagnostic decomposition, not a guarantee for a practical Transformer or DTA. Once Bayes Gap is small, more architecture cannot remove task noise/posterior variance; more informative context or better labels are required. |
| On Data Efficiency of Meta-learning (AISTATS, 2021) | Stability bounds and experiments show that performance depends separately on number of tasks and examples per task. It evaluates methods under a fixed total supervision budget and actively selects support points at meta-training time. | Direct support for reporting label budget and task count, rather than treating repeated episode resampling as new supervision. It does not justify a particular DTA architecture. |
| Near-Optimal Task Selection with Mutual Information (AISTATS, 2022) | Selects fixed candidate tasks using mutual information between latent task vectors in an implicit-process Bayesian model. Greedy selection has a \(1-1/e\)-type guarantee up to a meta-parameter entropy constant; inference uses SVGD particles. | The guarantee depends on the specified latent-task graphical model, conditional independence, known sample sets for every candidate task, and fixed candidate pool. A posterior-variance heuristic in another model does not inherit this guarantee. |
| DETAIL task-demonstration attribution (NeurIPS, 2024) | Treats hidden states as features of an internal kernel ridge regressor and applies influence functions to rank demonstrations. It helps curation/reordering, especially with corrupted demonstrations, and can transfer rankings from a white-box to a black-box LLM. | This is an attribution/selection method, not a learned predictor. Its kernel-ridge surrogate and internal-optimizer interpretation are assumptions; gains without corrupted demonstrations can be small or inconsistent. |
| Motion Attribution / Motive (ICML 2026) | Uses motion-masked diffusion loss gradients, common randomness, an identity-Hessian approximation, and projected gradient similarity to rank fine-tuning clips specifically by motion influence. | The main transferable evidence is methodological: attribution must isolate the behavior being claimed and validate rankings by intervention. Gradient similarity is an approximation, not proof of causality. |
| LUPI (JMLR, 2015) | Formalizes training-only privileged information and proposes similarity-control (SVM+) and teacher-student transfer mechanisms. | Establishes the paradigm, not a universal guarantee that privileged information transfers to a student. |
| Rethinking Knowledge Transfer in LUPI (2024 preprint) | Reproductions show several reported gains disappear with adequate training or sample size; replacing PI with constants can match TRAM, and four real-world datasets show no improvement over no-PI. | Any teacher/posterior auxiliary channel needs matched capacity, constant/random-information, training-time, and no-PI controls. Otherwise apparent transfer is not identifiable. |

## Broader papers read; indirect rather than architectural evidence

- The Flexibility Trap: arbitrary-order diffusion decoding can defer high-entropy logical forks,
  reducing Pass@k solution coverage despite competitive Pass@1. The result is specific to diffusion
  language-model reasoning and does not establish a DTA mechanism.
- AI-redesigned starting points and outcomes enhance protein evolution: ProteinMPNN/PROSS
  redesign changes the experimentally accessible fitness landscape before PACE. This is wet-lab
  evidence that starting-state quality can change later adaptation, not evidence for a neural
  affinity-prediction block.
- High-accuracy sampling for diffusion models and log-concave distributions: a primarily theoretical
  first-order rejection/proximal sampling framework obtains polylogarithmic dependence on target
  accuracy under explicit score-error and distributional assumptions; experimental evaluation is left
  for future work. It is not currently actionable for deterministic DTA regression.
- Developmental visual diet: chronological preprocessing that models acuity, colour, and especially
  contrast-sensitivity development increases shape bias and several robustness measures across vision
  models, with an accuracy/shape-bias trade-off. This is curriculum evidence in vision, not a DTA
  architecture result.
- Evo 2: a 7B/40B StripedHyena-2 genomic language model trained on trillions of bases supports
  million-base context and strong zero-shot genomic analyses. Protein DMS performance is competitive
  with general protein language models but below specialized state of the art and can saturate with
  scale. It supports careful pretrained-feature benchmarking, not guaranteed DTA gains.
- The Obfuscation Atlas: separates honest policies, blatant deception, obfuscated policies, and
  obfuscated activations. It shows detector penalties and representation drift have distinct causal
  pathways. For this project it reinforces the need for a failure taxonomy and intervention-based
  controls; it supplies no affinity predictor.
- Toward universal steering and monitoring: Recursive Feature Machines extract supervised concept
  directions from internal activations using kernel ridge regression and average gradient outer
  products. Steering and monitoring work across many concepts/models, but this is again kernel/feature
  prior art rather than a task-adaptive affinity estimator.
- RL-guided crystal generation: group-relative policy optimization with explicit creativity,
  stability, composition-diversity, and structure-diversity rewards moves a latent diffusion model
  along a novelty-validity frontier. Reported validity relies on computational reward models, not wet-lab
  synthesis in this paper.
- Octopus-inspired peripheral control: embedded suction-cup sensors and local controllers reduce
  communication and central computation in a physical soft robot. It demonstrates decentralized
  control under local sensing; it is not evidence that independently processed ligand/protein channels
  will improve DTA.

## Targeted set-to-function, matrix, and ranking papers supplied later

These notes use the papers themselves and exclude the architecture suggestions
in the supplied summary.

- **Neural Operator Processes (2026 preprint).** The paper combines
  neural-process conditioning with a fixed neural-operator decoder and compares
  SetConv pooling with query-aligned attention on GP regression and three PDEs.
  Preserving local context-query geometry helps the non-periodic Darcy problem,
  but attention is not uniformly best; benefits depend on PDE geometry and the
  probabilistic objective. Context budgets are 16–256 points, not four. A global
  latent can overwrite useful local geometry or absorb error into uncertainty,
  and the authors identify fixed-FNO geometry bias and training cost as
  limitations. This supports retaining query-local support geometry, but not a
  free operator decoder or a DTA performance claim.
  Primary source: https://arxiv.org/html/2606.22946v1
- **One Operator for Many Densities (2026 preprint).** The authors prove
  continuity and neural-operator approximation for a map from a joint density to
  its conditional over suitable compact density classes. Experiments are
  discretized correlated Gaussians and three-component Gaussian mixtures. The
  paper explicitly notes that density access and grids limit high-dimensional
  scalability. Four empirical support points do not provide the density input
  assumed by this theory.
  Primary source: https://arxiv.org/pdf/2605.06873
- **ANOVA Tensor Product Neural Networks (ICML 2025).** Tensor-product neural
  bases enforce component-wise sum-to-zero conditions and hence a unique,
  stable functional-ANOVA decomposition. The experiments concern interpretable
  tabular decomposition, not few-shot adaptation. The transferable result is
  structural identifiability: an interaction should not absorb an arbitrary
  main-effect offset. It supports the cross-fitted base and centered-residual
  contract but supplies no support-to-function inference rule.
  Primary source: https://proceedings.mlr.press/v267/park25d.html
- **R-learner / Quasi-oracle heterogeneous treatment effects (Biometrika
  2021).** The R-learner cross-fits conditional-outcome and propensity nuisances
  and then fits a residualized conditional treatment effect. Its quasi-oracle
  result is under a causal treatment model and kernel/regularity assumptions.
  DTA has no binary treatment assignment or matching propensity, so the theorem
  cannot transfer merely by renaming target-specific SAR as a treatment effect.
  Cross-fitted ligand nuisance estimation remains a sound analogy.
  Primary source: https://arxiv.org/pdf/1712.04912
- **SoftSort, NeuralNDCG, and PiRank.** SoftSort is a row-stochastic,
  temperature-controlled argsort relaxation; NeuralNDCG and PiRank use
  differentiable sorting to align learning with listwise ranking metrics.
  These are loss/operator components, not affinity predictors. PiRank's reported
  scaling benefit concerns large lists, whereas the Reinecke median is five
  query ligands per target. Introducing another temperature-sensitive listwise
  objective would require a separate ablation.
  Primary sources:
  https://proceedings.mlr.press/v119/prillo20a/prillo20a.pdf,
  https://arxiv.org/abs/2102.07831, and
  https://proceedings.neurips.cc/paper/2021/hash/b5200c6107fc3d41d19a2b66835c3974-Abstract.html
- **Sample-efficient inductive matrix completion with noise and inexact side
  information (2026 preprint).** Reduced sample complexity is proved for noisy
  low-rank IMC under effective side-information dimension, incoherence,
  conditioning, and observation assumptions; side-information misspecification
  is controlled through principal angles. Experiments are simulations and
  MovieLens. A scaffold/homology-blocked split is deliberately non-random, so
  this theory does not certify Reinecke, though it motivates a matrix baseline.
  Primary source: https://ar5iv.labs.arxiv.org/html/2605.17189v2
- **Orthogonal Inductive Matrix Completion.** OMIC jointly fits nuclear-norm
  components living in mutually orthogonal, prior-chosen side-information
  subspaces. It uniquely separates biases, communities, side-information
  effects, and residual low rank. It is not an unseen-target support encoder and
  uses static matrix directions/identities; its relevant lesson is again
  identifiable separation of main effects and interactions.
  Primary source: https://arxiv.org/pdf/2004.01653
- **Task-adaptive Neural Process for User Cold-Start Recommendation.** TaNP
  amortizes a latent task posterior, learns soft task clusters against a global
  pool, and modulates an adaptive decoder. This is close support-to-function
  prior art, but its task pool and generated decoder provide global paths that
  can bypass support-label information. Cold users with catalog items are also
  not simultaneous new-target/new-scaffold cases.
  Primary source: https://arxiv.org/pdf/2103.06137
- **Flow-based Adaptive Neural Process.** FANP adds a conditional normalizing
  flow and a modulation-augmented hypernetwork to TaNP. It shows greater latent
  flexibility on cold-start recommendation, not identifiability of a multimodal
  task posterior from four affinity observations. The extra flow/hypernetwork
  capacity would make the no-bypass contract harder to guarantee.
  Primary source:
  https://njuhugn.github.io/paper-conet/Towards%20Flexible%20and%20Adaptive%20Neural%20Process%20for%20Cold-Start%20Recommendation-Liu-tkde23.pdf
- **Transformers learn preconditioned gradient descent (NeurIPS 2023).** The
  result concerns linear self-attention without softmax on Gaussian linear
  regression. A one-layer global optimum implements one preconditioned gradient
  step; under sparsity assumptions, selected multilayer critical points
  correspond to iterative gradient algorithms. This is a mechanistic existence
  result, not a guarantee for nonlinear molecular regression with four noisy
  supports. Replacing an exact solve by an unrolled Transformer is not by itself
  an accuracy or novelty result.
  Primary source:
  https://proceedings.neurips.cc/paper_files/paper/2023/file/8ed3d610ea4b68e7afb30ea7d01422c6-Paper-Conference.pdf

## Evidence-level conclusion

The papers collectively strengthen the current failure localization rather than overturn it:

1. strict, label-blind chemical and target firewalls are more important than nominal "cold" identity
   splits;
2. any claim of target adaptation must survive label permutation, wrong-support, shuffled/random
   protein, ligand-only, and matched-capacity controls;
3. Bayesian task adaptation is well-established prior art and only helps when the shared task prior is
   identifiable;
4. adaptive task selection changes the supervision allocation, not the information content of a fixed
   rounded panel;
5. an architecture result is not credible if its apparent gain can be reproduced by longer training,
   extra capacity, privileged-channel constants, memorization, or a support-independent bypass.

No architecture is selected in this note. That decision remains contingent on the leakage and power
audit of the new Reinecke development panel.

## Follow-up: deep-kernel few-shot regression boundary

Three additional primary papers were read after SCGD/QACO failed:

- Patacchiola et al., NeurIPS 2020, learn a task-shared deep GP kernel by
  marginal likelihood (Deep Kernel Transfer):
  https://proceedings.neurips.cc/paper/2020/hash/b9cfe8b6042cf759dc4c0cccb27a6737-Abstract.html
- Chen et al., ICLR 2023, use bilevel implicit differentiation to adapt
  deep-kernel GPs for few-shot molecular property prediction:
  https://openreview.net/forum?id=KXRSh0sdVTP
- Falk et al., UAI 2022, learn a random-feature kernel distribution from
  regression tasks:
  https://proceedings.mlr.press/v180/falk22a.html

These works make a learned GP kernel, episodic marginal likelihood and
task-specific kernel fitting established prior art. More importantly, all
presume that task functions share a transferable covariance structure in the
chosen input representation. The Reinecke SI0 gate tested that premise before
implementing another kernel: Hodge-residual covariance alignment with Morgan
similarity was only +0.0271, below the frozen +0.03 minimum, and k=4
query-to-support maximum Tanimoto had median 0.185. Hence the papers justify the
form of a mature baseline but do not override the dataset-specific
identifiability failure; the planned learned kernel was correctly gated off.

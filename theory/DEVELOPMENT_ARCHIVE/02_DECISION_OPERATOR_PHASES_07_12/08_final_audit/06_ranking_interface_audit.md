# Ranking Interface Audit

## Marginals are insufficient

For two queries `Q={x_a,x_b}`, ranking depends on

$$
\delta=f(x_a)-f(x_b).
$$

The required identified object is the joint pushforward

$$
J_Q(O)=\{(f(x_a),f(x_b)):f\in I(O)\},
$$

or its decision-sufficient image

$$
\Delta_{ab}(O)=\{v_a-v_b:(v_a,v_b)\in J_Q(O)\}.
$$

Two families can have identical scalar marginals `[0,1]` at both queries:

$$
J_{\rm diag}=\{(t,t):t\in[0,1]\},
$$

$$
J_{\rm anti}=\{(t,1-t):t\in[0,1]\}.
$$

For the diagonal, `Delta_ab={0}` and the tie is identified. For the
anti-diagonal, `Delta_ab=[-1,1]` and both strict rankings remain possible. Scalar
marginal intervals cannot distinguish these cases.

## Population ranking object

When both signs remain admissible, pairwise 0-1 ranking needs

$$
p_{ab}=P(\delta>0\mid O)
$$

or an ambiguity interval for it. Graded ranking losses require the conditional
law of `delta`. Historical estimation additionally requires the cross-task,
likelihood, and historical-coverage assumptions stated in the shift audit.

The Phase-7 joint-object conclusion is correct. Its general output must use a
ranking-loss-specific robust floor rather than compare ranking loss numerically
to a value half-diameter.

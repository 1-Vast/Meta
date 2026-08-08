# Information Ceiling Audit

## Verified distinction

Historical population information may change a decision under ambiguity, but it
does not identify the hidden current member beyond current observations.

## Counterexample

Let two functions `f_-` and `f_+` agree exactly at every current support point,
but at the query satisfy

$$
f_-(x)=-1,\qquad f_+(x)=1.
$$

The current observation therefore identifies only

$$
I(O)=\{f_-,f_+\},\qquad J_x(O)=\{-1,1\}.
$$

Suppose historical tasks are IID from a declared population with

$$
P(f_+)=0.9,\qquad P(f_-)=0.1,
$$

and the history is sufficiently informative to estimate those frequencies. A
Bayes decision under squared error moves from the minimax center `0` to the
population mean `0.8`; a sign decision prefers `+`.

Nevertheless both current members remain observationally compatible. If the
true current member is `f_-`, the Bayes preference is wrong. The history changed
the action and its conditional risk, not the identified set.

Any report replacing `{-1,1}` by `{1}` or claiming an unconditional error below
the frozen minimax value fabricates current-task information.

# F4O preregistration: orthogonal bounded gauge-fixed section

Date frozen: 2026-08-08, after F4G localized the remaining transfer error to
location/amplitude coupling and before evaluating this orthogonal solve.
KCGS numeric outcomes remain unread.

## Change from F4G

For support curve values `s_i`, define `x_i=s_i-mean(s)`.  Replace the coupled
two-column ridge solve by the exactly orthogonal posterior

\[
 m=\frac{\sum_i r_i+\lambda_m\bar s}
          {k+\lambda_m},\qquad
 a=\Pi_{[0,2]}\frac{\sum_i x_i(r_i-\bar r)+\lambda_a}
                         {\sum_i x_i^2+\lambda_a},
\]

and predict

\[
    \hat r(P,L)=m+a\{s(P,L)-\bar s\}.
\]

Here `(m,a)=(mean(s),1)` is the support-free prior.  The centring makes the
location and shape columns orthogonal, so a wrong protein cannot improve its
location estimate merely by shifting its support-surface mean.  `m` is clipped
to `[-0.5,0.5]`; hence the two task coordinates are bounded.

All data sources, atlas construction, cold folds, candidate penalties, support
and query rules, controls, centred interaction endpoint, seeds, and gate are
identical to F4G.  Only this algebraic posterior is changed.  No development
outcome selects a penalty or bound.

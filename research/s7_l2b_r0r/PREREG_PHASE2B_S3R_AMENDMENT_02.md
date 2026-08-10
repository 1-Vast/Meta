# S3R preregistration amendment 02 - held-out B opening rule

Stage: `P1R2B-PHASE2B-S3R_REAL_STRUCTURAL_DIRECT_W`

Written before S3R implementation is committed and before any S3R score.

Held-out B is a scaffold-strict secondary development panel. It is not opened
unless module participation and all primary R1-R5 Gates pass. If the primary
stage fails, held-out B remains `NOT_OPENED_PRIMARY_DID_NOT_PASS`.

After a primary PASS, report held-out-B candidate AP, exact chance AP and their
difference. A negative difference is a sign reversal. It does not rewrite the
primary held-out-A verdict, but it blocks authorization of the I0 integration
stage and independent-confirmation escalation. It cannot be used to tune W,
change thresholds or rescue any primary Gate.

This rule replaces the parent document's ambiguous phrase that held-out B is
secondary but "must not" reverse sign.

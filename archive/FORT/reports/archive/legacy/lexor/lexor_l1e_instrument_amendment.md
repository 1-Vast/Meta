# LEXOR L1E Instrument Amendment

Date: 2026-07-27

## Purpose

`L1E` is a one-call provider connectivity and refusal-contract check authorized
directly by the user after the L0B11 metadata stop. It is separate from L1.
It does not claim extraction reliability, corpus feasibility, or affinity
utility, and it cannot authorize L2-L6 or F1-F4.

## Why It Does Not Weaken L0

L0B11 failed because public metadata cannot expose post-firewall query depth,
not because the provider was shown to be unnecessary or reliable. L1E therefore
uses no paper, supplement, raw measurement, structure, FORT record, or label.
Its single synthetic fixture contains no admissible measurement and requires an
abstention-only JSON response. A passing result proves only that the configured
endpoint can receive a bounded request and return a parseable refusal.

## Fixed Boundary

* one live request maximum; no retry;
* bounded output-token request and total-token budget;
* synthetic fixture only, with immutable prompt/fixture/firewall hashes;
* no numeric field is representable in the accepted response contract;
* no raw provider response, API key, endpoint value, or source data is written
  to an artifact;
* a result artifact records a response digest, usage, and pass/fail only;
* `StageAuthorization` remains unchanged and still refuses real L1 absent an
  L0 pass, blind fixture suite, and frozen firewall chain.

## Interpretation

`LEXOR_L1E_CONNECTIVITY_AND_REFUSAL_CONTRACT_PASS` means only that L1E may be
used as an engineering reference. It does not change the standing state
`L1_LOCKED`. Any failure is recorded as an instrument/configuration failure and
does not alter the L0B11 scientific conclusion.

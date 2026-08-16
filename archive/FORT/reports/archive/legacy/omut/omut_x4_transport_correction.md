# OpenMut `OMUT-X4` transport implementation correction

**Date:** 2026-07-28.
**Scope:** implementation correction under the frozen X4 transport rule; no
evidence threshold changed.

The first X4 execution classified any HTTP-200 parseable XML response as a
document. Both Elsevier text-mining URLs returned a
`full-text-retrieval-response` envelope containing only `coredata` metadata,
zero paragraph/section blocks, and no article body. This violated the frozen
requirement that the payload parse as a document.

The initial machine result reported 11 complete transports, two projected
construct fragments, zero accepted candidate/documents, and the already
negative four-component / three-family topology. It is superseded because
the transport count was wrong, although its topology verdict was not made
more favorable by the bug.

Correction:

- projection cache schema incremented so the old projections cannot be
  reused;
- a parsed response with zero document blocks is `invalid_document`;
- an offline regression test binds this behavior;
- X4 is rerun under the unchanged preregistration and evidence rules.

No source URL, construct rule, topology threshold, or outcome firewall was
changed.

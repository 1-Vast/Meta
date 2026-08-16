# DCST-R9 frozen-teacher atlas decision

Date: 2026-07-28  
Decision: `STOP_ATLAS_FAMILY_ADVANCE_TO_CONTINUOUS_CONTENT_ENERGY`

## Result

R9 reproduced the accepted R6 segment mechanism exactly:

- true centered alignment `0.05441`;
- target-destroyed `-0.01452`;
- ligand-destroyed `-0.05263`;
- target and ligand margins `0.06893` and `0.10704`.

The state audit copied 30 representation tensors from the R6 checkpoint,
copied no theta, and preserved the atlas centroids. The R6 checkpoint SHA-256
was `b29050444e69d18d9fe2879e229dc0253e03741f2b0d0b85d5b7ee93586f15ba`.

FTA privileged and FTA-NoPriv each certified `1/4` bands, so the frozen
privileged-attribution gate failed. Both activated zero-based band 3. The
privileged band confidence was only `0.10468`, whereas the no-privileged
confidence was `0.89957`. The atlas direction is therefore more consistent
with affinity-only structure than with privileged information.

Wall time was `166.037 s`; peak allocated CUDA memory was `533.6 MiB`. No new
downstream affinity label was loaded.

## Decision

The atlas family is stopped:

- learned soft roles were unidentifiable (R7);
- fixed hard roles preserved one privileged direction but damaged the segment
  mechanism under joint training (R8);
- freezing the mechanism restored that mechanism but the same atlas band was
  stronger without privileged labels (R9).

Discretizing diverse protein segments into eight global roles is the common
limitation. The next interface must use continuous ESM segment content
directly, without absolute position, learned role identity, or clustering.


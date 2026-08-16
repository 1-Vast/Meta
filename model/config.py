"""MetaSieve-DTA deployment configuration.

FROZEN MATHEMATICAL OBJECTS (must not change):
    V=[a_min,a_max], fixed mesh h, Delta_m, positive ridge mu, the band polytope,
    the assembly rule B(z), the operator A(F,z)=K(B(z)F(z)).

ENGINEERING CHOICES (declared, changeable, carry no theorem):
    d_z and the meaning of its blocks, the CSMO view bank, gate width/depth,
    anchor lattice, context quantiser, lambda_w, encoder widths.

The statistic dimension changed from the pre-biological realization (14 -> 28)
because the frontend now carries protein, ligand-graph and pair information that
did not previously exist. The frozen contract constrains z to be bounded,
finite-dimensional, deterministic, permutation invariant in S and free of query
labels -- it does not fix its dimension.
"""

from dataclasses import dataclass, asdict
import hashlib
import json


@dataclass(frozen=True)
class MetaSieveConfig:
    # ---- frozen value space and output grid (S-GRID) ---------------------
    a_min: float = 0.0
    a_max: float = 1.0
    M: int = 32                      # mesh intervals; grid has M+1 points

    # ---- frozen coefficient simplex and strong-convexity regularizer -----
    m: int = 7                       # anchors; Delta_m lives in R^{m+1}
    mu: float = 0.05                 # law-operator strong convexity; not few-shot ridge
    lambda_w: float = 0.35           # width weight inside L

    # ---- abstract statistic interface ------------------------------------
    # No biological construction is currently admitted to this interface.
    d_z: int = 28
    k0: float = 8.0                   # support-size scale constant
    s_max: float = 0.5                # max s.d. of a [0,1] label
    theta_dep: float = 0.5            # deployment id coordinate

    # ---- context map (finite codomain, frozen requirement) ---------------
    kappa_edges_y: tuple = (1.0 / 3.0, 2.0 / 3.0)   # 3 bins on z13
    kappa_edges_mass: tuple = (0.5,)                # 2 bins on z17
    # Transitional pre-P4 context quantizer. It has one reachable bin, so
    # continuous assay covariates cannot alter B before context_id -> B_c.
    kappa_edges_context_cont: tuple = ()

    # ---- population band estimator ---------------------------------------
    dkw_alpha: float = 0.10
    dkw_eps_min: float = 0.02
    dkw_n_min: int = 30
    dkw_eps_fallback: float = 0.25

    # ---- CSMO (mechanism frozen; view bank is an engineering choice) -----
    n_views: int = 6
    view_res: int = 6
    gate_width: int = 128
    gate_depth: int = 4
    gate_floor: float = 0.05

    seed: int = 20260804

    # ---- derived ---------------------------------------------------------
    @property
    def h(self) -> float:
        return (self.a_max - self.a_min) / self.M

    @property
    def D_V_val(self) -> float:
        return self.a_max - self.a_min

    @property
    def n_grid(self) -> int:
        return self.M + 1

    @property
    def band_dim(self) -> int:
        return 2 * self.n_grid

    @property
    def n_coef(self) -> int:
        return self.m + 1

    @property
    def n_context(self) -> int:
        return ((len(self.kappa_edges_y) + 1)
                * (len(self.kappa_edges_mass) + 1)
                * (len(self.kappa_edges_context_cont) + 1))

    @property
    def diam_simplex(self) -> float:
        return 2.0 ** 0.5

    @property
    def L_bar(self) -> float:
        return 1.0 + self.lambda_w

    @property
    def L_lip(self) -> float:
        return 1.0 + 2.0 * self.lambda_w

    def grid(self):
        import numpy as np
        return np.linspace(self.a_min, self.a_max, self.n_grid)

    def fingerprint(self) -> str:
        blob = json.dumps(asdict(self), sort_keys=True, default=str)
        return hashlib.sha256(blob.encode()).hexdigest()[:16]


DEFAULT = MetaSieveConfig()
assert DEFAULT.mu > 0.0, "The strong-convexity regularizer must be positive."


# ----------------------------------------------------------------------
# Scale profiles.
#
# The local workstation is used ONLY for construction and performance testing;
# real experiments run on separate hardware. So capacity is a declared knob, not
# a consequence of this machine. Every profile keeps the FROZEN mathematics
# identical (V, h, Delta_m, mu, the band polytope, the operator); they differ
# only in representation capacity, which is an engineering choice throughout.
# ----------------------------------------------------------------------
from dataclasses import replace as _replace     # noqa: E402

PROFILES = {
    "local": dict(gate_width=64, gate_depth=3),
    "base": dict(),
    "full": dict(gate_width=256, gate_depth=5, view_res=8),
}


def profile(name: str = "base", **overrides) -> MetaSieveConfig:
    """Return a config for a named scale profile, with optional overrides."""
    if name not in PROFILES:
        raise KeyError(f"unknown profile {name!r}; choose from {sorted(PROFILES)}")
    return _replace(DEFAULT, **{**PROFILES[name], **overrides})


# ----------------------------------------------------------------------
# z interface. Coordinates remain abstract until a biological statistic passes.
# ----------------------------------------------------------------------
Z_NAMES = [f"z{index}" for index in range(1, DEFAULT.d_z + 1)]
Z_INDEX = {name: index for index, name in enumerate(Z_NAMES)}

# CSMO view bank over abstract statistic coordinates. A future biological
# statistic must register its coordinate-to-view map before production use.
DEFAULT_VIEWS = [
    (8, 9), (10, 11), (0, 2), (17, 19), (12, 13), (4, 22),
]

assert len(Z_NAMES) == DEFAULT.d_z
assert all(0 <= a < DEFAULT.d_z and 0 <= b < DEFAULT.d_z for a, b in DEFAULT_VIEWS)

"""Generic geometry helper functions, not tied to a specific inflow/outflow
model component.

Translated from MATLAB source: GV_for_claude/utilities/geometry/
"""

import numpy as np
from numpy.typing import ArrayLike


def ellipse_perimeter(a: float | np.ndarray, c: float | np.ndarray) -> float | np.ndarray:
    """Ramanujan's approximation of the perimeter of an ellipse with semi-axes a, c."""
    # Shared by get_elliptical_cylinder_sa and s_spheroid, which each inlined it in MATLAB
    h = (a - c) ** 2 / (a + c) ** 2
    return np.pi * (a + c) * (1 + 3 * h / (10 + np.sqrt(4 - 3 * h)))


def get_elliptical_cylinder_sa(a: float | np.ndarray, c: float | np.ndarray, l: float =1, halved: bool=True) -> float | np.ndarray:
    """Ellipse perimeter approximation of Ramanujan.

    a, c = axes lengths. a is nominally semi-major but not critical
    l    = cylinder length. Default 1 to get perimeter only
    halved = bool. True [default] = half cylinder, total area/2
    """
    sa = ellipse_perimeter(a, c) * l

    if halved:
        sa = sa / 2

    return sa


# Translated from utilities/geometry/S_spheroid.m
def s_spheroid(a: ArrayLike, c: ArrayLike | None = None, n: int = 101, halved: bool = True, return_xz: bool = False) -> tuple:
    """Spheroid area and surface
    c = axis of symmetry, a = equatorial
    n = optional to set length of positional x and z vectors
    halved = bool. True [default] = hemispheroid, surface area/2

    OUT:
      S = surface area
      P = centered perimeter (approximation)
      x = horizontal position vector(s)
      z = height locations matching x

    Returns (S, P), or (S, P, x, z) when return_xz=True. x and z have shape (n, n_vecs).
    """
    # TRANSLATION NOTE: MATLAB's nargout-dependent outputs become (S, P) always, plus x, z when
    # return_xz=True (MATLAB's nargout >= 3). `halved` is new (MATLAB returned the full
    # spheroid area, except its c == 0 branch, which returned the hemispheroid's pi*a^2).
    if c is None:
        c = a
    a = np.asarray(a, dtype=float)
    c = np.asarray(c, dtype=float)

    assert a.size == c.size or a.size == 1 or c.size == 1, "a,c vectors must be same length, or one must be of length 1."
    a, c = np.broadcast_arrays(a, c)

    # TRANSLATION NOTE: MATLAB's min/max([a c],[],2) is row-wise only for column vectors (a row
    # vector collapses it to one global min/max); element-wise is the intended reading.
    def eccentricity(a: np.ndarray, c: np.ndarray) -> np.ndarray:
        return np.sqrt(1 - np.minimum(a, c) ** 2 / np.maximum(a, c) ** 2)

    def S_prolate(a: np.ndarray, c: np.ndarray) -> np.ndarray:
        return 2 * np.pi * a**2 * (1 + c / (a * eccentricity(a, c)) * np.arcsin(eccentricity(a, c)))

    def S_oblate(a: np.ndarray, c: np.ndarray) -> np.ndarray:
        return 2 * np.pi * a**2 + (np.pi * c**2 / eccentricity(a, c)) * np.log((1 + eccentricity(a, c)) / (1 - eccentricity(a, c)))

    # TRANSLATION NOTE: deliberate corrections - MATLAB's `if` on a vector only branched when ALL
    # elements matched (mixed vectors silently used the wrong formula), and a == c gave NaN
    # (zero eccentricity). Cases are now chosen element-wise, with a sphere case added. The
    # disc (c == 0) case is the full-spheroid limit 2*pi*a^2, so halved gives MATLAB's pi*a^2.
    S = np.empty(a.shape)
    disc = c == 0
    sphere = (a == c) & ~disc
    oblate = (a > c) & ~disc
    prolate = a < c
    S[disc] = 2 * np.pi * a[disc] ** 2
    S[sphere] = 4 * np.pi * a[sphere] ** 2
    S[oblate] = S_oblate(a[oblate], c[oblate])
    S[prolate] = S_prolate(a[prolate], c[prolate])
    if halved:
        S = S / 2
    S = S[()]  # 0-d array -> scalar for scalar inputs

    P = ellipse_perimeter(a, c)[()]

    if not return_xz:
        return S, P

    # Get roof heights
    # TRANSLATION NOTE: deliberate correction - MATLAB left a_vec/c_vec unassigned (an error) when
    # a and c were both scalar; scalars are now treated as length-1 vectors.
    a_vec = np.atleast_1d(a).ravel()
    c_vec = np.atleast_1d(c).ravel()

    nvecs = len(a_vec)
    x = np.zeros((n, nvecs))
    z = np.zeros((n, nvecs))

    for vi in range(nvecs):
        x[:, vi] = np.linspace(-a_vec[vi], a_vec[vi], n)
        # assuming vector centered in y (ie y=0)
        z[:, vi] = np.sqrt(c_vec[vi] ** 2 * (1 - x[:, vi] ** 2 / a_vec[vi] ** 2))

    return S, P, x, z

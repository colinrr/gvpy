"""Generic math helper functions, not tied to a specific inflow/outflow
model component.

Translated from MATLAB source: GV_for_claude/utilities/smoothUnion.m
"""

from collections.abc import Callable

import attrs
import numpy as np
from numpy.typing import ArrayLike

Kernel = Callable[[float | np.ndarray], float | np.ndarray]


def _quadratic_kernel() -> Kernel:
    return lambda x: (x * (2 + x) + 1) / 4


def _cubic_kernel() -> Kernel:
    return lambda x: (1 + 3 * x * (x + 1) - np.abs(x**3)) / 6


def _quartic_kernel() -> Kernel:
    return lambda x: (x + 1) ** 2 * (3 - x * (x - 2)) / 16


def _circular_kernel() -> Kernel:
    return lambda x: 1 + 0.5 * (x - np.sqrt(2 - x**2))


_KERNELS = {
    "quadratic": _quadratic_kernel,
    "cubic": _cubic_kernel,
    "quartic": _quartic_kernel,
    "circular": _circular_kernel,
}


def smooth_union(a: ArrayLike, b: ArrayLike, k: float, kernel: str = "quadratic") -> np.ndarray:
    """Return a smooth minimum of two curves a and b.

    k      = normalization distance (equals the maximum distance between
              smoothed and raw curves at the join where a=b)
    kernel = 'quadratic', 'cubic', 'quartic', or 'circular'

    Uses a smooth-minimum functions from a family of "clamped difference"
    kernels.
    See https://iquilezles.org/articles/smin/
    """
    a = np.asarray(a)
    b = np.asarray(b)
    assert a.shape == b.shape, "a and b vectors must match in size."
    # a and b should also be column vectors...

    if kernel not in _KERNELS:
        raise ValueError("Smoothing kernel not recognized.")
    g = _KERNELS[kernel]()

    k = k / g(0)
    # TRANSLATION NOTE: MATLAB used max([k - abs(a-b), zeros(size(a))],[],2) and
    # min([a b],[],2) - a concatenate-then-max/min-along-dim-2 idiom for elementwise max/min.
    # Translated to the equivalent np.maximum/np.minimum directly - same logic, no MATLAB
    # equivalent idiom needed in NumPy.
    h = np.maximum(k - np.abs(a - b), 0) / k
    smin = np.minimum(a, b) - k * g(h - 1.0)

    return smin


# Smooth top hat rise kernels: 0 -> 1 over x in [0, length]. Tails use the mirrored rise.
def _cosine_rise(x: np.ndarray, length: float, tol: float) -> np.ndarray:
    return 0.5 * (1 - np.cos(np.pi * x / length))


def _logistic_rise(x: np.ndarray, length: float, tol: float) -> np.ndarray:
    # k set so the curve runs from tol at x=0 to 1-tol at x=length
    k = 2 * np.log((1 - tol) / tol) / length
    return 1 / (1 + np.exp(-k * (x - length / 2)))


def _smoothstep_rise(x: np.ndarray, length: float, tol: float) -> np.ndarray:
    s = x / length
    return s**2 * (3 - 2 * s)


_TOP_HAT_KERNELS = {
    "cosine": _cosine_rise,
    "logistic": _logistic_rise,
    "smoothstep": _smoothstep_rise,
}


@attrs.frozen
class SmoothTopHat:
    """Piecewise smoothed top hat: smooth rise 0 -> 1, flat top at 1, smooth decay to 0.

    x = 0 is the start of the rise; the function ends at x_end, where the tail reaches 0
    (or tol, for the logistic kernel). Zero outside [0, x_end].

    rise_kernel/tail_kernel = 'cosine', 'logistic', or 'smoothstep'
    tol = logistic kernel value at its segment ends (its steepness k follows from length and
          tol), so a logistic segment jumps by tol where it meets the flat top
    y_max = dimensional maximum value (multiplies y)
    area  = dimensional total area; x is stretched so the area after scaling by y_max equals
            it. None leaves x unscaled.
    x_0   = dimensional x at the start of the rise
    """

    rise_length: float = 1.0
    flat_length: float = 2.0
    tail_length: float = 2.0
    rise_kernel: str = attrs.field(default="cosine", validator=attrs.validators.in_(_TOP_HAT_KERNELS))
    tail_kernel: str = attrs.field(default="logistic", validator=attrs.validators.in_(_TOP_HAT_KERNELS))
    tol: float = 1e-3
    y_max: float = 1.0
    area: float | None = None
    x_0: float = 0.0

    def x_end(self, dimensionalize: bool = True) -> float:
        """Return the x at which the function ends."""
        x_end = self.rise_length + self.flat_length + self.tail_length
        if dimensionalize:
            return self.x_0 + self.x_scale * x_end
        return x_end

    @property
    def x_scale(self) -> float:
        """Dimensional x per non-dimensional x."""
        if self.area is None:
            return 1.0
        return self.area / (self.y_max * self.get_area(dimensionalize=False))

    def __call__(self, x: ArrayLike, dimensionalize: bool = True) -> np.ndarray:
        x = np.asarray(x, dtype=float)
        if dimensionalize:
            return self.y_max * self((x - self.x_0) / self.x_scale, dimensionalize=False)
        x_flat = self.rise_length
        x_tail = self.rise_length + self.flat_length
        x_end = self.x_end(dimensionalize=False)
        rise = _TOP_HAT_KERNELS[self.rise_kernel]
        tail = _TOP_HAT_KERNELS[self.tail_kernel]
        return np.piecewise(
            x,
            [(x >= 0) & (x < x_flat), (x >= x_flat) & (x <= x_tail), (x > x_tail) & (x <= x_end)],
            [
                lambda x: rise(x, self.rise_length, self.tol),
                1.0,
                lambda x: tail(x_end - x, self.tail_length, self.tol),
                0.0,
            ],
        )

    def get_area(self, dimensionalize: bool = True) -> float:
        """Return the total area under the curve."""
        # Every kernel is point-symmetric about its segment midpoint, so each smooth segment
        # integrates to exactly length/2
        area = self.rise_length / 2 + self.flat_length + self.tail_length / 2
        if dimensionalize:
            return self.y_max * self.x_scale * area
        return area

    def get_xy(
        self, n: int = 101, dx: float | None = None, dimensionalize: bool = True
    ) -> tuple[np.ndarray, np.ndarray]:
        """Return x, y over [0, x_end]: n evenly spaced points, or a fixed step dx if given
        (last point then within dx/2 of x_end). With dimensionalize, x runs from x_0 to the
        dimensional x_end, and dx is dimensional."""
        x_start = self.x_0 if dimensionalize else 0.0
        x_end = self.x_end(dimensionalize=dimensionalize)
        if dx is None:
            x = np.linspace(x_start, x_end, n)
        else:
            x = x_start + np.arange(0, x_end - x_start + 0.5 * dx, dx)
        return x, self(x, dimensionalize=dimensionalize)

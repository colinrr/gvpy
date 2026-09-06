"""Generic math helper functions, not tied to a specific inflow/outflow
model component.

Translated from MATLAB source: GV_for_claude/utilities/smoothUnion.m
"""

import numpy as np


def _quadratic_kernel():
    return lambda x: (x * (2 + x) + 1) / 4


def _cubic_kernel():
    return lambda x: (1 + 3 * x * (x + 1) - np.abs(x**3)) / 6


def _quartic_kernel():
    return lambda x: (x + 1) ** 2 * (3 - x * (x - 2)) / 16


def _circular_kernel():
    return lambda x: 1 + 0.5 * (x - np.sqrt(2 - x**2))


_KERNELS = {
    "quadratic": _quadratic_kernel,
    "cubic": _cubic_kernel,
    "quartic": _quartic_kernel,
    "circular": _circular_kernel,
}


def smooth_union(a, b, k, kernel="quadratic"):
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

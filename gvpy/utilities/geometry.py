"""Generic geometry helper functions, not tied to a specific inflow/outflow
model component.

Translated from MATLAB source: GV_for_claude/utilities/geometry/
"""

import numpy as np


def get_elliptical_cylinder_sa(a, c, l=1, halved=True):
    """Ellipse perimeter approximation of Ramanujan.

    a, c = axes lengths. a is nominally semi-major but not critical
    l    = cylinder length. Default 1 to get perimeter only
    halved = bool. True [default] = half cylinder, total area/2
    """
    h = (a - c) ** 2 / (a + c) ** 2
    sa = np.pi * (a + c) * (1 + 3 * h / (10 + np.sqrt(4 - 3 * h))) * l

    if halved:
        sa = sa / 2

    return sa

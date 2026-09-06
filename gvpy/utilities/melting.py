"""Helper functions specific to ice-melting/heat-transfer geometry (the
central cauldron control-volume "vent" physics), as opposed to a specific
inflow/outflow drainage component.

Translated from MATLAB source: GV_for_claude/utilities/getEllipticHeatIntensity.m
"""

import numpy as np
from scipy.integrate import quad


def get_elliptic_heat_intensity(alpha):
    """Returns bar_alpha, the ratio of horizontal melt rate to average melt rate
    as a function of alpha.

    alpha = vertical/horizontal heat transfer ratio
    """
    # TODO: build in (or around):
    #   - the average intensity, values at both 0 and pi/2
    #   - u_melt --> da/dt
    #   -

    def int_fun(theta):
        return alpha / (np.sin(theta) ** 2 + alpha**2 * np.cos(theta) ** 2) ** (1 / 2)

    # TRANSLATION NOTE: MATLAB's integral() maps to scipy.integrate.quad here, which returns
    # (value, error_estimate) rather than just the value - only the value is used here, matching
    # MATLAB's single-output integral() call.
    integral_value, _ = quad(int_fun, 0, np.pi / 2)
    bar_alpha = 2 / np.pi * integral_value

    return bar_alpha

"""Subglacial drainage model component: standalone functions supporting
IceCauldron's pressure-wave/flood-drainage methods (get_l_lambda,
solve_delta_p).

Translated from MATLAB source: @IceCauldron/get_L_lambda.m's embedded local
function (get_x_equals_l_lambda), and @IceCauldron/solveDeltaP.m's embedded
local functions (objectiveDeltaP_1, Delta_P_from_Bernoulli_1,
Delta_P_from_U_tip).

TRANSLATION NOTE: per project direction, the ENTIRE subglacial drainage
component is explicitly in-development/test-only, not finalized physics -
e.g. IceCauldron.get_drainage_fluxes (the actual subglacial + inter-cauldron
drainage flux calculation) is not implemented in the MATLAB source at all,
only returning dummy zero fluxes. solveDeltaP itself only implements one
regime (x <= l_lambda) and explicitly errors otherwise. Bugs/incompleteness
here are preserved as-is per CLAUDE.md's in-development-components policy,
not silently fixed.
"""

import numpy as np
from scipy.optimize import fsolve


def get_x_equals_l_lambda(cauldron):
    """x_crit = get_x_equals_l_lambda(obj)
    Find the coordinate x at which x = l_lambda as defined by 1/4 * flexural
    wavelenth
    """
    h_i0, _, _, theta_i, theta_b = cauldron.get_elevation_profiles(0)

    A = np.tan(theta_b) - np.tan(theta_i)

    # [~,E_prime] = getFluxuralRigidity(obj,h_i0);
    l_const = np.pi / 2 * (cauldron.E_prime / (12 * cauldron.CONSTANTS.G * cauldron.CONSTANTS.RHO_ICE)) ** (1 / 4)

    # Polynomial in x to get l_lambda_crit
    def ff(x):
        return (x / l_const) ** (4 / 3) - x * A - h_i0

    _, _, l_lambda_0 = cauldron.get_l_lambda(0)  # Guess from starting val

    # TRANSLATION NOTE: MATLAB's fzero(ff, x0) finds a root near an initial guess x0 (not a
    # bracket) - scipy.optimize.fsolve is the closest match for that calling convention (unlike
    # brentq/bisect, which need a bracketing interval). Flagging per CLAUDE.md's ODE/solver
    # translation notes since this is a MATLAB->SciPy root-finder substitution.
    x_crit = fsolve(ff, l_lambda_0)[0]

    assert 0 <= x_crit <= cauldron.l_G, "Critical x to match l_lambda does not fall within bounds."

    return x_crit


def objective_delta_p_1(cauldron, u, P_0, rho_f, x):
    """d_Delta_P = objectiveDeltaP_1(obj,u,P_0,rho_f,x)
    Case of early propagating crack tip:
      -> u_1 = u_tip
      -> (bar_h_s / h_s(L_s)) ~ 1 (average crack opening ~ inlet opening)
      -> u_0 = 4 * u_tip * (bar_h_s / h_s(L_s)) = 4 * u_tip
    """
    # TRANSLATION NOTE: MATLAB references h_i and l here, but neither is a parameter of this
    # function nor defined locally within it - a genuine bug in the MATLAB source (local
    # functions in a MATLAB function file don't share workspace with solveDeltaP's h_i/l).
    # Preserved as-is: h_i and l are undefined names below, so this raises NameError if ever
    # called, matching MATLAB's "Undefined function or variable" error. Not guessing at what was
    # intended (e.g. threading x through get_elevation_profiles/get_l_lambda) per the translation
    # skill's guidance on genuinely unclear broken code.

    # ---- Use Tsai and Rice to find del_h, DeltaP ----
    Delta_P_TR = delta_p_from_u_tip(cauldron, u, rho_f, x)
    del_h_TR = cauldron.get_del_h(Delta_P_TR, h_i, l)  # noqa: F821

    # ---- Now use Bernoulli equation ----
    Delta_P_B = delta_p_from_bernoulli_1(cauldron, u, P_0, rho_f, x, del_h_TR)

    d_Delta_P = abs(Delta_P_B - Delta_P_TR)

    return d_Delta_P


def delta_p_from_bernoulli_1(cauldron, u, P_0, rho_f, x, h_s):
    """Delta_P_B = Delta_P_from_Bernoulli_1(obj,u,P_0,rho_f,x, h_s)
    Case of early propagating crack tip:
      -> u_1 = u_tip
      -> (bar_h_s / h_s(L_s)) ~ 1 (average crack opening ~ inlet opening)
      -> u_0 = 4 * u_tip * (bar_h_s / h_s(L_s)) = 4 * u_tip
    """
    # Get ice thickness, elevations, static pressures at beginning and end of profile
    h_i, _, z_b, _, _ = cauldron.get_elevation_profiles(np.array([0, x]))
    P_i = cauldron.CONSTANTS.RHO_ICE * cauldron.CONSTANTS.G * h_i
    # Elevation drops - ice and bedrock
    delta_h_i = h_i[0] - h_i[1]
    delta_z_b = z_b[0] - z_b[1]

    # Friction factor
    # --> coefficient
    # f_0 = obj.Constants.g .* obj.manning_roughness.^2 ./ (2^(2/3));
    f_0 = cauldron.CONSTANTS.G * cauldron.manning_roughness**2 * (2 ** (7 / 3))
    # --> full factor
    f = f_0 * x / h_s ** (4 / 3)

    # Excess pressure at source
    DeltaP_0 = P_0 - P_i[0]

    # Delta P from Bernoulli
    Delta_P_B = (
        DeltaP_0
        + cauldron.CONSTANTS.RHO_ICE * cauldron.CONSTANTS.G * delta_h_i
        + cauldron.CONSTANTS.G * rho_f * delta_z_b
        + rho_f / 2 * (15 - 5 / 2 * f) * u**2
    )

    return Delta_P_B


def delta_p_from_u_tip(cauldron, u_tip, rho, x):
    """Delta_P = Delta_P_from_U_tip(obj,u_tip,rho,x)
    This function for crack tip velocity is from Tsai & Rice (2012) (Journal
    of Applied Mechanics). It is known to be valid for L_hi_ratio <= 5. It is
    unknown how valid the relation is above that point.
     dP   = driving fluid overpressure
     rho  = fluid density
     x    = position of crack tip along glacier profile
    """
    #     [h_i, ~, ~, ~,~] = obj.getElevationProfiles(x);
    #     [~,E_prime] = obj.getFluxuralRigidity(h_i);

    l, l_hi_ratio, _ = cauldron.get_l_lambda(x)

    phi = 5.13 + 0.64 * l_hi_ratio + 0.94 * l_hi_ratio**2

    k = (cauldron.manning_roughness / 0.038) ** 6  # Nikuradse roughness height

    Delta_P = (u_tip / (phi * rho ** (-1 / 2) * (cauldron.E_prime) ** (-2 / 3) * (l / k) ** (1 / 6))) ** (6 / 7)

    return Delta_P

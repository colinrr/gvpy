"""Physical constants used throughout the glaciovolcanic model.

Translated from MATLAB source: @IceCauldron folder property `Constants`,
class defined in GV_for_claude/ThermoConstants.m
"""


class ThermoConstants:
    # WATER
    RHO_W = 1000  # (kg/m3)  density of liquid water
    RHO_ICE = 917  # (kg/m3)  density of glacier ice
    L_ICE = 334e3  # (J/kg)   Ice/water latent heat of fusion
    C_W = 4184  # (J/kg/K) water heat capacity
    G = 9.81  # (m/s) gravity

    # MAGMA PARAMS
    C_M = 1350

    # rho_t   = 1000; % Tomasson 1996 value for tephra on Myrdalsandur, in pretty good agreement with model results
    # T_m     = 1450;    % K

    # rho_m   = 2600; % We need this density for the appropriate conduit heat/mass flux. Can use this density or a reduced tephra density for flood discharge mass
    #                 # Magma basalt lava-- unfragmented (kg/m^3); (Gudmundsson, 2003)

    # Can use a second bulk tephra density for volume in vault later? This
    # is global value, in that case
    # C.rho_t   = 1000; % Tomasson 1996 value for tephra on Myrdalsandur, in pretty good agreement with model results


# TRANSLATION NOTE: applied UPPER_SNAKE_CASE per CLAUDE.md's PEP8 Formatting section, since
# these are true constants (never reassigned) - name identity preserved (rho_w -> RHO_W, etc.).
# CLAUDE.md's Project Context notes this will likely become part of a larger config ecosystem
# later; kept as a direct, minimal translation for now.

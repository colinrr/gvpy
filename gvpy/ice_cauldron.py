"""Central IceCauldron control-volume class: parameters, initialization, and
methods for the model's core control-volume geometry and bookkeeping.

Translated from MATLAB source: @IceCauldron/IceCauldron.m

Additional methods from the @IceCauldron classdef folder will be folded into
this same module as they're translated (see CLAUDE.md's Project Structure
section - the classdef folder maps to a single ice_cauldron.py module). This
pass covers IceCauldron.m itself, plus get_elevation_profiles (needed at
construction time - see TRANSLATION NOTE in __attrs_post_init__).

TRANSLATION NOTE (general, applies throughout this file): MATLAB header
comments that simply restate a function's call signature (e.g.
"% V = get_CV_control_volume(obj,a)") are dropped when converting the header
comment block to a Python docstring, since that information (name, inputs,
outputs) is already visible in the `def` line itself. Descriptive text beyond
the signature echo is preserved.
"""

import json
from typing import ClassVar, Optional, Sequence, Union

import attrs
import numpy as np
import xarray as xr
from scipy.optimize import minimize, minimize_scalar

from .plume import predict_forest

# SUBGLACIAL DRAINAGE (IN DEVELOPMENT): supports get_l_lambda/solve_delta_p below - the entire
# subglacial drainage component is explicitly in-development/test-only, not finalized physics.
from .subglacial import (
    delta_p_from_bernoulli_1,
    delta_p_from_u_tip,
    get_x_equals_l_lambda,
    objective_delta_p_1,
)
from .utilities.constants import ThermoConstants
from .utilities.geometry import get_elliptical_cylinder_sa
from .utilities.math_helpers import smooth_union  # used by get_l_lambda (subglacial drainage, in development)
from .utilities.melting import get_elliptic_heat_intensity
from .utilities.validators import (
    allow_none,
    attrs_range_inclusive,
    attrs_validate_nonnegative,
    attrs_validate_pos_or_nan,
    attrs_validate_positive,
)

Number = Union[float, int]
PerCauldronInput = Union[Number, Sequence[Number]]

# TRANSLATION NOTE: property validation is now handled via attrs.field(validator=...) below,
# using the validator functions in utilities/validators.py (translated from IceCauldron.m's
# local validators plus MATLAB's built-in mustBeNonnegative/mustBePositive). This resolves the
# design decision flagged in an earlier translation pass (hand-written __post_init__ checks vs.
# a validation library): attrs was chosen and added as a project dependency (see pyproject.toml)
# since MATLAB validates most of this properties block (~28 fields) declaratively, one validator
# per field - attrs.field(validator=...) preserves that same per-field, declarative structure
# and line-for-line traceability to the MATLAB source, which hand-written checks would not.


@attrs.define
class IceCauldron:
    """Cauldron class to hold simulation parameters.
     - geometry
     - heat transfer and material properties
     - physical constants
     - model switches
    """

    CONSTANTS: ClassVar = ThermoConstants
    # TRANSLATION NOTE: applied UPPER_SNAKE_CASE per CLAUDE.md's PEP8 rules for true constants
    # (name identity preserved: Constants -> CONSTANTS).

    # ---------------------------
    # CAULDRON SPECIFIC PARAMS
    # ---------------------------
    # ---- Cauldron Geometry ------
    n_cauldrons: int        = 1  # Number of discrete cauldrons
    geometry: str           = attrs.field(default="cylinder", # Cauldron geometry approximation - half cylinder or half spheroid. SPHEROID NYI
                                          validator=attrs.validators.in_(("cylinder", "spheroid")))  # {'cylinder','spheroid'}
    L_n: PerCauldronInput   = attrs.field(default   = 1000.0, # (m) Length of fissure segments in each cauldron
                                          validator = attrs_validate_pos_or_nan)
    G_n: PerCauldronInput   = attrs.field(default   = 300.0, # (m) Average cauldron ice thickness.
                                          validator = attrs_validate_nonnegative)
    # - note 270-300m from Larsen 2021-- ice thickness, east side
    a_0: Number             = attrs.field(default   = 10.0,
                                          validator = attrs_validate_nonnegative)  # (m) Initial radius of bottom heat transfer surface
    z_b_n: PerCauldronInput = 750.0  # (m) RELATIVE bedrock elevation at cauldron vent
    # c_0 double      % Vertical axis initial value (m)
    # a double        %
    # c double

    # ---- Drainage Glacier Geometry ----
    L_d_n: PerCauldronInput = attrs.field(default   = 2500.0, # (m) subglacial (horizontal) path length from cauldron n to main drainage inlet
                                          validator = attrs_validate_positive)
    l_G: Number             = attrs.field(default   = 1.4e4, # (m) Horizontal length of outlet glacier path
                                          validator = attrs_validate_nonnegative)
    z_b_0: Number           = attrs.field(default   = 350.0, # (m) Bedrock RELATIVE elevation of main drainage outlet
                                          validator = attrs_validate_nonnegative)
    # NOTE: final bedrock elevation at terminus is taken as 0
    h_i_0: Number           = attrs.field(default   = 550.0, # (m) Ice thickness above main drainage inlet
                                          validator = attrs_validate_nonnegative)
    h_i_f: Number           = attrs.field(default   = 0.0, # (m) Ice thickness above main drainage outlet (glacier terminus)
                                          validator = attrs_validate_nonnegative)
    b: Number               = attrs.field(default   = 1.5e3, # (m) Average width of the subglacial drainage sheet
                                          validator = attrs_validate_nonnegative)

    # ---- Ice Properties ------
    T_f: Number = attrs.field(default   = 273.15, # (K) melting temperature of ice, assuming no pressure effect
                              validator = attrs_validate_nonnegative)
    T_i: Number = attrs.field(default   = 273.15, # (K) Initial ice temperature
                              validator = attrs_validate_nonnegative)
    A: Number   = attrs.field(default   = 2.3e-24, # (TODO) Glen parameter
                              validator = attrs_validate_positive)
    E: Number   = attrs.field(default   = 5e9, # (Pa) Young's modulus of glacier ice. [0.9-10]e9 is a good range
                              validator = attrs_validate_nonnegative)
    nu: Number  = attrs.field(default   = 0.3, # Poisson's ratio for glacier ice. Values 0.3 to 0.4 are frequently used
                              validator = attrs_validate_nonnegative)

    # ----- Magma flux & heat transfer ---------
    # Estimated total airborne tephra volume (used as a stopping condition)
    V_t_max: Number      = attrs.field(default   = 0.95e9 * 1.2 * 0.5,
                                       validator = attrs_validate_nonnegative)

    Q_decay_frac: Number = attrs.field(default   = 0.3, # (-) Remaining fraction of airborne tephra after which to force plume flux decay
                                       validator = attrs_validate_nonnegative)
    Q_0: Number          = attrs.field(default   = 5e7, # (kg/s) Total magma discharge rate
                                       validator = attrs_validate_nonnegative)
    alpha: Number        = attrs.field(default   = 1.0, # (-) % Vertical/radial ratio of heat transfer rate
                                       validator = attrs_validate_positive)
    f_i: Number          = attrs.field(default   = 0.6, # (-) Parameterized heat transfer efficiency: magma->ice
                                       validator = attrs_range_inclusive(0, 1))
    T_m: Number          = attrs.field(default   = 1450.0, # (K) Magma temperature
                                       validator = attrs_validate_nonnegative)
    rho_m: Number        = attrs.field(default   = 2600.0, # (kg/m3) Magma density
                                       validator = attrs_validate_nonnegative)
    rho_p: Number        = attrs.field(default   = 1000.0, # (kg/m3) Erupted pyroclast density
                                       validator = attrs_validate_nonnegative)
    rho_r: Number        = attrs.field(default   = 2500.0, # (kg/m3) Lithic rock density
                                       validator = attrs_validate_nonnegative)

    # ----- Flood/discharge properties ---------
    fixed_phi_p: Number       = attrs.field(default   = 0.4, # (-) Fixed volume fraction of pyroclasts in flood (flood_density_model = 2,3 only)
                                            validator = attrs_range_inclusive(0, 1))
    manning_roughness: Number = attrs.field(default   = 0.033, # (m^(-1/3) s) Manning's friction coefficient (0.01 - 0.1 is a reasonable range)
                                            validator = attrs_validate_nonnegative)

    # ----- Calculated/determined internally - generally don't touch -----
    open_cauldron: PerCauldronInput = False  # Whether or not cauldron roof is open
    ice_free_cauldron: PerCauldronInput = False  # Whether initial ice mass has been completely melted

    # TRANSLATION NOTE: open_cauldron/ice_free_cauldron are per-cauldron mutable simulation
    # state that changes during the ODE solve via events (see CLAUDE.md's Data Structures
    # section) - eventually these should move to a separate solver-state structure rather than
    # living on IceCauldron itself. Keeping them here as plain broadcast NumPy arrays for now
    # since several not-yet-translated methods (get_u_melt, get_f_i,
    # check_open/ice_free_cauldron_conditions, get_material_heights) read them directly as
    # obj.openCauldron / obj.iceFreeCauldron - revisit once those are translated and the
    # solver-state design is worked out.

    chi: Optional[float]        = attrs.field(init=False, default=None)  # (m3/kg) MAX volume ice meltable per kg magma
    t_Qdecay: Optional[float]   = attrs.field(init=False, default=None)  # (s) Time at which eruption mass flux begins to decay
    # TRANSLATION NOTE: TODO t_Qdecay is declared but never assigned within IceCauldron.m's
    # constructor or its other methods - presumably set elsewhere (gvMain.m / glaciovolcano.m,
    # not yet translated). Left as None here; flagging rather than guessing at its assignment
    # logic.
    tau_ND: Optional[float]      = attrs.field(init=False, default=None)  # (s) Non-dimensional time scale
    V_ND: Optional[float]        = attrs.field(init=False, default=None)  # (m^3) Non-dimensional volume scale
    L_ND: Optional[float]        = attrs.field(init=False, default=None)  # (m) Non-dimensional length scale
    theta_i: Optional[float]     = attrs.field(init=False, default=None)  # (rad) Slope of outlet glacier ice surface
    theta_b: Optional[float]     = attrs.field(init=False, default=None)  # (rad) Slope of outlet glacier bed surface
    k_nikuradse: Optional[float] = attrs.field(init=False, default=None)  # (m) Nikuradse roughness height
    E_prime: Optional[float]     = attrs.field(init=False, default=None, validator=allow_none(attrs_validate_nonnegative))  # Flexural parameter - Young's modulus adjusted with poisson's ratio
    # TRANSLATION NOTE: E_prime carries MATLAB's {mustBeNonnegative} validator, but unlike the
    # input parameters above, it's computed (not user-supplied) and starts as None until
    # __attrs_post_init__ assigns it - allow_none() lets it skip validation while unset, then
    # validates normally once __attrs_post_init__ sets its real value (attrs.define validates on
    # every subsequent assignment by default, not just at construction).

    fastest_cauldron_index: Optional[int] = attrs.field(init=False, default=None)
    slowest_cauldron_index: Optional[int] = attrs.field(init=False, default=None)

    # ---------------------------------------
    #    SIMULATION PARAMETERS AND SWITCHES
    # ---------------------------------------
    non_dimensionalized: bool = False  # If we can get a non-dim version running
    t_final: Number = 3 * 3600  # End time (TODO: in normalized units of tau_m - prob becomes dynamic with Q history)
    timesteps: int = 2001  # Output timesteps

    # Switch for flood fluid density model (see IceCauldron.drainageDensity.m)
    flood_density_model: int = attrs.field(default=4, validator=attrs_range_inclusive(1, 6))  # {1,...,6}

    # Turn on/off ice flow into cauldrons - ice flow velocity fixed to 0
    ice_inflow_mode: str = "off"  # "fixed-glen", "off"

    # Turn on/off sup. flow out of cauldrons
    supraglacial_drainage_mode: str = "overflow"  # "overflow", "off"

    # If true, all plume material fluxes (water, pyroclasts) and dimensions will be fixed to 0
    disable_plume_flux: bool = True
    terminal_condition: str = "default"  # "default", "t_final"
    plume_emulator_json: str = "hydroplume_emulator_2025-04-22.json"  # JSON containing random forest emulator parameters
    plume_emulator: Optional[dict] = attrs.field(init=False, default=None)  # Structure array containing random forest parameters

    # units struct

    # ---- Hidden (bookkeeping) fields ----
    # TRANSLATION NOTE: these were MATLAB `properties (Hidden)` - per-instance properties with
    # static content never mutated per-instance anywhere in the source. Translated as ClassVar
    # tuples (shared, immutable) rather than per-instance dataclass fields, matching actual
    # usage. per_cauldron_fields lists real attribute names, so its entries are renamed to
    # match this file's snake_case identifiers (open_cauldron, ice_free_cauldron). By contrast,
    # vector_solution_vars/vector_derived_vars/etc. are symbolic solver-variable-name labels,
    # not real Python attributes on this class (e.g. 'horizontalMeltingRate_n' has no
    # corresponding self.horizontalMeltingRate_n) - these are left in their original MATLAB
    # casing since they must stay consistent with matching string keys in not-yet-translated
    # code (e.g. getVarLabels.m uses the same camelCase keys).
    per_cauldron_fields: ClassVar[tuple] = ("L_n", "G_n", "z_b_n", "L_d_n", "open_cauldron", "ice_free_cauldron")
    # Vars directly solved in ode integrator
    vector_solution_vars: ClassVar[tuple] = ("a_n", "c_n", "V_ice_n", "V_w_n", "V_p_n")  # Per cauldron
    scalar_solution_vars: ClassVar[tuple] = ()  # 1 each
    # Vars derived from ode solution
    vector_derived_vars: ClassVar[tuple] = (
        "V_cavity_n",
        "horizontalMeltingRate_n",
        "verticalMeltingRate_n",
        "f_i_n",
        "V_w_plus_p",
        "V_cum",
        "V_CV_n",
        "H_w_n",
        "H_i_n",
        "H_p_n",
        "H_cum_n",
    )
    scalar_derived_vars: ClassVar[tuple] = ()
    # Helpers determined internally
    # (fastest_cauldron_index, slowest_cauldron_index are real per-instance fields, above)

    params: Optional[xr.Dataset] = attrs.field(init=False, default=None)

    def __attrs_post_init__(self):
        """Initiate IceCauldron object."""
        # TRANSLATION NOTE: MATLAB's constructor uses a generic reflective loop over all class
        # properties to assign user-supplied opts, falling back to class defaults. Python's
        # dataclass-generated __init__ already handles "user value or default" assignment for
        # every field via normal keyword arguments - that portion of the MATLAB loop needs no
        # Python equivalent. What follows is the per-cauldron broadcast/validation behavior
        # from that same loop, applied explicitly per field, plus the secondary-constants
        # calculations that followed it in the MATLAB source.
        # (attrs equivalent of dataclasses' __post_init__ is __attrs_post_init__.)

        # ---- Build per-cauldron parameter Dataset (user-supplied + computed secondary) ----
        cauldron = np.arange(self.n_cauldrons)

        data = {}
        for name in ("L_n", "G_n", "z_b_n", "L_d_n"):
            value = np.asarray(getattr(self, name), dtype=float)
            if value.ndim == 0:
                value = np.full(self.n_cauldrons, value)
            else:
                self.validate_vector_length(value, self.n_cauldrons, "IceCauldron", f"input parameter: {name}")
            data[name] = ("cauldron", value)
        self.params = xr.Dataset(data, coords={"cauldron": cauldron})

        # State fields - broadcast to per-cauldron arrays, kept separate from self.params (see
        # TRANSLATION NOTE on open_cauldron/ice_free_cauldron above).
        for name in ("open_cauldron", "ice_free_cauldron"):
            value = np.asarray(getattr(self, name))
            if value.ndim == 0:
                value = np.full(self.n_cauldrons, value)
            else:
                self.validate_vector_length(value, self.n_cauldrons, "IceCauldron", f"input parameter: {name}")
            setattr(self, name, value)

        # ---- Calculate secondary constants ----

        # Calculate cauldron magma flux and check geometry
        L_n = self.params["L_n"].values
        self.params["Q_n"] = ("cauldron", self.Q_0 * (L_n / L_n.sum()))
        if self.geometry == "cylinder":
            self.params["Q_per_L"] = ("cauldron", self.params["Q_n"].values / L_n)
        elif self.geometry == "spheroid":
            raise NotImplementedError("Spheroid geometry is not currently implemented")
            # TRANSLATION NOTE: MATLAB set obj.Q_per_L = nan and obj.L_n = nan immediately
            # before this error() call. Since MATLAB's error() halts execution and the object
            # is never returned, those assignments have no observable effect - dropped as dead
            # code rather than translated, since keeping them would add lines with no behavior
            # change. Flagging per Strict Rule #3 since this is a (behavior-neutral)
            # simplification, not a pure syntax conversion.

        C = self.CONSTANTS

        # Chi: m3 meltable ice per kg erupted magma
        self.chi = C.C_M * (self.T_m - self.T_i) / (C.RHO_ICE * C.L_ICE)
        self.params["melt_vol_per_s"] = ("cauldron", self.chi * self.params["Q_n"].values)

        _, _, _, self.theta_i, self.theta_b = self.get_elevation_profiles(1)
        # TRANSLATION NOTE: get_elevation_profiles is translated here (ahead of the planned
        # dependency order for the remaining @IceCauldron method files) because the MATLAB
        # constructor calls it directly (obj.getElevationProfiles(1)) - it's a hard dependency
        # for __post_init__ to run at all. All other @IceCauldron methods remain deferred to
        # the next translation pass as planned.

        self.get_dimensional_scales()

        # Friction factors
        self.k_nikuradse = (self.manning_roughness / 0.038) ** 6

        # Flexural parameter
        self.E_prime = self.E / (1 - self.nu**2)

        if not self.disable_plume_flux:
            self.plume_emulator = self.load_emulator()

    # TODO: build constants and units functions
    def get_dimensional_scales(self):
        """Calculate main dimensional scales for equation non-dimensionalization.
        TODO: non-dimensionalization not yet implemented
        Current setup:
          - Let V_ND = V_melt : volume of semi-cylindrical cavity w/
            vertical axis equal to ice thickness and horizontal axis
            set by the heat transfer ratio, alpha
          - tau_ND = tau_melt : the time to melt out the above volume
            at constant melt rate
          - L_ND  = G_0 : ice thickness over the cauldron (= cavity
            height)
          - For multi-cauldron case, tau_melt is taken as the minimum
            across all cauldrons, with corresponding volumes and length
            scales
        """
        G_n = self.params["G_n"].values
        L_n = self.params["L_n"].values
        Q_n = self.params["Q_n"].values

        V_melt = np.pi / 2 * G_n**2 * 1 / self.alpha * L_n  # Ignore any a_0 for now
        tau_melt = V_melt / (self.f_i * self.chi * Q_n)

        self.params["V_melt"] = ("cauldron", V_melt)
        self.params["tau_melt"] = ("cauldron", tau_melt)

        self.fastest_cauldron_index = int(np.argmin(tau_melt))
        self.slowest_cauldron_index = int(np.argmax(tau_melt))
        # TRANSLATION NOTE: MATLAB's [~,idx] = min(...)/max(...) returns a 1-based index;
        # np.argmin/np.argmax return 0-based. This index is used directly, 0-based, throughout
        # the Python translation - a deliberate, consistent convention change (not a bug), but
        # flagging since any MATLAB-side documentation assuming 1-based indexing here no
        # longer applies.

        self.V_ND = V_melt[self.fastest_cauldron_index]
        self.tau_ND = tau_melt[self.fastest_cauldron_index]
        self.L_ND = G_n[self.fastest_cauldron_index]

    def load_emulator(self):
        print("Loading plume emulator json...")
        with open(self.plume_emulator_json) as f:
            emulator = json.load(f)
        return emulator

    # ---------- FUNCTIONS USED IN INTEGRATION SOLVER ---------- #
    # These should return column vectors when nCauldrons > 1
    #  '-> watch when gv vector properties are called
    def get_initial_conditions(self):
        """Build solver initial conditions for:
        a_n
        c_n
        V_ice_n
        V_water_n
        V_pyroc_n
        E_water_n
        q_d
        h_s
        L_s
        """
        ic = np.concatenate(
            [
                self.a_0 * np.ones(self.n_cauldrons),  # a_n
                np.zeros(self.n_cauldrons),  # c_n
                self.get_cv_control_volume(np.full(self.n_cauldrons, self.a_0)),  # V_ice_n
                np.zeros(self.n_cauldrons),  # V_w_n
                np.zeros(self.n_cauldrons),  # V_pyroc_n
                # np.zeros(self.n_cauldrons),  # E_water_n - TODO
                # 0,                          # q_d
                # 0,                          # h_s
                # 0,                          # L_s
            ]
        )

        _, ns = self.get_solution_indices()
        assert len(ic) == ns, "Length of initial conditions vector does not match IceCauldron solution vars."

        return ic

    # ----------------------- CAULDRON-BASED EVENTS -----------------------
    def check_open_cauldron_conditions(self, c_n):
        # Brief dummy condition - roof drops below X% thickness
        min_roof_fraction = 0.6

        event_val = self.params["G_n"].values * min_roof_fraction - c_n
        is_term = np.ones_like(c_n)
        direction = -np.ones_like(c_n)

        # Cancel if already open
        event_val[self.open_cauldron] = 99

        return event_val, is_term, direction

    def check_ice_free_cauldron_conditions(self, V_ice_n, a_n):
        # Brief dummy condition - roof drops below 25% thickness
        event_val = V_ice_n / self.get_cv_control_volume(a_n) - 1e-5
        is_term = np.ones_like(V_ice_n)
        direction = -np.ones_like(V_ice_n)

        # Cancel if already ice free
        event_val[self.ice_free_cauldron] = 99

        return event_val, is_term, direction

    # --------------------------- EQUATIONS -------------------------------
    # --> Basic geometry and some related equations will go here. Most
    # equations however will have their own functions in the @IceCauldron
    # methods directory
    def get_cv_control_volume(self, a):
        """Get total volume (m3) of control volume."""
        a = np.asarray(a, dtype=float)
        # For cylindrical geometry:
        return (self.params["G_n"].values * self.params["L_n"].values) * 2 * a

    def get_cv_ice_area(self):
        """Get area of cauldron control volume (m2) vertical ICE surface
        for ice inflow and open cauldron melting.
        """
        # Excludes cauldron ENDS:
        return self.params["G_n"].values * self.params["L_n"].values * 2

    # -------------------- VENT/MELTING & ICE INFLOW -----------------------
    # Translated from @IceCauldron/get_u_ice.m, get_f_i.m, get_u_melt.m, get_da_dt.m
    def get_u_ice(self):
        """[u_ice_0, u_ice_bar] = get_u_ice(obj)
        Estimate sliding+creep inflow velocity into cauldron.
          u_ice_bar = average inflow velocity across the CV boundary (m/s)
          u_ice_0   = velocity at glacier base (i.e. sliding) (m/s)

        Behaviour is set by IceCauldron.ice_inflow_mode.
        """
        # TRANSLATION NOTE: MATLAB's switch has no otherwise/default case - if ice_inflow_mode is
        # neither 'fixed-glen' nor 'off', u_ice_0/u_ice_bar are left undefined in MATLAB (an error
        # at use). Preserved here via if/elif with no else - Python raises UnboundLocalError on
        # return if neither branch matches, the equivalent failure mode.
        # Also: MATLAB's obj.G_n' (transpose to row vector) has no Python equivalent needed here -
        # self.params["G_n"].values is already a flat 1-D array (see CLAUDE.md's Data Structures
        # section on resolving orientation ambiguity structurally).
        if self.ice_inflow_mode == "fixed-glen":
            # V0: crudely assume a maximum strain rate fitting
            # within Glen's law over a length scale ~2*G_n
            # - gives O(10's to ~100 m/day)
            # - defined with negative velocity for inflow into
            # cauldron
            max_strain_rate = 1e-6
            u_ice_bar = -max_strain_rate * self.params["G_n"].values * 2
            u_ice_0 = -u_ice_bar
        elif self.ice_inflow_mode == "off":
            u_ice_0 = 0
            u_ice_bar = 0

        return u_ice_0, u_ice_bar

    def get_f_i(self, a=None, c=None, t=None, events=None):
        """f_i_n = get_f_i(obj,a,c,t,events)
        Calculate the geometry-dependent heat transfer efficiency for
        ice melting.

        a,c = scalar or vector cauldron dimensions
              [n values x m cauldrons]
              Must be included together
              If these are not included, defaults to max f_i value

        iceFreeIndex = last index of a/c for each cauldron for which
                      the iceFree condition is false

         - TODO: may rewrite f_i at some point to be a function of
         water volume/ heat conservation
        """
        # TRANSLATION NOTE: MATLAB dispatches on nargin (0/2/4/5 extra args) to pick between the
        # default, per-step ODE, and full-output modes, and explicitly errors if called with
        # exactly 3 extra args (a,c,t but no events). Translated using None-default sentinels:
        # `t is None` covers both MATLAB's nargin<4 and isempty(t) cases (Python doesn't
        # distinguish "not passed" from "empty" the same way); `t` provided without `events`
        # falls through to the same "Could not parse input args" error as MATLAB.
        if a is None:
            return self.f_i * np.ones((self.n_cauldrons, 1))

        if t is None:
            full_output = False
            assert a.shape[0] == self.n_cauldrons, "a_n is not the correct size for ode solver."
            assert c.shape[0] == self.n_cauldrons, "c_n is not the correct size for ode solver."
        elif events is not None:
            full_output = True
            assert a.shape[0] == len(t) and c.shape[0] == len(t), (
                "a, c, and t vectors must be of equal length for full output."
            )
        else:
            raise ValueError("Could not parse input args.")

        f_i_n = np.zeros(a.shape)

        # Get which cauldrons are ice free at which time steps
        if full_output:
            # TRANSLATION NOTE: `events` is expected to carry the same information as MATLAB's
            # events table (columns indexName, cauldronIndex, t) - the Python ODE-events
            # machinery this depends on hasn't been translated yet, so this branch is untested
            # and its exact `events` structure (DataFrame? dict of arrays?) is not yet settled.
            ice_free = np.zeros((a.shape[0], self.n_cauldrons), dtype=bool)
            f_i_n = np.zeros((a.shape[0], self.n_cauldrons))
            for ci in range(self.n_cauldrons):
                event_idx = (events["indexName"] == "iceFreeCauldron") & (events["cauldronIndex"] == ci)
                ice_free[:, ci] = t >= events["t"][event_idx]
        else:
            ice_free = np.zeros(a.shape, dtype=bool)
            # Events timing not included, so can only assume f_i based
            # on current IceFree cauldron status
            ice_free[self.ice_free_cauldron] = True

            #                 f_i_n(:,~obj.iceFreeCauldron) = obj.f_i; % Max value for now...
            #                 f_i_n(:,obj.iceFreeCauldron)  = obj.f_i .* (c./(c+a));

        # ---- Calculate f_i_n value here -----
        f_i_n[~ice_free] = self.f_i
        f_i_n[ice_free] = self.f_i * (c[ice_free] / (c[ice_free] + self.alpha * a[ice_free]))

        return f_i_n

    def get_u_melt(self, a, c, f_i, open_idx=None, ice_free_idx=None):
        """[u_melt, v_melt] = get_u_melt(obj,a,c,f_i,openIdx,iceFreeIdx)
        IN:
          obj    = IceCauldron
          a, c  = cauldron horizontal, vertical dimensions
          f_i   = melting efficiency
        These two inputs for vectorized output:
          openIdx = optional index at which a,c are open in each cauldron
                  -> either [nCauldrons x 1] or ==size(a)
          iceFreeIdx = optional index at which a,c are ice free in each cauldron

        OUT:
          u_melt = horizontal melt rate
          v_melt = vertical melt rate

         - TODO: make radial melting area dependent on water height conditions?
               -> At least set u_melt to 0 once h_w is 0?
               -> Must balance with ice inflow, however
        """
        # TRANSLATION NOTE: MATLAB's checkVectorOrientation-based transpose handling (acTranspose)
        # is dropped here, per this project's convention of resolving vector-orientation ambiguity
        # structurally rather than with defensive orientation-checking code - see
        # validate_vector_length's TRANSLATION NOTE and CLAUDE.md's Data Structures section.
        # Callers are expected to pass a, c, f_i with consistent, already-aligned shapes.
        assert a.shape == c.shape, "u_melt: Check a_n, c_n sizes."
        assert a.shape == f_i.shape, "u_melt: Check f_i_n size."

        # Get open/ice-free booleans
        if open_idx is not None:
            if np.shape(open_idx) == a.shape:
                i_open = open_idx
            else:
                i_open = np.zeros(a.shape, dtype=bool)
                # Will fail if openIdx values are too large or shape is wrong
                for ii in range(self.n_cauldrons):
                    i_open[ii, open_idx[ii] :] = True
        else:
            # TRANSLATION NOTE: MATLAB's repmat(obj.openCauldron', size(a,2)) tiles the
            # per-cauldron open_cauldron row vector size(a,2) times along BOTH dimensions
            # (MATLAB's repmat(V,n) with scalar n tiles n times in every dimension) - genuinely
            # unclear whether this shape was intended to always match a's actual shape in the
            # general multi-cauldron case, or only happens to work out for n_cauldrons==1.
            # Translated literally via np.tile rather than guessing a "corrected" broadcast.
            i_open = np.tile(self.open_cauldron, (self.n_cauldrons, self.n_cauldrons))
        if ice_free_idx is not None:
            if np.shape(ice_free_idx) == a.shape:
                i_ice_free = ice_free_idx
            else:
                i_ice_free = np.zeros(a.shape, dtype=bool)
                # Will fail if openIdx values are too large or shape is wrong
                for ii in range(self.n_cauldrons):
                    i_ice_free[ii, ice_free_idx[ii] :] = True
        else:
            i_ice_free = np.tile(self.ice_free_cauldron, (self.n_cauldrons, self.n_cauldrons))

        u_melt = np.zeros(a.shape)
        v_melt = np.zeros(a.shape)
        # --- CLOSED CAULDRON CONDITION - ELLIPTICAL CASE ----
        if np.any(~i_open):
            # Using elliptical geometry
            # bar_u_melt = (f_i(~obj.openCauldron) .* obj.chi .* obj.Q_per_L(~obj.openCauldron)') ./ getEllipticalCylinderSA(a(iOpen),c(~obj.openCauldron,:),1,true);
            # bar_alpha  = getEllipticHeatIntensity(obj.alpha);
            try:
                bar_u_melt = (f_i * self.chi * self.params["Q_per_L"].values) / get_elliptical_cylinder_sa(a, c, 1, True)
                bar_alpha = get_elliptic_heat_intensity(self.alpha)

                u_melt[~i_open] = bar_u_melt[~i_open] / bar_alpha
                v_melt[~i_open] = self.alpha * u_melt[~i_open]
            except:  # noqa: E722 - TRANSLATION NOTE: MATLAB's bare `catch ME` (unused ME) catches everything; preserved as a bare except.
                # TRANSLATION NOTE: MATLAB's catch body is literally `faafo` - undefined-function
                # placeholder dev text (not valid MATLAB command syntax with args, just a bare
                # call), preserved as-is per CLAUDE.md's in-development-components policy. Will
                # raise NameError if ever triggered, same as MATLAB raising an "undefined
                # function" error.
                faafo()

        # --- OPEN CAULDRON CONDITION ----
        if np.any(i_open):
            # Using rectangular geometry
            try:
                # u_melt(obj.openCauldron) = (f_i(obj.openCauldron) .* obj.chi .* obj.Q_per_L(obj.openCauldron)') .* ...
                #     1./(c(obj.openCauldron) + obj.alpha .* a(obj.openCauldron));
                # v_melt(obj.openCauldron) = obj.alpha .* u_melt(obj.openCauldron);

                u_open = (f_i * self.chi * self.params["Q_per_L"].values) * 1 / (2 * (c + self.alpha * a))
                v_open = self.alpha * u_melt

                u_melt[i_open] = u_open[i_open]
                v_melt[i_open] = v_open[i_open]
            except:  # noqa: E722 - see TRANSLATION NOTE above
                faafo()

        # --- ICE FREE CAULDRON ---> stop vertical melting, geometric adjustment shifts to f_i
        # ifc = or(and(iOpen,iIceFree),c>=obj.G_n');
        ifc = i_open & i_ice_free
        if np.any(ifc):
            u_melt_if = (f_i * self.chi * self.params["Q_per_L"].values) * (1 / (2 * c))
            u_melt[ifc] = u_melt_if[ifc]
            v_melt[ifc] = 0

        return u_melt, v_melt

        #             bar_u_melt = (f_i .* obj.chi .* obj.Q_per_L') ./ getEllipticalCylinderSA(a,c,1,true);
        #
        #             bar_alpha  = getEllipticHeatIntensity(obj.alpha);
        #             u_melt(~obj.openCauldron) = bar_u_melt ./ bar_alpha;
        #             v_melt(~obj.openCauldron) = obj.alpha .* u_melt;
        #             end
        #             if any(obj.openCauldron)

        #                 bar_u_melt = (f_i .* obj.chi .* obj.Q_per_L') ./ (obj.get_CV_ice_area + );
        #             end

    def get_da_dt(self, a_n, c_n, f_i, u_ice_0):
        """[da_dt, dc_dt] = get_da_dt(gv,a_n, c_n, f_i, u_ice_0)
        Get cauldron axes growth rates
        """
        u_melt, v_melt = self.get_u_melt(a_n, c_n, f_i)

        da_dt = u_melt - u_ice_0

        # Set vertical cavity growth rate
        dc_dt = v_melt  # Ignores roof closure
        # dc_dt(c_n == gv.G_n') = 0; % TEMP: Zero once ice thickness is reached - probably breaks continuity for now

        return da_dt, dc_dt

    def get_material_heights(self, V_ice_n, V_w_n, V_p_n, V_cavity_n, a_n, c_n):
        """[H_i_n, H_w_n, H_p_n, H_cum_n] = getMaterialHeights(obj,V_ice_n,V_w_n,V_p_n,V_cavity_n,a_n,c_n)
        Return height of the various material phases
        """
        # TRANSLATION NOTE: MATLAB's checkVectorOrientation-based transpose handling and the
        # size(...)==nCauldrons-style assertions are dropped here, per this project's convention
        # (see get_u_melt's TRANSLATION NOTE and validate_vector_length) of resolving vector shape
        # consistency structurally rather than with defensive MATLAB-style orientation checks.
        # Callers are expected to pass all six arrays with matching, already-aligned shapes.
        assert (
            a_n.shape == c_n.shape == V_ice_n.shape == V_w_n.shape == V_p_n.shape == V_cavity_n.shape
        ), "Input variable dimensions do not match."

        n_steps = a_n.shape[0]

        H_i_n = np.zeros(a_n.shape)
        H_w_n = np.zeros(a_n.shape)
        H_p_n = np.zeros(a_n.shape)
        H_cum_n = np.zeros(a_n.shape)
        for nc in range(self.n_cauldrons):
            if not self.open_cauldron[nc]:
                # Closed cauldron case - elliptical
                for ti in range(n_steps):
                    # NOTE: fminbnd functions WILL prevent cumulative water and pyroclast
                    # heights greater than cavity height, c_n

                    # Pyroclast pile height
                    # TRANSLATION NOTE: MATLAB's fminbnd -> scipy.optimize.minimize_scalar with
                    # method="bounded" (similar golden-section/parabolic search; default
                    # tolerances differ slightly from MATLAB's fminbnd defaults - flagging per
                    # CLAUDE.md's ODE/solver translation notes since this is a MATLAB->SciPy
                    # optimizer substitution, not a like-for-like default).
                    def height_fun(h, ti=ti, nc=nc):
                        return abs(
                            np.arcsin(h / c_n[ti, nc])
                            + h / c_n[ti, nc] * (1 - (h / c_n[ti, nc]) ** 2) ** (1 / 2)
                            - np.pi / 2 * V_p_n[ti, nc] / V_cavity_n[ti, nc]
                        )

                    H_p_n[ti, nc] = minimize_scalar(height_fun, bounds=(0, c_n[ti, nc]), method="bounded").x

                    # Water depth
                    def height_fun(h, ti=ti, nc=nc):
                        return abs(
                            np.arcsin(h / c_n[ti, nc])
                            + h / c_n[ti, nc] * (1 - (h / c_n[ti, nc]) ** 2) ** (1 / 2)
                            - np.pi / 2 * (V_p_n[ti, nc] + V_w_n[ti, nc]) / V_cavity_n[ti, nc]
                        )

                    H_w_n[ti, nc] = minimize_scalar(height_fun, bounds=(0, c_n[ti, nc]), method="bounded").x - H_p_n[ti, nc]

                    H_i_n[ti, nc] = self.params["G_n"].values[nc] - c_n[ti, nc]
                assert np.all(H_w_n >= 0), "Water height returned negative, conditions check required."
            else:
                # Open cauldron case - rectangular

                # NOTE: these functions WILL NOT bound the cumulative heights
                # -> !! Currently the open cauldron transition may produce a
                # discontinuous change in heights and therefore water pressures
                A_c = 2 * a_n[:, nc] * self.params["L_n"].values[nc]
                H_p_n[:, nc] = V_p_n[:, nc] / A_c
                H_w_n[:, nc] = V_w_n[:, nc] / A_c
                H_i_n[:, nc] = V_ice_n[:, nc] / A_c

        H_cum_n = H_p_n + H_w_n + H_i_n
        return H_i_n, H_w_n, H_p_n, H_cum_n

    # ------------------------------- PLUME ---------------------------------
    # Translated from @IceCauldron/getPlumeFluxes.m
    def get_plume_fluxes(self, Ze, Q, a, L):
        """[h_m, r_C, qw0, qwC, qs0, qsC] = getPlumeFluxes(obj, Ze, Q, a, L)
        Given matched vector of Ze, Q, return 6 plume flux/length values:
        Input:
          -> Ze : (m) Current water depth above vent
          -> Q  : (kg/s) Current mass eruption rate at vent
          -> a  : (m) Cauldron radius/half-width
          -> L  : (m) Cauldron (fissure segment) length

        Output:
          -> Qs0  : (kg/s) Mass flux of pyroclasts into the subaerial plume base
          -> QsC  : (kg/s) Mass flux of pyroclasts at collapse point or plume max height
          -> Qw0  : (kg/s) Mass flux of water into the subaerial plume base
          -> QwC  : (kg/s) Mas flux of water at collapse point or plume max height
          -> Hmax : (m) Max plume height
          -> Rmax : (m) Plume radius at collapse point or max height

         - Todo - vectorize as needed - at least for multiple cauldrons?
         - Todo - load random forest metadata in IceCauldron, raise warning if Ze/Q go
               outside training bounds
        """
        # TRANSLATION NOTE: a and L are unused in the active MATLAB code path (only referenced in
        # the fully commented-out footprint/collapse-fraction section below) - preserved as
        # unused parameters to match the MATLAB signature.
        #
        # TRANSLATION NOTE: per project direction, this whole random-forest emulator path is an
        # in-development/test mode, not finalized physics - see gvpy/plume.py's module docstring.
        # The disable_plume_flux "off" branch (zeros) is the reliable path for now.
        #
        # TRANSLATION NOTE: MATLAB's two branches produce different shapes (zeros(6,1) column vs.
        # predictForest's (1,6) row), but linear-indexes both identically as fluxes(1)..fluxes(6).
        # Unified here to a flat length-6 array for both branches, consistent with this project's
        # convention of not carrying MATLAB's row/column distinction into Python.
        if self.disable_plume_flux:
            fluxes = np.zeros(6)
        else:
            log_q = np.log10(Q)
            fluxes = predict_forest(self.plume_emulator, np.array([[Ze, log_q]]))[0]

        h_m = fluxes[0]  # Plume height (m)
        r_c = fluxes[1]  # Plume radius at collapse point or buoyancy (m)
        qw0 = fluxes[2]  # Water flux at plume source (kg/s)
        qwc = fluxes[3]  # Water flux at collapse point/buoyancy (kg/s)
        qs0 = fluxes[4]  # Pyroclast flux at plume source (kg/s)
        qsc = fluxes[5]  # Pyroclast flux at collapse point/buoyancy (kg/s)

        # Todo - figure out how the collapse radius part and net flux needs to work
        # again...?

        # Rectangular cauldron footprint
        # if strcmp(obj.geometry, 'cylinder')
        #     cpoly = polyshape([-a -a a a],[-L L L -L]/2);
        # else
        #     error('Only cylindrical cauldron geometries are implemented for plume flux')
        # end

        # Circular plume footprint
        # xc = 0;
        # yc = 0;
        # n = 100;
        # theta = (0:n-1)*(2*pi/n);
        # x = xc + r_C.*cos(theta);
        # y = yc + r_C.*sin(theta);
        # ppoly = polyshape(x,y);
        #
        # clps_frac = area(intersect(cpoly,ppoly))./area(ppoly);
        # clps_frac(clps_frac>1) = 1;
        #
        # qw_clps = -qwC.*clps_frac;
        # qs_clps = -qsC.*clps_frac;
        #
        # qwNet = qw0 + qw_clps;
        # qsNet = qs0 + qs_clps;
        # Qv    = fluxSurfaces.Qv;
        # hwv   = h_m;

        return h_m, r_c, qw0, qwc, qs0, qsc

    # --------------------------- SUPRAGLACIAL DRAINAGE -----------------------
    # Translated from @IceCauldron/get_q_s.m
    def get_q_s(self, V_fluid_n, V_cavity_n, dV_fluid_n, dV_cavity_n):
        """q_s_n = get_q_s(obj,V_fluid_n, V_cavity_n, dV_fluid_n, dV_cavity_n)
        Get supra-glacial drainage term
        """
        # TRANSLATION NOTE: MATLAB's checkVectorOrientation-based transpose handling is dropped
        # here, per this project's convention (see get_u_melt's TRANSLATION NOTE) of resolving
        # vector shape consistency structurally rather than with defensive orientation-checking.
        # Also dropped: MATLAB's separate `any(size(V_fluid_n)==obj.nCauldrons)` check (whose
        # error message oddly refers to "a_n", not V_fluid_n - a copy-paste artifact in the
        # MATLAB source, preserved only as this note since the check itself is being simplified
        # away) - combined into one shape-equality assert across all four inputs, consistent with
        # get_material_heights' equivalent simplification.
        if self.supraglacial_drainage_mode == "overflow":
            assert V_fluid_n.shape == V_cavity_n.shape == dV_fluid_n.shape == dV_cavity_n.shape, (
                "Input variable dimensions do not match."
            )

            q_s_n = np.zeros(V_fluid_n.shape)
            # overflow = (V_fluid_n >= obj.G_n & dV_fluid_n > dV_cavity_n);
            # overflow = ( abs(V_fluid_n - V_cavity_n) <= V_tol  & dV_fluid_n > dV_cavity_n); % alternative Volume tolerance overflow
            overflow = (V_fluid_n + dV_fluid_n) / (V_cavity_n + dV_cavity_n) > 1  # Alternative volume+volume change overflow
            if np.any(overflow):
                try:
                    q_s_n[overflow] = dV_fluid_n[overflow] - dV_cavity_n[overflow]
                except:  # noqa: E722 - TRANSLATION NOTE: MATLAB's bare `catch ME` (unused ME) catches everything; preserved as a bare except.
                    # TRANSLATION NOTE: MATLAB's catch body is literally `faafo fix it` -
                    # undefined-function placeholder dev text using MATLAB *command syntax*
                    # (equivalent to calling faafo('fix','it')). Preserved as-is per CLAUDE.md's
                    # in-development-components policy; translated to explicit call syntax since
                    # MATLAB command syntax isn't valid Python (a mechanical syntax conversion,
                    # not a logic change) - would raise NameError if ever triggered, matching
                    # MATLAB's "undefined function" error. In practice this is unreachable under
                    # normal, already-shape-validated inputs (confirmed via test) - the try body
                    # is a plain, non-failing numpy assignment - so this is dead/defensive code,
                    # not a bug that fires in ordinary use.
                    faafo("fix", "it")

        elif self.supraglacial_drainage_mode == "off":  # Drainage disabled
            q_s_n = np.zeros(V_fluid_n.shape)

        return q_s_n

    # ------------------- SUBGLACIAL DRAINAGE (IN DEVELOPMENT) ---------------
    # Translated from @IceCauldron/getFluxuralRigidity.m, get_L_lambda.m, get_U_tip.m,
    # get_del_h.m, drainageDensity.m, getDrainageFluxes.m, solveDeltaP.m
    #
    # FLAG: per project direction, the ENTIRE subglacial drainage component below is explicitly
    # in-development/test-only, not finalized physics (see gvpy/subglacial.py's module
    # docstring). get_drainage_fluxes - the actual subglacial + inter-cauldron drainage flux
    # calculation - is not implemented in the MATLAB source at all (dummy zeros only), and
    # solve_delta_p only implements one regime, explicitly erroring otherwise. Bugs/incompleteness
    # here are preserved as-is per CLAUDE.md's in-development-components policy, not fixed.
    def get_fluxural_rigidity(self, h):
        """D = getFluxuralRigidity(obj,h)
        [SUBGLACIAL DRAINAGE - IN DEVELOPMENT]
        Get flexural rigidity of a beam of glacier ice with thickness h (m).
        INPUT:
          h  : scalar or array of arbitrary size. Beam thickness in meters
        OUTPUT:
          D       : (Pa . m^3) Flexural rigidity
          E_prime : reduced young's modulus, E / (1 - nu^2), where nu is
                      Poisson's ratio
        """
        # TRANSLATION NOTE: "Fluxural" is a typo in the original MATLAB function name (should be
        # "Flexural") - preserved verbatim in this method's name for traceability to the MATLAB
        # source, per this project's precedent of preserving source typos (see
        # get_elevation_profiles' "traingular" docstring typo).
        # E_prime =  obj.E ./ (1 - obj.nu.^2);
        D = h**3 / 12 * self.E_prime

        return D

    def get_l_lambda(self, x, smooth_transition=False):
        """[l, l_hi_ratio, l_lambda] = get_L_lambda(obj,x,smoothTransition)
        [SUBGLACIAL DRAINAGE - IN DEVELOPMENT]
        Estimate Pressure wave dimension as 1/4*(ice flexural wavelength), as a
        function of down-glacier coordinate x
          Input:
              x           : position in meters along glacier
              smoothTransition : T/F - apply linear combination over a transition
                              length to allow l to smoothly approach l_lambda
          Output:
              l           : actual flexural length scale
                            (equal to x when x < l_lambda, otherwise = l_lambda)
              l_hi_ratio  : l_lambda/(ice thickness) ratio
              l_lambda    : actual 1/4 flexural wavelength at position x
        """
        # TRANSLATION NOTE: MATLAB's checkVectorOrientation-based transpose handling is dropped
        # here, per this project's convention (see get_u_melt's TRANSLATION NOTE) of resolving
        # vector shape consistency structurally rather than with defensive orientation-checking.
        x = np.asarray(x, dtype=float)
        nx = x.size  # unused directly below - only referenced by the commented-out dead code at the end

        # The fraction of the critical transition l_lambda over which to apply smoothing
        default_smooth_fraction = 1 / 5

        h_i, _, _, _, _ = self.get_elevation_profiles(x)

        D = self.get_fluxural_rigidity(h_i)

        # A) One quarter flexural wavelength formulation
        # l_lambda = 2*pi.*(D./(obj.Constants.g*obj.Constants.rho_ice)).^(1/4) / 4;
        # B) Flexural parameter formulation
        l_lambda = (4 * D / (self.CONSTANTS.G * self.CONSTANTS.RHO_ICE)) ** (1 / 4)

        # True l flexural length (equal to x when x < l_lambda)
        l = np.minimum(x, l_lambda)

        if smooth_transition:
            x_crit = get_x_equals_l_lambda(self)

            l = smooth_union(x, l_lambda, x_crit * default_smooth_fraction, "cubic")

            # Quick plot the effects of different smoother kernels
            #     figure
            #     plot(x,x,'k',x,l_lambda,'k--')
            #     set(gca,'ColorOrderIndex',1);
            #     hold on
            #     kernels = {'quadratic','cubic','quartic','circular'};
            #     for ii = 1:length(kernels)
            #         l_smooth{ii} = smoothUnion(x,l_lambda,x_crit/5,kernels{ii});
            #         plot(x,l_smooth{ii})
            #     end
            #     legend([{'x','l_\lambda'} kernels])

        l_hi_ratio = l / h_i

        # ---- ALL THIS COMMENTED STUFF vv ----
        #   Because l_lambda is non-linear in h, thought it might be important to
        #   average it over a certain range of x. Turns out, not really. Therefore
        #   use simpler function above.
        #  Below:
        #  -> The rigidity function in h_i^3 is averaged over +/- 1 h_i,
        #    bounded from x = 0 to l_G
        #
        # dz_i = obj.h_i_0 - obj.h_i_f;
        # dz_b = obj.z_b_0;
        #
        # ddz = dz_i - dz_b;
        # ddzl = ddz ./ obj.l_G;
        # h0 = obj.h_i_0;
        #
        # min_x = max([repmat(0,size(x)) x-h_i],[],2);
        # max_x = min([ x+h_i repmat(obj.l_G,size(x))],[],2);
        #
        # % flexural constant
        # [~,E_prime] = obj.getFluxuralRigidity(h_i);
        # C = pi/2 .* (E_prime / (12*obj.Constants.g*obj.Constants.rho_ice) ).^(1/4);
        #
        # % int_h = @(x) h0.^3.*x + ddzl.^3.*x.^4/4 + h0.*ddzl.^2.*x.^3 + 3/2.*h0.^2.*ddzl.*x.^2;
        #
        # % l_lambda = 1./(max_x-min_x) .* C .* (int_h(max_x) - int_h(min_x));
        #
        # int_fun = @(xx) C .* obj.getElevationProfiles(xx).^(3/4);
        #
        # l_lambda = zeros(size(x));
        # for ii = 1:nx
        #     l_lambda(ii) = 1./(max_x(ii)-min_x(ii)) .* integral(int_fun,min_x(ii),max_x(ii));
        # end

        return l, l_hi_ratio, l_lambda

    def get_u_tip(self, dP, rho, h_i, l):
        """U_tip = get_U_tip(obj,dP,rho,h_i,l)
        [SUBGLACIAL DRAINAGE - IN DEVELOPMENT]
        This function for crack tip velocity is from Tsai & Rice (2012) (Journal
        of Applied Mechanics). It is known to be valid for L_hi_ratio <= 5. It is
        unknown how valid the relation is above that point.
         dP   = driving fluid overpressure
         rho  = fluid density
         x    = position of crack tip along glacier profile
        """
        # Deciding for now to input h_i,l, instead of calculating internally, for
        # more flexible usage
        #     [h_i, ~, ~, ~,~] = obj.getElevationProfiles(x);
        #     [~,E_prime] = obj.getFluxuralRigidity(h_i);
        #
        #     [l, l_hi_ratio] = get_L_lambda(obj,x);
        l_hi_ratio = l / h_i

        phi = 5.13 + 0.64 * l_hi_ratio + 0.94 * l_hi_ratio**2

        # TRANSLATION NOTE: this recomputes the Nikuradse roughness height locally rather than
        # reusing self.k_nikuradse (already computed identically in __attrs_post_init__) -
        # preserved as-is, matching the MATLAB source's own redundant computation.
        k = (self.manning_roughness / 0.038) ** 6  # Nikuradse roughness height

        U_tip = phi * rho ** (-1 / 2) * (self.E_prime) ** (-2 / 3) * (l / k) ** (1 / 6) * (dP) ** (7 / 6)
        U_tip = np.where(np.asarray(dP) < 0, 0.0, U_tip)

        return U_tip

    def get_del_h(self, dP, h_i, l):
        """del_h = get_del_h(obj,dP,h_i,l)
        [SUBGLACIAL DRAINAGE - IN DEVELOPMENT]
        Get vertical ice deflection due to elastic pressure.
          - dP  : peak overpressure (occurs at position l)
          - h_i : average ice thickness above the crack
          - l   : length of crack that is *elastically* supported
        """
        l_arr = np.asarray(l)
        dP_arr = np.asarray(dP)
        assert l_arr.shape == dP_arr.shape or l_arr.size == 1 or dP_arr.size == 1, (
            "Input x and delta P sizes must match, or at least one must be scalar."
        )

        # Model 1 --- fixed elastic strain + buoyancy
        # epsilon = 1e-3;
        # del_h = l_lambda .* gv.epsilon;

        # Model 2 --- buoyant force --> 0, elastic strain only
        #  --> following Tsai and Rice (2012)
        #     [l, l_hi_ratio]     = obj.get_L_lambda(x);
        #     [~,E_prime]         = obj.getFluxuralRigidity(obj.getElevationProfiles(x));
        l_hi_ratio = l / h_i

        w_hat = 1.72 + 0.89 * l_hi_ratio**2
        del_h = w_hat * dP * l / self.E_prime

        return del_h

    def drainage_density(self):
        """[rho_f,phi_w,phi_p,chi_f] = drainageDensity(obj)
        [SUBGLACIAL DRAINAGE - IN DEVELOPMENT]
        Get density and mass fractions of draining fluid, depending on choice of
        model:
        Par = parameter struct for glaciovolcano script. Must include field
        "flood_density_model".
        FLOOD_DENSITY_MODEL CHOICES:
          1 = pure water only
          2 = fixed fraction of tephra in outflow to give constant drainage
              density somewhat greater than water (requires specified input value)
          3 = same as 2, but includes a fixed fraction of lithics (~10 vol% in Tomasson 1996)
          4 = [DEFAULT] tephra fraction proportional to melt/magma production rates
              (ie gamma, Lambda - this is ~ the assumption of Tomasson 1996)
          5 = same as 3, but includes a fixed fraction of lithics (~10 vol% in Tomasson 1996)
          6 = outflow tephra fraction proportional to relative amounts in
              cauldron? Easy to get unrealistic for low water remaining
          ? = phi proportional to pressure, flow rate, fill volumes?

        OUT:
          rho_f   = builk density of flood fluid
          phi_w   = volume fraction of water in discharge
          phi_t   = volume fraction of tephra in discharge
          chi_f   = water/tephra mass ratio of floodwater
        """
        C = self.CONSTANTS

        # rho_r = 2500; % Lithics density

        # TRANSLATION NOTE: MATLAB's switch has no otherwise/default case, and case 6 is only a
        # comment stub (`% case 6 % Based on pressure/discharge/proportions in cauldron?`), not
        # an implemented branch - preserved as-is per this component's in-development policy.
        # Selecting flood_density_model==6 (a value the attrs validator otherwise allows, 1-6)
        # falls through with phi_p/nr_nt undefined, raising UnboundLocalError below - matching
        # MATLAB's equivalent undefined-variable failure.
        if self.flood_density_model == 1:
            rho_f = C.RHO_W
            phi_w = 1
            phi_p = 0
            chi_f = np.nan

        elif self.flood_density_model == 2:
            phi_p = self.fixed_phi_p
            phi_r = 0.0
            phi_w = 1 - phi_p - phi_r
            chi_f = phi_w * C.RHO_W / (phi_p * self.rho_p)
            nr_nt = 0

        elif self.flood_density_model == 3:
            phi_p = self.fixed_phi_p
            phi_r = 0.1
            phi_w = 1 - phi_p - phi_r
            chi_f = phi_w * C.RHO_W / (phi_p * self.rho_p)
            nr_nt = phi_r * self.rho_r / (phi_p * self.rho_p)

        elif self.flood_density_model == 4:
            # chi_f = obj.gamma .* obj.cm .* obj.deltaT ./ obj.Li;
            chi_f = self.chi * C.RHO_W  # Idealized water/tephra mass ratio
            nr_nt = 0

        elif self.flood_density_model == 5:
            # chi_f = obj.gamma .* obj.cm .* obj.deltaT ./ obj.Li;
            chi_f = self.chi * C.RHO_W
            nr_nt = 1 / 3

        # case 6 % Based on pressure/discharge/proportions in cauldron?

        if self.flood_density_model != 1:
            nt = (1 + chi_f + nr_nt) ** (-1)
            nw = nt * chi_f
            nr = nt * nr_nt

            rho_f = (nw / C.RHO_W + nt / self.rho_p + nr / self.rho_r) ** (-1)
            phi_w = rho_f * nw / C.RHO_W
            phi_p = rho_f * nt / self.rho_p

        return rho_f, phi_w, phi_p, chi_f

    def get_drainage_fluxes(self, phi_w, P_w, H_w_n, H_cum_n):
        """[q_d_n,q_c_n] = getDrainageFluxes(obj,phi_w,P_w,H_w_n,H_cum_n)
        [SUBGLACIAL DRAINAGE - IN DEVELOPMENT]
        Calculate the three flood discharge terms
        """
        # TRANSLATION NOTE: not implemented in the MATLAB source either - this only returns dummy
        # zero fluxes and ignores all of its inputs (phi_w, P_w, H_w_n, H_cum_n). The actual
        # subglacial + inter-cauldron drainage flux calculation remains to be written.
        q_d_n = np.zeros(self.n_cauldrons)  # dummy for now
        q_c_n = np.zeros(self.n_cauldrons)  # dummy for now

        # q_s_n = zeros(gv.nCauldrons,1); % dummy for now

        return q_d_n, q_c_n

    def solve_delta_p(self, x, P_0, rho_f, u_prev=None, smooth_l_lambda=True):
        """[Delta_P,u_1,h_lambda] = solveDeltaP(obj,x,P_0,rho_f,u_prev,smooth_l_lambda)
        [SUBGLACIAL DRAINAGE - IN DEVELOPMENT]
        Iterative solution to find a deltaP with frictional pressure loss that
        satisfies momentum and volume conservation equations.

        -> assuming known, scalar P0 and l for now
        """
        # global smooth_l_lambda % maybe temporary
        default_u_guess = 1  # (m/s) Initial guess for an average flow velocity that satisfies pressure losses

        # Get ice thickness, elevations, static pressures at beginning and end of profile
        # [h_i, z_i, z_b, theta_i, theta_b] = obj.getElevationProfiles([0 x]);
        h_i, _, _, _, _ = self.get_elevation_profiles(x)
        # P_i = gv.Constants.rho_ice .* gv.Constants.g .* h_i;

        # Flexural wavelength parameter
        l, _, l_lambda = self.get_l_lambda(x, smooth_l_lambda)

        # TRANSLATION NOTE: MATLAB's `opts` struct (TolX, FunValCheck) configures fminsearch;
        # FunValCheck (error on NaN/Inf objective values) has no direct scipy.optimize equivalent
        # and is dropped here - flagging per CLAUDE.md's ODE/solver translation notes.
        # opts.PlotFcns = @optimplotfval;

        if x <= l_lambda:
            # Case of early propagating crack tip:
            # -> u_1 = u_tip
            # -> (bar_h_s / h_s(L_s)) ~ 1 (average crack opening ~ inlet opening)
            # -> u_0 = 4 * u_tip * (bar_h_s / h_s(L_s)) = 4 * u_tip
            # ->
            if u_prev is None:
                u_max = self.get_u_tip(P_0, rho_f, h_i, l)  # noqa: F841 - unused, matches MATLAB (u_max feeds the commented-out fminbnd call only)
                u_guess = default_u_guess
            else:
                u_max = u_prev * 2  # noqa: F841 - see note above
                u_guess = u_prev

            # dP_max = Delta_P_from_Bernoulli_1(obj,0,P_0,rho_f,x,10); % The 10 needs a fix...
            def dP_fun(u):
                return objective_delta_p_1(self, u, P_0, rho_f, x)

            # u_1 = fminbnd(dP_fun,0,u_max,opts);
            # TRANSLATION NOTE: MATLAB's fminsearch (unconstrained Nelder-Mead) ->
            # scipy.optimize.minimize's Nelder-Mead method, with TolX -> xatol. Flagging per
            # CLAUDE.md's ODE/solver translation notes since this is a MATLAB->SciPy optimizer
            # substitution, not a like-for-like default.
            result = minimize(dP_fun, u_guess, method="Nelder-Mead", options={"xatol": 1e-5})
            u_1 = result.x[0]

            Delta_P = delta_p_from_u_tip(self, u_1, rho_f, x)
            h_lambda = self.get_del_h(x, Delta_P, h_i, l)
            # Delta_P_check = Delta_P_from_Bernoulli_1(obj,u_1,P_0,rho_f,x, del_h);

            # ---- TEMPORARY QC PLOTS FOR TESTING/DEV ----
            # TRANSLATION NOTE: this diagnostic plotting block (MATLAB figure/plot/get(gca,...)
            # calls, plus a typo referencing `hi` instead of `h_i`) is dev/debug-only scaffolding,
            # not core physics, and depends on MATLAB-specific plotting (ismembertol, rgba2rgb)
            # that has no counterpart yet in this Python translation's (untouched) plotting
            # utilities. Left untranslated rather than porting to matplotlib as a tangent to this
            # pass - flagging per CLAUDE.md's "flag if out of scope for this pass" guidance rather
            # than silently dropping.
            #     test_x_values = [10 100 500 1000 1500 2000 2030 2050 2300];
            # test_x_values = [1800, 2030, 2300]  # values around which dP solutions diverge sharply
            # x_check = ismembertol(test_x_values,x,1e-3)
            # if any(x_check):
            #     u_test = linspace(0,200,2001)
            #     dP_test1 = Delta_P_from_U_tip(obj,u_test,rho_f,x)
            #     del_h_test   = obj.get_del_h(x,dP_test1,hi,l)
            #     dP_test2 = Delta_P_from_Bernoulli_1(obj,u_test,P_0,rho_f,x, del_h_test)
            #     co = get(gca,'ColorOrder')
            #     co = [co; rgba2rgb(co,0.5)]
            #     figure(8)
            #     hold on
            #     plot(u_test,dP_test1,'Color',co(x_check,:),'DisplayName',sprintf('U_{tip}, x = %.0f',x))
            #     plot(u_test,dP_test2,'--','Color',co(x_check,:),'DisplayName',sprintf('B, x = %.0f',x))
        else:
            raise NotImplementedError("Regime of x > l_lambda not yet implemented.")

        return Delta_P, u_1, h_lambda

    # ---- HELPER FUNCTIONS to get lists of solver-calculated variables ---
    def get_solution_indices(self):
        """Quick helper function to tell which indices in ode solver
        solutions y correspond to which variables, since number of
        cauldrons can vary.

        indices = dict with keys corresponding to variable names and
              values equal to a slice giving the range of each variable
              in the ode solution vector
        n = number of solution variables
        """
        # TRANSLATION NOTE: MATLAB built a struct of 1-based index ranges; translated here to a
        # dict of 0-based slice objects instead, matching the concatenation order used in
        # get_initial_conditions.
        n_vector_vars = len(self.vector_solution_vars)
        n_scalar_vars = len(self.scalar_solution_vars)

        n_vector_indices = self.n_cauldrons * n_vector_vars

        indices = {}
        for fi, name in enumerate(self.vector_solution_vars):
            start = fi * self.n_cauldrons
            indices[name] = slice(start, start + self.n_cauldrons)
        for fi, name in enumerate(self.scalar_solution_vars):
            indices[name] = n_vector_indices + fi
        # - todo: Any non-per-cauldron fields should be added in here -

        n = n_vector_indices + n_scalar_vars

        return indices, n

    def solution_vars(self):
        """Direct solver vars."""
        return self.vector_solution_vars + self.scalar_solution_vars

    def derived_vars(self):
        """Secondary vars derived from solution."""
        return self.vector_derived_vars + self.scalar_derived_vars

    def time_series_vars(self):
        # Secondary vars derived from solution
        # TRANSLATION NOTE: this comment is carried over verbatim from the MATLAB source, but
        # appears to be a copy-paste of derived_vars' comment above it - this function actually
        # returns ALL time series vars (solution + derived), not just the secondary/derived
        # ones. Flagging the apparent doc/comment mismatch rather than silently correcting it.
        return self.solution_vars() + self.derived_vars()

    def vector_vars(self):
        """Per-cauldron variables."""
        return self.vector_solution_vars + self.vector_derived_vars

    def scalar_vars(self):
        """Single variables."""
        return self.scalar_solution_vars + self.scalar_derived_vars

    # ---------------------------------------------------------------------
    def get_elevation_profiles(self, x):
        """Return linear (average) elevation profiles and slope angles of triangular
        glacier surface and bedrock.

        OUTPUT:
          H_i     : glacier thickness (m)
          z_i     : ice surface profile (m)
          z_b     : bedrock surface profile (m)
          theta_i : ice surface slope (radians)
          theta_b : bedrock surface slope (radians)
        """
        x = np.asarray(x, dtype=float)

        z_i_0 = self.z_b_0 + self.h_i_0

        z_i_f = self.h_i_f  # Final ice elevation (usually 0)

        # Bedrock elevation profile
        z_b = self.z_b_0 + (0 - self.z_b_0) * x / self.l_G

        # Glacier elevation profile
        z_i = z_i_0 + (z_i_f - z_i_0) * x / self.l_G

        # Get ice thickness profile
        h_i = z_i - z_b

        theta_i = np.arctan((z_i_0 - z_i_f) / self.l_G)
        theta_b = np.arctan(self.z_b_0 / self.l_G)

        return h_i, z_i, z_b, theta_i, theta_b

    @staticmethod
    def validate_vector_length(prop, n_cauldrons, func_name, var_name):
        """Quick function to validate vector property lengths."""
        # TRANSLATION NOTE: MATLAB's validateattributes also checked that prop was numeric or
        # logical, and required a [1 x n_cauldrons] row-vector shape specifically. Translated
        # here to a generic 1-D shape check only (no row/column distinction, no type-class
        # check), consistent with dropping checkVectorOrientation - see CLAUDE.md's Data
        # Structures section on resolving MATLAB's orientation ambiguity structurally.
        prop = np.asarray(prop)
        if prop.shape != (n_cauldrons,):
            raise ValueError(f"{func_name}: {var_name} must have shape ({n_cauldrons},), got {prop.shape}.")

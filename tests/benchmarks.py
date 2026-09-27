"""Analytical benchmarks for comparison against completed gv_main runs.

Translated from MATLAB source: benchmarking/benchmark_functions/cylindricalCauldronGrowthBenchmark.m

TRANSLATION NOTE: the MATLAB function's (dat, benchmark, makePlots) -> result_table signature is
restructured per project direction: the core calculation (MATLAB struct C) is the
CylindricalCauldronGrowthBenchmark class, the RMSE check against `benchmark` (errorBenchmark.m) is
replaced by the placeholder test_cylindrical_cauldron_growth, and the makePlots section by the
placeholder CylindricalCauldronGrowthBenchmark.plot.
"""

from dataclasses import dataclass, field
from typing import Callable

import numpy as np
import xarray as xr

from gvpy.gv_main import GVResult
from gvpy.utilities.constants import ThermoConstants


@dataclass
class CylindricalCauldronGrowthBenchmark:
    """Simple cauldron growth at a constant melt rate (MATLAB struct C), built from a completed run.

    dat  : completed gv_main run - supplies the IceCauldron (dat.gv) and output times
    data : xarray Dataset - benchmark curves with dims (time, cauldron), per-cauldron values with
           dim (cauldron)
    """

    dat:             GVResult
    test_name:       str = field(init=False)
    f_ice_area:      Callable = field(init=False)
    f_cavity_radius: Callable = field(init=False)
    V_f_a:           Callable = field(init=False)
    data:            xr.Dataset = field(init=False)

    # Reset fields to make for simple cauldron growth
    # TRANSLATION NOTE: MATLAB's deprecated allow_ice_inflow = false -> ice_inflow_mode = "off".
    # alpha and ice_inflow_mode aren't used in the calc; kept for comparison against the run's
    # IceCauldron settings.
    alpha:           float = field(init=False, default=1.0)
    ice_inflow_mode: str = field(init=False, default="off")
    a0:              float = field(init=False, default=0.0)

    def __post_init__(self):
        gv = self.dat.gv
        # TRANSLATION NOTE: MATLAB's %i on a (possibly non-integer) alpha -> :g.
        self.test_name = f"simpleCylinder: alpha = {gv.alpha:g}"

        nv = self.dat.data.sizes["time"]  # noqa: F841

        L_n = gv.params["L_n"]  # (m)
        G_n = gv.params["G_n"]  # (m)

        # TRANSLATION NOTE: MATLAB's switch has no otherwise case. Python: IceCauldron's "spheroid" 
        # geometry option raises NotImplementedError at construction anyway.
        if gv.geometry == "cylinder":
            self.f_ice_area = lambda a: np.pi * a * L_n  # Ignores ends
            self.f_cavity_radius = lambda V: (V * (2 / (np.pi * L_n))) ** (1 / 2)
            self.V_f_a = lambda a: np.pi * a**2 * L_n / 2

            V_max = 2 * G_n * G_n * L_n  # (becomes a rectangular box)

        elif gv.geometry == "spheroid":
            # Not currently implemented - same situation as in ice_cauldron
            raise NotImplementedError("Spheroid geometry is not currently implemented")  
            c0 = np.nan  # noqa: F841
            l0 = np.nan  # noqa: F841
            plate_model = 6  # noqa: F841

            self.f_ice_area = lambda a: 2 / 3 * np.pi * a**2
            self.f_cavity_radius = lambda V: (V * (3 / (2 * np.pi))) ** (1 / 3)
            self.V_f_a = lambda a: 2 / 3 * np.pi * a**3

            V_max = np.pi * G_n**2 * G_n  # (Becomes a vertical cylinder)

        else:
            raise ValueError(f"Unsupported geometry: {gv.geometry}")

        # %% Evolution calcs
        S0 = self.f_ice_area(self.a0)
        V0 = self.V_f_a(self.a0)

        V_melt_max = V_max - V0
        melt_vol_per_s = gv.f_i * gv.params["melt_vol_per_s"]  # (m3/s)

        t_melt = V_melt_max / melt_vol_per_s  # (s)
        # TRANSLATION NOTE: MATLAB's repmat(dat.t,[1 gv.nCauldrons]) is replaced by the run's time
        # coordinate - xarray broadcasts it against the per-cauldron values below.
        t_vec = self.dat.data["time"]
        #     C.t_vec      = zeros(nv,gv.nCauldrons);
        #     for ci = 1:gv.nCauldrons
        #         C.t_vec(:,ci) = linspace(0,C.t_melt(ci),nv);
        #     end
        #     C.t_vec      = linspace(0,C.t_melt,nv);

        # Total, unbounded ice volume melted at constant melt rate
        V_of_t = lambda t: V0 + melt_vol_per_s * t  # noqa: E731

        V_water = V_of_t(t_vec) * ThermoConstants.RHO_ICE / ThermoConstants.RHO_W
        V_cavity = V_of_t(t_vec)
        r = self.f_cavity_radius(V_cavity)
        # TRANSLATION NOTE: with a0 = 0, r = 0 at t = 0, so dr_dt is inf there (NumPy also emits a
        # divide-by-zero RuntimeWarning; MATLAB is silent).
        dr_dt = melt_vol_per_s / self.f_ice_area(r)

        V_ice = (2 * r * G_n * L_n) - V_cavity

        # Stop each cauldron's curves (NaN) once it melts through, c_n == G_n (c = r for a half
        # cylinder) - the cylindrical benchmark is no longer meaningful past that point.
        is_closed = r <= G_n
        V_water, V_cavity, r, dr_dt, V_ice = (x.where(is_closed) for x in (V_water, V_cavity, r, dr_dt, V_ice))

        self.data = xr.Dataset(
            {
                "S0": S0,
                "V0": V0,
                "V_max": V_max,
                "V_melt_max": V_melt_max,
                "melt_vol_per_s": melt_vol_per_s,
                "t_melt": t_melt,
                "V_water": V_water,
                "V_cavity": V_cavity,
                "r": r,
                "dr_dt": dr_dt,
                "V_ice": V_ice,
            }
        ).transpose("time", "cauldron")

    def plot(self):
        """Plot benchmark curves over the simulation output (MATLAB's makePlots section)."""
        # TODO: translate once compatible with the run-results dashboard plotting.
        pass


def test_cylindrical_cauldron_growth(dat: GVResult):
    """Compare a completed run against CylindricalCauldronGrowthBenchmark (MATLAB: errorBenchmark on
    cavity volume)."""
    # TODO: define test criteria (MATLAB RMSE thresholds deliberately not reused).
    pass

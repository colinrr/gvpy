"""Right-hand side of the model's ODE system, integrated by gv_main.py.

Translated from MATLAB source: glaciovolcano.m
"""


import numpy as np

from .ice_cauldron import IceCauldron
from .utilities.solver_helpers import assign_integrated_values


def glaciovolcano(t: float, y: np.ndarray, gv: IceCauldron, return_par: bool = False) -> np.ndarray | tuple:
    """[ydot, par] = glaciovolcano(t, y, gv)
    glaciovolcano defines the system of equations (ideally in
    nondimensional form )
    to solve for the accumulation of meltwater during a
    subglacial eruption.

    Returns ydot, or (ydot, par) when return_par=True (MATLAB's nargout > 1).
    """
    # TRANSLATION NOTE: MATLAB's optional second output `par` (nargout > 1) becomes the
    # `return_par` flag, and `par` is a dict. Per-cauldron values are read from `gv.params`
    # (broadcast to arrays) rather than the raw gv.G_n / gv.L_n fields, which are stored as
    # passed (scalar or list). `global step_ct` is dropped: it is only referenced by the
    # commented-out plotting hack below, never by active code. (t is unused - autonomous system.)

    # TODO: Build input/setup for u_ice_0
    # TODO: implement non-dimensionalization
    #

    #     if t>0.000001
    #         step_ct = 0;
    #     end
    # =========== CAULDRON EVOLUTION =============

    # Get integrated values from vector a
    C = gv.CONSTANTS
    int_vars = assign_integrated_values(y, gv)
    a_n = int_vars["a_n"]
    c_n = int_vars["c_n"]
    V_ice_n = int_vars["V_ice_n"]
    V_w_n = int_vars["V_w_n"]
    V_p_n = int_vars["V_p_n"]
    V_cavity_n = gv.get_cv_control_volume(a_n) - V_ice_n

    # Get heat transfer partitioning efficiencies
    f_i_n = gv.get_f_i(a_n, c_n)
    u_ice_0, u_ice_bar = gv.get_u_ice()

    # Cavity horizontal and vertical growth rate(s)
    da_dt, dc_dt = gv.get_da_dt(a_n, c_n, f_i_n, u_ice_0)

    # =============== DRAINAGE ===============
    # Material fill heights in column to set water pressure
    H_i_n, H_w_n, H_p_n, H_cum_n = gv.get_material_heights(V_ice_n, V_w_n, V_p_n, V_cavity_n, a_n, c_n)

    # Drainage fluid pressure and density in cauldron
    P_w = 0  # dummy for now
    phi_w = 0.4
    phi_p = 1 - phi_w

    # Density model for drainage outflow
    #     [rho_f,phi_w,phi_p,chi_f] = gv.drainageDensity;
    _, phi_w, phi_p, _ = gv.drainage_density()
    q_d_n, q_c_n = gv.get_drainage_fluxes(phi_w, P_w, H_w_n, H_cum_n)

    # =============== CONS. OF VOLUME EQN'S ===============

    # Ice in control volume(s)
    # A_cv_n = gv.G_n .* gv.L_n .* 2; % Cylindrical geom
    Q_n = gv.params["Q_n"].values
    dV_i_melt_dt = f_i_n * gv.chi * Q_n
    dV_ice_dt = -dV_i_melt_dt + gv.get_cv_ice_area() * (da_dt - u_ice_bar)

    # Water in cavity(s)
    dV_w_dt = C.RHO_ICE / C.RHO_W * dV_i_melt_dt - phi_w * (q_d_n + q_c_n)

    # Pyroclast volume in cavity
    dV_p_dt = Q_n / gv.rho_p - phi_p * (q_d_n + q_c_n)

    # ---> Iterate here to compensate volume flux with supra-glacial
    # drainage
    #     q_s_n = gv.get_q_s(H_cum_n, (dV_ice_dt+dV_w_dt+dV_p_dt - gv.get_CV_ice_area .* (da_dt)) );
    q_s_n = gv.get_q_s((V_w_n + V_p_n), V_cavity_n, (dV_w_dt + dV_p_dt), dV_i_melt_dt)
    dV_w_dt = dV_w_dt - phi_w * q_s_n
    dV_p_dt = dV_p_dt - phi_p * q_s_n

    # Temp hack to plot steps
    #     if and(t>200, gv.supraglacial_drainage_mode ~=  'off')
    # %         step_ct = step_ct + 1;
    #         figure(4)
    # %         plot(step_ct,V_ice_n,'.')
    #         plot(t,V_w_n./gv.get_CV_control_volume(a_n),'bo')
    #         hold on
    # %         figure(4)
    #         plot(t,V_p_n./gv.get_CV_control_volume(a_n),'r.')
    #         plot(t,V_ice_n./gv.get_CV_control_volume(a_n),'c.')
    #         hold on
    #     end

    ydot = np.concatenate([da_dt, dc_dt, dV_ice_dt, dV_w_dt, dV_p_dt])

    # Quick check to make sure we're keeping things up-to-date
    #   - maybe not strictly necessary since we separately check IC's
    _, ns = gv.get_solution_indices()
    try:
        assert len(ydot) == ns, "Length of solution vector does not match IceCauldron solution vars."
    except:  # noqa: E722 - TRANSLATION NOTE: MATLAB's bare `catch ME` (unused ME) catches everything; preserved as a bare except.
        # TRANSLATION NOTE: MATLAB's catch body is literally `faafo` - undefined-function
        # placeholder dev text, preserved as-is per CLAUDE.md's in-development-components policy
        # (would raise NameError if the assert ever failed).
        faafo()

    # Additional fields for reconstruction at end of simulation
    if return_par:
        par = {
            "dV_w_dt": dV_w_dt,
            "dV_p_dt": dV_p_dt,
            "dV_ice_dt": dV_ice_dt,
            "dV_cavity_dt": dV_i_melt_dt,
            "q_s_n": q_s_n,
            "phi_w": phi_w,
            "phi_p": phi_p,
            "phi_r": 1 - phi_p - phi_w,
            # Flux terms, exposed for plotting
            "Q_n": Q_n,
            "u_ice_bar": u_ice_bar,
            "q_d_n": q_d_n,
            "q_c_n": q_c_n,
        }
        return ydot, par

    return ydot

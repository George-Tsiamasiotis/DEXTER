import numpy as np
import pytest
import dexter as dex
from math import isfinite

from semver import Version

LCFS = dex.MagneticFlux.Toroidal(0.05)


def test_flute():
    mode = dex.FluteMode(1e-4, LCFS, 3, 2, 0)
    _test_mode_base(mode)
    assert mode.lcfs == dex.MagneticFlux.Toroidal(0.05)
    assert mode.phase == 0
    assert mode.epsilon == 1e-4


def test_nc(nc_flute_mode: dex.NcFluteMode):
    _test_mode_base(nc_flute_mode)
    assert nc_flute_mode.interp_type == "Cubic"
    assert isinstance(nc_flute_mode.path, str)
    assert isinstance(nc_flute_mode.netcdf_version, Version)
    assert nc_flute_mode.phase_method == "Interpolation"
    assert nc_flute_mode.analytical_threshold_index == 3
    assert nc_flute_mode.phase_average is None
    assert isinstance(nc_flute_mode.psi_array, np.ndarray)
    assert isinstance(nc_flute_mode.psip_array, np.ndarray)
    assert isinstance(nc_flute_mode.alpha_array, np.ndarray)
    assert isinstance(nc_flute_mode.phase_array, np.ndarray)


def _test_mode_base(mode: dex.ModeObject):

    mode.__repr__()
    mode.__str__()

    assert mode.machine_type in ["Numerical", "Analytical"]
    assert mode.psi_state in ["Good", "Bad"]
    assert mode.psip_state in ["Good", "Bad"]

    assert mode.m == 3
    assert mode.n == 2

    flux = 0.02
    fluxes = np.linspace(0.01, 0.04, 10)
    theta = zeta = t = 1.2
    thetas = zetas = ts = np.linspace(0.01, np.pi, 10)

    try:

        mode.eval_amplitude(psi=flux, theta=theta, zeta=zeta, t=t)
        mode.eval_amplitude(psi=fluxes, theta=thetas, zeta=zetas, t=ts)
        mode.eval_amplitude(psip=flux, theta=theta, zeta=zeta, t=t)
        mode.eval_amplitude(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        mode.eval_phase(psi=flux, theta=theta, zeta=zeta, t=t)
        mode.eval_phase(psi=fluxes, theta=thetas, zeta=zetas, t=ts)
        mode.eval_phase(psip=flux, theta=theta, zeta=zeta, t=t)
        mode.eval_phase(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        mode.eval_m(psi=flux, theta=theta, zeta=zeta, t=t)
        mode.eval_m(psi=fluxes, theta=thetas, zeta=zetas, t=ts)
        mode.eval_m(psip=flux, theta=theta, zeta=zeta, t=t)
        mode.eval_m(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        mode.eval_deriv_flux(psi=flux, theta=theta, zeta=zeta, t=t)
        mode.eval_deriv_flux(psi=fluxes, theta=thetas, zeta=zetas, t=ts)
        mode.eval_deriv_flux(psip=flux, theta=theta, zeta=zeta, t=t)
        mode.eval_deriv_flux(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        mode.eval_deriv_theta(psi=flux, theta=theta, zeta=zeta, t=t)
        mode.eval_deriv_theta(psi=fluxes, theta=thetas, zeta=zetas, t=ts)
        mode.eval_deriv_theta(psip=flux, theta=theta, zeta=zeta, t=t)
        mode.eval_deriv_theta(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        mode.eval_deriv_zeta(psi=flux, theta=theta, zeta=zeta, t=t)
        mode.eval_deriv_zeta(psi=fluxes, theta=thetas, zeta=zetas, t=ts)
        mode.eval_deriv_zeta(psip=flux, theta=theta, zeta=zeta, t=t)
        mode.eval_deriv_zeta(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

        mode.eval_deriv_t(psi=flux, theta=theta, zeta=zeta, t=t)
        mode.eval_deriv_t(psi=fluxes, theta=thetas, zeta=zetas, t=ts)
        mode.eval_deriv_t(psip=flux, theta=theta, zeta=zeta, t=t)
        mode.eval_deriv_t(psip=fluxes, theta=thetas, zeta=zetas, t=ts)

    except Exception as e:
        if not "[D] EvalError" in str(e):
            raise RuntimeError(
                f"only testing the vectorized functions here (error: {e})"
            )

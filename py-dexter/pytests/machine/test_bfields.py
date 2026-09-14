import numpy as np
import pytest
import dexter as dex
from math import isfinite

from semver import Version


def test_lar():
    bfield = dex.LarBfield()
    _test_bfield_base(bfield)


def test_nc(nc_bfield: dex.NcBfield):
    _test_bfield_base(nc_bfield)
    assert isinstance(nc_bfield.path, str)
    assert isinstance(nc_bfield.netcdf_version, Version)
    assert nc_bfield.interp_type == "Bicubic"
    assert isfinite(nc_bfield.baxis)
    assert isfinite(nc_bfield.padding)
    assert isinstance(nc_bfield.padding_theta, float)
    assert isinstance(nc_bfield.psi_array, np.ndarray)
    assert isinstance(nc_bfield.psip_array, np.ndarray)
    assert isinstance(nc_bfield.theta_array, np.ndarray)
    assert isinstance(nc_bfield.theta_array_padded, np.ndarray)
    b_array = nc_bfield.b_array
    assert isinstance(b_array, np.ndarray)
    assert b_array.shape == nc_bfield.shape
    b_array_padded = nc_bfield.b_array_padded
    assert isinstance(b_array_padded, np.ndarray)
    assert b_array_padded.shape == nc_bfield.shape_padded


def _test_bfield_base(bfield: dex.BfieldObject):

    bfield.__repr__()
    bfield.__str__()

    assert bfield.machine_type in ["Numerical", "Analytical"]
    assert bfield.psi_state in ["Good", "Bad"]
    assert bfield.psip_state in ["Good", "Bad"]

    flux = 0.02
    fluxes = np.linspace(0, 0.04, 10)
    theta = 1.2
    thetas = np.linspace(0, np.pi, 10)

    try:

        bfield.eval_b(psi=flux, theta=theta)
        bfield.eval_b(psi=fluxes, theta=thetas)
        bfield.eval_b(psip=flux, theta=theta)
        bfield.eval_b(psip=fluxes, theta=thetas)

        bfield.eval_deriv_flux(psi=flux, theta=theta)
        bfield.eval_deriv_flux(psi=fluxes, theta=thetas)
        bfield.eval_deriv_flux(psip=flux, theta=theta)
        bfield.eval_deriv_flux(psip=fluxes, theta=thetas)

        bfield.eval_deriv_theta(psi=flux, theta=theta)
        bfield.eval_deriv_theta(psi=fluxes, theta=thetas)
        bfield.eval_deriv_theta(psip=flux, theta=theta)
        bfield.eval_deriv_theta(psip=fluxes, theta=thetas)

    except Exception as e:
        if not "[D] EvalError" in str(e):
            raise RuntimeError(
                f"only testing the vectorized functions here (error: {e})"
            )

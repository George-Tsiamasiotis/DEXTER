import numpy as np
import pytest
import dexter as dex
from math import isfinite

from semver import Version


def test_lar():
    current = dex.LarCurrent()
    _test_current_base(current)


def test_nc(nc_current: dex.NcCurrent):
    _test_current_base(nc_current)
    assert nc_current.interp_type == "Cubic"
    assert isinstance(nc_current.path, str)
    assert isinstance(nc_current.netcdf_version, Version)
    assert isinstance(nc_current.psi_array, np.ndarray)
    assert isinstance(nc_current.psip_array, np.ndarray)
    assert isinstance(nc_current.i_array, np.ndarray)
    assert isinstance(nc_current.g_array, np.ndarray)


def _test_current_base(current: dex.CurrentObject):

    current.__repr__()
    current.__str__()

    assert current.machine_type in ["Numerical", "Analytical"]
    assert current.psi_state in ["Good", "Bad"]
    assert current.psip_state in ["Good", "Bad"]

    methods = [
        current.eval_g,
        current.eval_i,
        current.eval_g_deriv,
        current.eval_i_deriv,
    ]

    flux = 0.02
    fluxes = np.linspace(0, 0.04, 10)

    try:

        current.eval_g(psi=flux)
        current.eval_g(psi=fluxes)
        current.eval_g(psip=flux)
        current.eval_g(psip=fluxes)

        current.eval_i(psi=flux)
        current.eval_i(psi=fluxes)
        current.eval_i(psip=flux)
        current.eval_i(psip=fluxes)

        current.eval_g_deriv(psi=flux)
        current.eval_g_deriv(psi=fluxes)
        current.eval_g_deriv(psip=flux)
        current.eval_g_deriv(psip=fluxes)

        current.eval_i_deriv(psi=flux)
        current.eval_i_deriv(psi=fluxes)
        current.eval_i_deriv(psip=flux)
        current.eval_i_deriv(psip=fluxes)

    except Exception as e:
        if not "[D] EvalError" in str(e):
            raise RuntimeError(
                f"only testing the vectorized functions here (error: {e})"
            )

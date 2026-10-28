import numpy as np
import dexter as dex

from semver import Version


def test_lar():
    current = dex.LarCurrent()
    _test_current_base(current)


def test_nc(nc_current: dex.NcCurrent):
    _test_current_base(nc_current)
    assert nc_current.interp_type == "Cubic"
    assert isinstance(nc_current.path, str)
    assert isinstance(nc_current.netcdf_version, Version)
    assert isinstance(nc_current.psit_array, np.ndarray)
    assert isinstance(nc_current.psip_array, np.ndarray)
    assert isinstance(nc_current.i_array, np.ndarray)
    assert isinstance(nc_current.g_array, np.ndarray)


def _test_current_base(current: dex.CurrentObject):

    current.__repr__()
    current.__str__()

    assert current.machine_type in ["Numerical", "Analytical"]
    assert current.psit_state in ["Good", "Bad"]
    assert current.psip_state in ["Good", "Bad"]

    flux = 0.02
    fluxes = np.linspace(0, 0.04, 10)

    try:

        current.eval_g(psit=flux)
        current.eval_g(psit=fluxes)
        current.eval_g(psip=flux)
        current.eval_g(psip=fluxes)

        current.eval_i(psit=flux)
        current.eval_i(psit=fluxes)
        current.eval_i(psip=flux)
        current.eval_i(psip=fluxes)

        current.eval_g_deriv(psit=flux)
        current.eval_g_deriv(psit=fluxes)
        current.eval_g_deriv(psip=flux)
        current.eval_g_deriv(psip=fluxes)

        current.eval_i_deriv(psit=flux)
        current.eval_i_deriv(psit=fluxes)
        current.eval_i_deriv(psip=flux)
        current.eval_i_deriv(psip=fluxes)

        current.eval_g_deriv2(psit=flux)
        current.eval_g_deriv2(psit=fluxes)
        current.eval_g_deriv2(psip=flux)
        current.eval_g_deriv2(psip=fluxes)

        current.eval_i_deriv2(psit=flux)
        current.eval_i_deriv2(psit=fluxes)
        current.eval_i_deriv2(psip=flux)
        current.eval_i_deriv2(psip=fluxes)

    except Exception as e:
        if not "[D] EvalError" in str(e):
            raise RuntimeError(
                f"only testing the vectorized functions here (error: {e})"
            )

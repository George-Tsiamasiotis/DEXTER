import numpy as np
import pytest
import dexter as dex
from math import isfinite

from semver import Version


def test_lar():
    geometry = dex.LarGeometry(2, 1.75, 0.5)
    _test_geometry_base(geometry)


def test_nc(nc_geometry: dex.NcGeometry):
    nc_geometry.__repr__()
    nc_geometry.__str__()
    _test_geometry_base(nc_geometry)
    assert isinstance(nc_geometry.path, str)
    assert isinstance(nc_geometry.netcdf_version, Version)
    assert nc_geometry.interp1d_type == "Cubic"
    assert nc_geometry.interp2d_type == "Bicubic"
    assert nc_geometry.psi_last is not None
    assert nc_geometry.psip_last is not None
    assert isfinite(nc_geometry.psi_last.value)
    assert isfinite(nc_geometry.psip_last.value)
    assert isinstance(nc_geometry.psi_array, np.ndarray)
    assert isinstance(nc_geometry.psip_array, np.ndarray)
    assert isinstance(nc_geometry.theta_array, np.ndarray)
    assert isinstance(nc_geometry.r_array, np.ndarray)
    rlab_array = nc_geometry.rlab_array
    assert isinstance(rlab_array, np.ndarray)
    assert rlab_array.shape == nc_geometry.shape
    zlab_array = nc_geometry.zlab_array
    assert isinstance(zlab_array, np.ndarray)
    assert zlab_array.shape == nc_geometry.shape
    jacobian_array = nc_geometry.jacobian_array
    assert isinstance(jacobian_array, np.ndarray)
    assert jacobian_array.shape == nc_geometry.shape


def _test_geometry_base(geometry: dex.GeometryObject):

    geometry.__repr__()
    geometry.__str__()

    assert geometry.machine_type in ["Numerical", "Analytical"]
    assert geometry.psi_state in ["Good", "Bad"]
    assert geometry.psip_state in ["Good", "Bad"]

    assert isfinite(geometry.baxis)
    assert isfinite(geometry.raxis)
    assert isfinite(geometry.zaxis)
    assert isfinite(geometry.rgeo)
    assert isfinite(geometry.rlast)
    assert isinstance(geometry.rlab_last, np.ndarray)
    assert isinstance(geometry.zlab_last, np.ndarray)

    r = 0.01
    rs = np.linspace(0, 0.02, 10)
    flux = 0.02
    fluxes = np.linspace(0, 0.04, 10)
    theta = 1.2
    thetas = np.linspace(0, np.pi, 10)
    try:

        geometry.eval_r(psi=flux)
        geometry.eval_r(psi=fluxes)
        geometry.eval_r(psip=flux)
        geometry.eval_r(psip=fluxes)

        geometry.eval_psi_of_r(r=r)
        geometry.eval_psi_of_r(r=rs)
        geometry.eval_psip_of_r(r=r)
        geometry.eval_psip_of_r(r=rs)

        geometry.eval_rlab(psi=flux, theta=theta)
        geometry.eval_rlab(psi=fluxes, theta=thetas)
        geometry.eval_rlab(psip=flux, theta=theta)
        geometry.eval_rlab(psip=fluxes, theta=thetas)

        geometry.eval_zlab(psi=flux, theta=theta)
        geometry.eval_zlab(psi=fluxes, theta=thetas)
        geometry.eval_zlab(psip=flux, theta=theta)
        geometry.eval_zlab(psip=fluxes, theta=thetas)

        geometry.eval_jacobian(psi=flux, theta=theta)
        geometry.eval_jacobian(psi=fluxes, theta=thetas)
        geometry.eval_jacobian(psip=flux, theta=theta)
        geometry.eval_jacobian(psip=fluxes, theta=thetas)

    except Exception as e:
        if not "[D] EvalError" in str(e):
            raise RuntimeError(
                f"only testing the vectorized functions here (error: {e})"
            )

import pytest
import numpy as np
from math import isfinite
from dexter.equilibrium import _TOROIDAL_TEST_NETCDF_PATH, _POLOIDAL_TEST_NETCDF_PATH
from dexter import LarGeometry, NcGeometry


def test_lar_geometry():
    geometry = LarGeometry(baxis=2, raxis=1.75, rwall=0.5)
    assert geometry.equilibrium_type == "Analytical"
    assert geometry.baxis == 2
    assert geometry.raxis == 1.75
    assert geometry.rwall == 0.5
    assert isfinite(geometry.psi_wall)
    assert np.all(np.isfinite(geometry.rlab_wall))
    assert np.all(np.isfinite(geometry.zlab_wall))
    assert isfinite(geometry.r_of_psi(0.01))
    assert isfinite(geometry.psi_of_r(0.01))
    assert isfinite(geometry.rlab_of_psi(0.01, 3.14))
    assert isfinite(geometry.zlab_of_psi(0.01, 3.14))
    with pytest.raises(BaseException):
        geometry.r_of_psip(0.01)
    with pytest.raises(BaseException):
        geometry.psip_of_r(0.01)
    with pytest.raises(BaseException):
        geometry.rlab_of_psip(0.01, 3.14)
    with pytest.raises(BaseException):
        geometry.zlab_of_psip(0.01, 3.14)
    with pytest.raises(BaseException):
        geometry.jacobian_of_psi(0.01, 3.14)
    with pytest.raises(BaseException):
        geometry.jacobian_of_psip(0.01, 3.14)
    assert isinstance(geometry.__str__(), str)
    assert isinstance(geometry.__repr__(), str)


def test_nc_geometry_getters(nc_geometry: NcGeometry):
    assert isinstance(nc_geometry.path, str)
    assert isinstance(nc_geometry.netcdf_version, str)
    assert nc_geometry.equilibrium_type == "Numerical"
    assert isinstance(nc_geometry.interp1d_type, str)
    assert isinstance(nc_geometry.interp1d_type, str)
    assert isfinite(nc_geometry.baxis)
    assert isfinite(nc_geometry.raxis)
    assert isfinite(nc_geometry.zaxis)
    assert isfinite(nc_geometry.rgeo)
    assert isfinite(nc_geometry.rwall)
    assert len(nc_geometry.shape) == 2
    assert isfinite(nc_geometry.psi_wall)
    assert isfinite(nc_geometry.psip_wall)
    assert nc_geometry.psi_state == "Good"
    assert nc_geometry.psip_state == "Good"
    assert nc_geometry.psi_array.ndim == 1
    assert nc_geometry.psip_array.ndim == 1
    assert nc_geometry.theta_array.ndim == 1
    assert nc_geometry.r_array.ndim == 1
    assert nc_geometry.rlab_array.ndim == 2
    assert nc_geometry.zlab_array.ndim == 2
    assert nc_geometry.jacobian_array.ndim == 2
    assert isinstance(nc_geometry.__str__(), str)
    assert isinstance(nc_geometry.__repr__(), str)


def test_nc_geometry_eval(nc_geometry: NcGeometry):
    _test_geometry_evals(nc_geometry)


def test_toroidal_nc_geometry():
    nc_geometry = NcGeometry(_TOROIDAL_TEST_NETCDF_PATH, "Steffen", "Bicubic")
    assert nc_geometry.psi_state == "Good"
    assert nc_geometry.psip_state == "Bad"
    assert isfinite(nc_geometry.psip_of_psi(0.01))
    assert isfinite(nc_geometry.r_of_psi(0.01))
    assert isfinite(nc_geometry.psi_of_r(0.01))
    assert isfinite(nc_geometry.psip_of_r(0.01))
    assert isfinite(nc_geometry.rlab_of_psi(0.01, 3.14))
    assert isfinite(nc_geometry.zlab_of_psi(0.01, 3.14))
    assert isfinite(nc_geometry.jacobian_of_psi(0.01, 3.14))
    with pytest.raises(BaseException):
        nc_geometry.psi_of_psip(0.01)
    with pytest.raises(BaseException):
        nc_geometry.r_of_psip(0.01)
    with pytest.raises(BaseException):
        nc_geometry.rlab_of_psip(0.01, 3.14)
    with pytest.raises(BaseException):
        nc_geometry.zlab_of_psip(0.01, 3.14)
    with pytest.raises(BaseException):
        nc_geometry.jacobian_of_psip(0.01, 3.14)


def test_poloidal_nc_geometry():
    nc_geometry = NcGeometry(_POLOIDAL_TEST_NETCDF_PATH, "Steffen", "Bicubic")
    assert nc_geometry.psi_state == "Bad"
    assert nc_geometry.psip_state == "Good"
    assert isfinite(nc_geometry.psi_of_psip(0.01))
    assert isfinite(nc_geometry.r_of_psip(0.01))
    assert isfinite(nc_geometry.psi_of_r(0.01))
    assert isfinite(nc_geometry.psip_of_r(0.01))
    assert isfinite(nc_geometry.rlab_of_psip(0.01, 3.14))
    assert isfinite(nc_geometry.zlab_of_psip(0.01, 3.14))
    assert isfinite(nc_geometry.jacobian_of_psip(0.01, 3.14))
    with pytest.raises(BaseException):
        nc_geometry.psip_of_psi(0.01)
    with pytest.raises(BaseException):
        nc_geometry.r_of_psi(0.01)
    with pytest.raises(BaseException):
        nc_geometry.rlab_of_psi(0.01, 3.14)
    with pytest.raises(BaseException):
        nc_geometry.zlab_of_psi(0.01, 3.14)
    with pytest.raises(BaseException):
        nc_geometry.jacobian_of_psi(0.01, 3.14)


def _test_geometry_evals(geometry: NcGeometry):
    r = 0.02
    psi = 0.01
    psip = 0.015
    theta = 3.14
    assert isfinite(geometry.psip_of_psi(psi))
    assert isfinite(geometry.psi_of_psip(psip))
    assert isfinite(geometry.r_of_psi(psi))
    assert isfinite(geometry.r_of_psip(psip))
    assert isfinite(geometry.psi_of_r(r))
    assert isfinite(geometry.psip_of_r(r))
    assert isfinite(geometry.rlab_of_psi(psi, theta))
    assert isfinite(geometry.rlab_of_psip(psip, theta))
    assert isfinite(geometry.zlab_of_psi(psi, theta))
    assert isfinite(geometry.zlab_of_psip(psip, theta))
    assert isfinite(geometry.jacobian_of_psi(psi, theta))
    assert isfinite(geometry.jacobian_of_psip(psip, theta))

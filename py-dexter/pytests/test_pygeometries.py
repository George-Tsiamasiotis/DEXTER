import pytest
from dexter.equilibrium import _TOROIDAL_TEST_NETCDF_PATH, _POLOIDAL_TEST_NETCDF_PATH
from dexter import LarGeometry, NcGeometry


def test_lar_geometry():
    geometry = LarGeometry(baxis=2, raxis=1.75, rwall=0.5)
    assert isinstance(geometry.__str__(), str)
    assert isinstance(geometry.__repr__(), str)
    assert isinstance(geometry.r_of_psi(0.01), float)
    assert isinstance(geometry.psi_of_r(0.01), float)
    assert isinstance(geometry.rlab_of_psi(0.01, 3.14), float)
    assert isinstance(geometry.zlab_of_psi(0.01, 3.14), float)
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


def test_nc_geometry_getters(nc_geometry: NcGeometry):
    assert isinstance(nc_geometry.__str__(), str)
    assert isinstance(nc_geometry.__repr__(), str)
    assert isinstance(nc_geometry.path, str)
    assert isinstance(nc_geometry.netcdf_version, str)
    assert nc_geometry.equilibrium_type == "Numerical"
    assert isinstance(nc_geometry.interp1d_type, str)
    assert isinstance(nc_geometry.interp1d_type, str)
    assert isinstance(nc_geometry.baxis, float)
    assert isinstance(nc_geometry.raxis, float)
    assert isinstance(nc_geometry.zaxis, float)
    assert isinstance(nc_geometry.rgeo, float)
    assert isinstance(nc_geometry.rwall, float)
    assert len(nc_geometry.shape) == 2
    assert isinstance(nc_geometry.psi_wall, float)
    assert isinstance(nc_geometry.psip_wall, float)
    assert nc_geometry.psi_state == "Good"
    assert nc_geometry.psip_state == "Good"
    assert nc_geometry.psi_array.ndim == 1
    assert nc_geometry.psip_array.ndim == 1
    assert nc_geometry.theta_array.ndim == 1
    assert nc_geometry.r_array.ndim == 1
    assert nc_geometry.rlab_array.ndim == 2
    assert nc_geometry.zlab_array.ndim == 2
    assert nc_geometry.jacobian_array.ndim == 2


def test_nc_geometry_eval(nc_geometry: NcGeometry):
    _test_geometry_evals(nc_geometry)


def test_toroidal_nc_geometry():
    nc_geometry = NcGeometry(_TOROIDAL_TEST_NETCDF_PATH, "Steffen", "Bicubic")
    assert nc_geometry.psi_state == "Good"
    assert nc_geometry.psip_state == "Bad"
    assert isinstance(nc_geometry.psip_of_psi(0.01), float)
    assert isinstance(nc_geometry.r_of_psi(0.01), float)
    assert isinstance(nc_geometry.psi_of_r(0.01), float)
    assert isinstance(nc_geometry.psip_of_r(0.01), float)
    assert isinstance(nc_geometry.rlab_of_psi(0.01, 3.14), float)
    assert isinstance(nc_geometry.zlab_of_psi(0.01, 3.14), float)
    assert isinstance(nc_geometry.jacobian_of_psi(0.01, 3.14), float)
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
    assert isinstance(nc_geometry.psi_of_psip(0.01), float)
    assert isinstance(nc_geometry.r_of_psip(0.01), float)
    assert isinstance(nc_geometry.psi_of_r(0.01), float)
    assert isinstance(nc_geometry.psip_of_r(0.01), float)
    assert isinstance(nc_geometry.rlab_of_psip(0.01, 3.14), float)
    assert isinstance(nc_geometry.zlab_of_psip(0.01, 3.14), float)
    assert isinstance(nc_geometry.jacobian_of_psip(0.01, 3.14), float)
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
    assert isinstance(geometry.psip_of_psi(psi), float)
    assert isinstance(geometry.psi_of_psip(psip), float)
    assert isinstance(geometry.r_of_psi(psi), float)
    assert isinstance(geometry.r_of_psip(psip), float)
    assert isinstance(geometry.psi_of_r(r), float)
    assert isinstance(geometry.psip_of_r(r), float)
    assert isinstance(geometry.rlab_of_psi(psi, theta), float)
    assert isinstance(geometry.rlab_of_psip(psip, theta), float)
    assert isinstance(geometry.zlab_of_psi(psi, theta), float)
    assert isinstance(geometry.zlab_of_psip(psip, theta), float)
    assert isinstance(geometry.jacobian_of_psi(psi, theta), float)
    assert isinstance(geometry.jacobian_of_psip(psip, theta), float)

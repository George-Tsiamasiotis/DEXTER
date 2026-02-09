import pytest
from dexter.equilibrium import _TOROIDAL_TEST_NETCDF_PATH, _POLOIDAL_TEST_NETCDF_PATH
from dexter import LarBfield, NcBfield


def test_lar_bfield():
    bfield = LarBfield()
    assert isinstance(bfield.__str__(), str)
    assert isinstance(bfield.__repr__(), str)
    assert bfield.equilibrium_type == "Analytical"
    assert isinstance(bfield.b_of_psi(0.01, 3.14), float)
    assert isinstance(bfield.db_dpsi(0.01, 3.14), float)
    assert isinstance(bfield.db_of_psi_dtheta(0.01, 3.14), float)
    with pytest.raises(BaseException):
        bfield.b_of_psip(0.01, 3.14)
    with pytest.raises(BaseException):
        bfield.db_dpsip(0.01, 3.14)
    with pytest.raises(BaseException):
        bfield.db_of_psip_dtheta(0.01, 3.14)


def test_nc_bfield_getters(nc_bfield: NcBfield):
    assert isinstance(nc_bfield.__str__(), str)
    assert isinstance(nc_bfield.__repr__(), str)
    assert isinstance(nc_bfield.path, str)
    assert isinstance(nc_bfield.netcdf_version, str)
    assert nc_bfield.equilibrium_type == "Numerical"
    assert isinstance(nc_bfield.interp_type, str)
    assert isinstance(nc_bfield.baxis, float)
    assert len(nc_bfield.shape) == 2
    assert isinstance(nc_bfield.psi_wall, float)
    assert isinstance(nc_bfield.psip_wall, float)
    assert nc_bfield.psi_state == "Good"
    assert nc_bfield.psip_state == "Good"
    assert nc_bfield.psi_array.ndim == 1
    assert nc_bfield.psip_array.ndim == 1
    assert nc_bfield.theta_array.ndim == 1
    assert nc_bfield.b_array.ndim == 2


def test_nc_bfield_eval(nc_bfield: NcBfield):
    psi = 0.01
    psip = 0.015
    theta = 3.14
    assert isinstance(nc_bfield.b_of_psi(psi, theta), float)
    assert isinstance(nc_bfield.b_of_psip(psip, theta), float)
    assert isinstance(nc_bfield.db_dpsi(psi, theta), float)
    assert isinstance(nc_bfield.db_dpsip(psip, theta), float)
    assert isinstance(nc_bfield.db_of_psi_dtheta(psi, theta), float)
    assert isinstance(nc_bfield.db_of_psip_dtheta(psip, theta), float)


def test_toroidal_nc_bfield():
    nc_bfield = NcBfield(_TOROIDAL_TEST_NETCDF_PATH, "Bicubic")
    assert nc_bfield.psi_state == "Good"
    assert nc_bfield.psip_state == "Bad"
    assert isinstance(nc_bfield.b_of_psi(0.01, 3.14), float)
    assert isinstance(nc_bfield.db_dpsi(0.01, 3.14), float)
    assert isinstance(nc_bfield.db_of_psi_dtheta(0.01, 3.14), float)
    with pytest.raises(BaseException):
        nc_bfield.b_of_psip(0.01, 3.14)
    with pytest.raises(BaseException):
        nc_bfield.db_dpsip(0.01, 3.14)
    with pytest.raises(BaseException):
        nc_bfield.db_of_psip_dtheta(0.01, 3.14)


def test_poloidal_nc_bfield():
    nc_bfield = NcBfield(_POLOIDAL_TEST_NETCDF_PATH, "Bicubic")
    assert nc_bfield.psi_state == "Bad"
    assert nc_bfield.psip_state == "Good"
    assert isinstance(nc_bfield.b_of_psip(0.01, 3.14), float)
    assert isinstance(nc_bfield.db_dpsip(0.01, 3.14), float)
    assert isinstance(nc_bfield.db_of_psip_dtheta(0.01, 3.14), float)
    with pytest.raises(BaseException):
        nc_bfield.b_of_psi(0.01, 3.14)
    with pytest.raises(BaseException):
        nc_bfield.db_dpsi(0.01, 3.14)
    with pytest.raises(BaseException):
        nc_bfield.db_of_psi_dtheta(0.01, 3.14)

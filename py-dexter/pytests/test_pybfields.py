import pytest
from dexter.equilibrium import _TOROIDAL_TEST_NETCDF_PATH, _POLOIDAL_TEST_NETCDF_PATH
from dexter import LarBfield, NcBfield


def test_lar_bfield():
    bfield = LarBfield()
    assert isinstance(bfield.__str__(), str)
    assert isinstance(bfield.__repr__(), str)
    assert bfield.equilibrium_type == "Analytical"
    args = (0.01, 3.14)
    assert isinstance(bfield.b_of_psi(*args), float)
    assert isinstance(bfield.db_dpsi(*args), float)
    assert isinstance(bfield.db_of_psi_dtheta(*args), float)
    with pytest.raises(BaseException):
        bfield.b_of_psip(*args)
    with pytest.raises(BaseException):
        bfield.db_dpsip(*args)
    with pytest.raises(BaseException):
        bfield.db_of_psip_dtheta(*args)


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
    args = (0.01, 3.14)
    assert isinstance(nc_bfield.b_of_psi(*args), float)
    assert isinstance(nc_bfield.b_of_psip(*args), float)
    assert isinstance(nc_bfield.db_dpsi(*args), float)
    assert isinstance(nc_bfield.db_dpsip(*args), float)
    assert isinstance(nc_bfield.db_of_psi_dtheta(*args), float)
    assert isinstance(nc_bfield.db_of_psip_dtheta(*args), float)


def test_toroidal_nc_bfield():
    nc_bfield = NcBfield(_TOROIDAL_TEST_NETCDF_PATH, "Bicubic")
    assert nc_bfield.psi_state == "Good"
    assert nc_bfield.psip_state == "Bad"
    args = (0.01, 3.14)
    assert isinstance(nc_bfield.b_of_psi(*args), float)
    assert isinstance(nc_bfield.db_dpsi(*args), float)
    assert isinstance(nc_bfield.db_of_psi_dtheta(*args), float)
    with pytest.raises(BaseException):
        nc_bfield.b_of_psip(*args)
    with pytest.raises(BaseException):
        nc_bfield.db_dpsip(*args)
    with pytest.raises(BaseException):
        nc_bfield.db_of_psip_dtheta(*args)


def test_poloidal_nc_bfield():
    nc_bfield = NcBfield(_POLOIDAL_TEST_NETCDF_PATH, "Bicubic")
    assert nc_bfield.psi_state == "Bad"
    assert nc_bfield.psip_state == "Good"
    args = (0.01, 3.14)
    assert isinstance(nc_bfield.b_of_psip(*args), float)
    assert isinstance(nc_bfield.db_dpsip(*args), float)
    assert isinstance(nc_bfield.db_of_psip_dtheta(*args), float)
    with pytest.raises(BaseException):
        nc_bfield.b_of_psi(*args)
    with pytest.raises(BaseException):
        nc_bfield.db_dpsi(*args)
    with pytest.raises(BaseException):
        nc_bfield.db_of_psi_dtheta(*args)

import pytest
from dexter import NcCurrent, _TOROIDAL_TEST_NETCDF_PATH, _POLOIDAL_TEST_NETCDF_PATH


def test_nc_current_getters(nc_current: NcCurrent):
    assert isinstance(nc_current.path, str)
    assert isinstance(nc_current.typ, str)
    assert isinstance(nc_current.psi_wall, float)
    assert isinstance(nc_current.psip_wall, float)
    assert nc_current.psi_state == "Good"
    assert nc_current.psip_state == "Good"
    assert nc_current.psi_array.ndim == 1
    assert nc_current.psip_array.ndim == 1
    assert nc_current.g_array.ndim == 1
    assert nc_current.i_array.ndim == 1
    assert len(nc_current) == len(nc_current.psi_array)
    assert isinstance(nc_current.__str__(), str)
    assert isinstance(nc_current.__repr__(), str)


def test_nc_current_eval(nc_current: NcCurrent):
    psi = 0.01
    psip = 0.015
    assert isinstance(nc_current.g_of_psi(psi), float)
    assert isinstance(nc_current.g_of_psip(psip), float)
    assert isinstance(nc_current.i_of_psi(psi), float)
    assert isinstance(nc_current.i_of_psip(psip), float)
    assert isinstance(nc_current.dg_dpsi(psi), float)
    assert isinstance(nc_current.dg_dpsip(psip), float)
    assert isinstance(nc_current.di_dpsi(psi), float)
    assert isinstance(nc_current.di_dpsip(psip), float)


def test_toroidal_nc_current():
    nc_current = NcCurrent(_TOROIDAL_TEST_NETCDF_PATH, "Steffen")
    assert nc_current.psi_state == "Good"
    assert nc_current.psip_state == "Bad"
    assert isinstance(nc_current.g_of_psi(0.01), float)
    assert isinstance(nc_current.i_of_psi(0.01), float)
    assert isinstance(nc_current.dg_dpsi(0.01), float)
    assert isinstance(nc_current.di_dpsi(0.01), float)
    with pytest.raises(BaseException):
        nc_current.g_of_psip(0.01)
    with pytest.raises(BaseException):
        nc_current.i_of_psip(0.01)
    with pytest.raises(BaseException):
        nc_current.dg_dpsip(0.01)
    with pytest.raises(BaseException):
        nc_current.di_dpsip(0.01)


def test_poloidal_nc_current():
    nc_current = NcCurrent(_POLOIDAL_TEST_NETCDF_PATH, "Steffen")
    assert nc_current.psi_state == "Bad"
    assert nc_current.psip_state == "Good"
    assert isinstance(nc_current.g_of_psip(0.01), float)
    assert isinstance(nc_current.i_of_psip(0.01), float)
    assert isinstance(nc_current.dg_dpsip(0.01), float)
    assert isinstance(nc_current.di_dpsip(0.01), float)
    with pytest.raises(BaseException):
        nc_current.g_of_psi(0.01)
    with pytest.raises(BaseException):
        nc_current.i_of_psi(0.01)
    with pytest.raises(BaseException):
        nc_current.dg_dpsi(0.01)
    with pytest.raises(BaseException):
        nc_current.di_dpsi(0.01)

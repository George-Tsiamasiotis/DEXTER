import pytest
from math import isfinite
from dexter.equilibrium import _TOROIDAL_TEST_NETCDF_PATH, _POLOIDAL_TEST_NETCDF_PATH
from dexter import Current, NcCurrent, LarCurrent


def test_lar_current():
    current = LarCurrent()
    assert isinstance(current.__str__(), str)
    assert isinstance(current.__repr__(), str)
    assert current.equilibrium_type == "Analytical"
    _test_current_evals(current)


def test_nc_current_getters(nc_current: NcCurrent):
    assert isinstance(nc_current.__str__(), str)
    assert isinstance(nc_current.__repr__(), str)
    assert nc_current.equilibrium_type == "Numerical"
    assert isinstance(nc_current.path, str)
    assert isinstance(nc_current.netcdf_version, str)
    assert isinstance(nc_current.equilibrium_type, str)
    assert isinstance(nc_current.interp_type, str)
    assert isfinite(nc_current.psi_wall)
    assert isfinite(nc_current.psip_wall)
    assert nc_current.psi_state == "Good"
    assert nc_current.psip_state == "Good"
    assert nc_current.psi_array.ndim == 1
    assert nc_current.psip_array.ndim == 1
    assert nc_current.g_array.ndim == 1
    assert nc_current.i_array.ndim == 1


def test_nc_current_eval(nc_current: NcCurrent):
    _test_current_evals(nc_current)


def test_toroidal_nc_current():
    nc_current = NcCurrent(_TOROIDAL_TEST_NETCDF_PATH, "Steffen")
    assert nc_current.psi_state == "Good"
    assert nc_current.psip_state == "Bad"
    assert isfinite(nc_current.g_of_psi(0.01))
    assert isfinite(nc_current.i_of_psi(0.01))
    assert isfinite(nc_current.dg_dpsi(0.01))
    assert isfinite(nc_current.di_dpsi(0.01))
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
    assert isfinite(nc_current.g_of_psip(0.01))
    assert isfinite(nc_current.i_of_psip(0.01))
    assert isfinite(nc_current.dg_dpsip(0.01))
    assert isfinite(nc_current.di_dpsip(0.01))
    with pytest.raises(BaseException):
        nc_current.g_of_psi(0.01)
    with pytest.raises(BaseException):
        nc_current.i_of_psi(0.01)
    with pytest.raises(BaseException):
        nc_current.dg_dpsi(0.01)
    with pytest.raises(BaseException):
        nc_current.di_dpsi(0.01)


def _test_current_evals(current: Current):
    psi = 0.01
    psip = 0.015
    assert isfinite(current.g_of_psi(psi))
    assert isfinite(current.g_of_psip(psip))
    assert isfinite(current.i_of_psi(psi))
    assert isfinite(current.i_of_psip(psip))
    assert isfinite(current.dg_dpsi(psi))
    assert isfinite(current.dg_dpsip(psip))
    assert isfinite(current.di_dpsi(psi))
    assert isfinite(current.di_dpsip(psip))

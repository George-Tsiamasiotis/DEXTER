import pytest
import numpy as np
from math import isfinite
from dexter.equilibrium import _TOROIDAL_TEST_NETCDF_PATH, _POLOIDAL_TEST_NETCDF_PATH
from dexter import Current, NcCurrent, LarCurrent


def test_lar_current():
    current = LarCurrent()
    assert current.equilibrium_type == "Analytical"
    assert current.psi_state == "Good"
    assert current.psip_state == "Good"
    _test_current_vectorized_evals(current)
    assert isinstance(current.__str__(), str)
    assert isinstance(current.__repr__(), str)


def test_nc_current_getters(nc_current: NcCurrent):
    assert isinstance(nc_current.path, str)
    assert isinstance(nc_current.netcdf_version, str)
    assert nc_current.equilibrium_type == "Numerical"
    assert isinstance(nc_current.interp_type, str)
    assert isfinite(nc_current.psi_last)
    assert isfinite(nc_current.psip_last)
    assert nc_current.psi_state == "Good"
    assert nc_current.psip_state == "Good"
    assert nc_current.psi_array.ndim == 1
    assert nc_current.psip_array.ndim == 1
    assert nc_current.g_array.ndim == 1
    assert nc_current.i_array.ndim == 1
    assert isinstance(nc_current.__str__(), str)
    assert isinstance(nc_current.__repr__(), str)


def test_nc_current_eval(nc_current: NcCurrent):
    _test_current_vectorized_evals(nc_current)


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


def _test_current_vectorized_evals(current: Current):
    methods = [
        current.g_of_psi,
        current.g_of_psip,
        current.i_of_psi,
        current.i_of_psip,
        current.dg_dpsi,
        current.dg_dpsip,
        current.di_dpsi,
        current.di_dpsip,
    ]

    # 0D evaluations
    flux = 1e-4
    for method in methods:
        assert isfinite(method(flux))
        assert isinstance(method(flux), float)

    # 1D Evaluations
    fluxes = np.linspace(1e-5, 1e-4, 5)
    for method in methods:
        assert method(fluxes).ndim == 1
        assert isinstance(method(fluxes), np.ndarray)

    # 4D Evaluations
    grid = np.random.random([2] * 4) * 1e-5
    assert grid.ndim == 4
    for method in methods:
        assert method(grid).ndim == 4
        assert isinstance(method(grid), np.ndarray)

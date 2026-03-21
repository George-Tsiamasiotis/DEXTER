import pytest
import numpy as np
from math import isfinite
from math import pi as PI
from dexter.equilibrium import _TOROIDAL_TEST_NETCDF_PATH, _POLOIDAL_TEST_NETCDF_PATH
from dexter import LarBfield, NcBfield


def test_lar_bfield():
    bfield = LarBfield()
    assert bfield.equilibrium_type == "Analytical"
    args = (0.01, 3.14)
    assert isfinite(bfield.b_of_psi(*args))
    assert isfinite(bfield.db_dpsi(*args))
    assert isfinite(bfield.db_of_psi_dtheta(*args))
    with pytest.raises(BaseException):
        bfield.b_of_psip(*args)
    with pytest.raises(BaseException):
        bfield.db_dpsip(*args)
    with pytest.raises(BaseException):
        bfield.db_of_psip_dtheta(*args)
    assert isinstance(bfield.__str__(), str)
    assert isinstance(bfield.__repr__(), str)


def test_nc_bfield_getters(nc_bfield: NcBfield):
    assert isinstance(nc_bfield.path, str)
    assert isinstance(nc_bfield.netcdf_version, str)
    assert nc_bfield.equilibrium_type == "Numerical"
    assert isinstance(nc_bfield.interp_type, str)
    assert isfinite(nc_bfield.baxis)
    assert len(nc_bfield.shape) == 2
    assert isfinite(nc_bfield.psi_wall)
    assert isfinite(nc_bfield.psip_wall)
    assert nc_bfield.psi_state == "Good"
    assert nc_bfield.psip_state == "Good"
    assert nc_bfield.psi_array.ndim == 1
    assert nc_bfield.psip_array.ndim == 1
    assert nc_bfield.theta_array.ndim == 1
    assert nc_bfield.b_array.ndim == 2
    assert isinstance(nc_bfield.__str__(), str)
    assert isinstance(nc_bfield.__repr__(), str)


def test_nc_bfield_eval(nc_bfield: NcBfield):
    _test_bfield_vectorized_evals(nc_bfield)


def test_toroidal_nc_bfield():
    nc_bfield = NcBfield(_TOROIDAL_TEST_NETCDF_PATH, "Bicubic")
    assert nc_bfield.psi_state == "Good"
    assert nc_bfield.psip_state == "Bad"
    args = (0.01, 3.14)
    assert isfinite(nc_bfield.b_of_psi(*args))
    assert isfinite(nc_bfield.db_dpsi(*args))
    assert isfinite(nc_bfield.db_of_psi_dtheta(*args))
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
    assert isfinite(nc_bfield.b_of_psip(*args))
    assert isfinite(nc_bfield.db_dpsip(*args))
    assert isfinite(nc_bfield.db_of_psip_dtheta(*args))
    with pytest.raises(BaseException):
        nc_bfield.b_of_psi(*args)
    with pytest.raises(BaseException):
        nc_bfield.db_dpsi(*args)
    with pytest.raises(BaseException):
        nc_bfield.db_of_psi_dtheta(*args)


def _test_bfield_vectorized_evals(bfield: NcBfield):
    methods = [
        bfield.b_of_psi,
        bfield.b_of_psip,
        bfield.db_dpsi,
        bfield.db_dpsip,
        bfield.db_of_psi_dtheta,
        bfield.db_of_psip_dtheta,
    ]

    # 0D evaluations
    flux = 1e-5
    theta = PI
    for method in methods:
        assert isfinite(method(flux, theta))
        assert isinstance(method(flux, theta), float)

    # 1D Evaluations
    fluxes = np.linspace(1e-5, 1e-4, 5)
    thetas = np.linspace(0, PI, 5)
    for method in methods:
        assert method(fluxes, thetas).ndim == 1
        assert isinstance(method(fluxes, thetas), np.ndarray)

    # 4D Evaluations
    flux_grid = np.random.random([2] * 4) * 1e-5
    theta_grid = np.random.random([2] * 4) * PI
    assert flux_grid.ndim == 4
    for method in methods:
        assert method(flux_grid, theta_grid).ndim == 4
        assert isinstance(method(flux_grid, theta_grid), np.ndarray)

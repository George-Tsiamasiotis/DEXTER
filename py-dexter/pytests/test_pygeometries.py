import pytest
import numpy as np
from math import isfinite
from math import pi as PI
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
    _test_geometry_1d_vectorized_evals(nc_geometry)
    _test_geometry_2d_vectorized_evals(nc_geometry)


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


def _test_geometry_1d_vectorized_evals(geometry: NcGeometry):
    methods = [
        geometry.psip_of_psi,
        geometry.psi_of_psip,
        geometry.r_of_psi,
        geometry.r_of_psip,
        geometry.psi_of_r,
        geometry.psip_of_r,
    ]

    # 0D evaluations
    fluxr = 1e-5
    for method in methods:
        assert isfinite(method(fluxr))
        assert isinstance(method(fluxr), float)

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


def _test_geometry_2d_vectorized_evals(geometry: NcGeometry):
    methods = [
        geometry.rlab_of_psi,
        geometry.rlab_of_psip,
        geometry.zlab_of_psi,
        geometry.zlab_of_psip,
        geometry.jacobian_of_psi,
        geometry.jacobian_of_psip,
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

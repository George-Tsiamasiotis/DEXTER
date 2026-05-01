from collections.abc import Callable
import pytest
import numpy as np
from math import isfinite
from math import pi as PI

from dexter.equilibrium import _TEST_NETCDF_PATH as netcdf_path
from dexter.equilibrium import _TOROIDAL_TEST_NETCDF_PATH as toroidal_netcdf_path
from dexter.equilibrium import _POLOIDAL_TEST_NETCDF_PATH as poloidal_netcdf_path

from dexter import Harmonic, CosHarmonic, NcHarmonic, LastClosedFluxSurface


def test_cos_harmonic_toroidal_lcfs():
    lcfs = LastClosedFluxSurface("Toroidal", 0.45)
    harmonic = CosHarmonic(1e-3, lcfs, 3, 2, 0)
    assert harmonic.equilibrium_type == "Analytical"
    assert harmonic.psi_state == "Good"
    assert harmonic.psip_state == "Bad"
    assert harmonic.epsilon == 1e-3
    assert harmonic.lcfs.kind == "Toroidal"
    assert harmonic.lcfs.value == 0.45
    assert harmonic.m == 3
    assert harmonic.n == 2
    assert harmonic.phase == 0
    _test_toroidal_harmonic_vectorized_evals(harmonic)
    assert isinstance(harmonic.__str__(), str)
    assert isinstance(harmonic.__repr__(), str)


def test_cos_harmonic_poloidal_lcfs():
    lcfs = LastClosedFluxSurface("Poloidal", 0.45)
    harmonic = CosHarmonic(1e-3, lcfs, 3, 2, 0)
    assert harmonic.psi_state == "Bad"
    assert harmonic.psip_state == "Good"
    assert harmonic.equilibrium_type == "Analytical"
    assert harmonic.epsilon == 1e-3
    assert harmonic.lcfs.kind == "Poloidal"
    assert harmonic.lcfs.value == 0.45
    assert harmonic.m == 3
    assert harmonic.n == 2
    assert harmonic.phase == 0
    _test_poloidal_harmonic_vectorized_evals(harmonic)
    assert isinstance(harmonic.__str__(), str)
    assert isinstance(harmonic.__repr__(), str)


def test_nc_harmonic_getters(nc_harmonic: NcHarmonic):
    assert nc_harmonic.equilibrium_type == "Numerical"
    assert isinstance(nc_harmonic.path, str)
    assert isinstance(nc_harmonic.netcdf_version, str)
    assert isinstance(nc_harmonic.equilibrium_type, str)
    assert isinstance(nc_harmonic.interp_type, str)
    assert isinstance(nc_harmonic.m, int)
    assert isinstance(nc_harmonic.n, int)
    assert nc_harmonic.phase_method == "Zero"
    assert nc_harmonic.analytical_threshold_index == 3
    assert isfinite(nc_harmonic.psi_last)
    assert isfinite(nc_harmonic.psip_last)
    assert nc_harmonic.psi_state == "Good"
    assert nc_harmonic.psip_state == "Good"
    assert nc_harmonic.psi_array.ndim == 1
    assert nc_harmonic.psip_array.ndim == 1
    assert nc_harmonic.alpha_array.ndim == 1
    assert nc_harmonic.phase_array.ndim == 1
    assert isinstance(nc_harmonic.__str__(), str)
    assert isinstance(nc_harmonic.__repr__(), str)


def test_nc_harmonic_eval(nc_harmonic: NcHarmonic):
    _test_toroidal_harmonic_vectorized_evals(nc_harmonic)
    _test_poloidal_harmonic_vectorized_evals(nc_harmonic)


def test_nc_harmonic_phase_methods():
    harmonic = NcHarmonic(netcdf_path, "Cubic", 3, 2, "Zero")
    assert harmonic.phase_method == "Zero"
    assert harmonic.phase_of_psi(1e-5, 0, 0, 0) == 0
    _test_toroidal_harmonic_vectorized_evals(harmonic)
    _test_poloidal_harmonic_vectorized_evals(harmonic)

    harmonic = NcHarmonic(netcdf_path, "Cubic", 3, 2, "Average")
    assert harmonic.phase_method == "Average"
    _test_toroidal_harmonic_vectorized_evals(harmonic)
    _test_poloidal_harmonic_vectorized_evals(harmonic)

    harmonic = NcHarmonic(netcdf_path, "Cubic", 3, 2, "Interpolation")
    assert harmonic.phase_method == "Interpolation"
    _test_toroidal_harmonic_vectorized_evals(harmonic)
    _test_poloidal_harmonic_vectorized_evals(harmonic)

    # harmonic = NcHarmonic(poloidal_netcdf_path, "Cubic", 3, 2, "Resonance")
    # assert harmonic.phase_method == "Resonance"
    # # Not all evals work in this one of course

    harmonic = NcHarmonic(netcdf_path, "Cubic", 3, 2, ("Custom", 10))
    assert harmonic.phase_method == "Custom(10.0)"
    assert harmonic.phase_of_psi(1e-5, 0, 0, 0) == 10
    _test_toroidal_harmonic_vectorized_evals(harmonic)
    _test_poloidal_harmonic_vectorized_evals(harmonic)

    with pytest.raises(BaseException):
        NcHarmonic(netcdf_path, "Cubic", 3, 2, "not a method")  # type: ignore
    with pytest.raises(BaseException):
        NcHarmonic(netcdf_path, "Cubic", 3, 2, ("Not Custom", 10))  # type: ignore
    with pytest.raises(BaseException):
        NcHarmonic(netcdf_path, "Cubic", 3, 2, ("Not Custom"))  # type: ignore
    with pytest.raises(BaseException):
        NcHarmonic(netcdf_path, "Cubic", 3, 2, ("Zero", 10))  # type: ignore


def test_toroidal_nc_harmonic():
    harmonic = NcHarmonic(toroidal_netcdf_path, "Cubic", 3, 2, "Interpolation")
    assert harmonic.psi_state == "Good"
    assert harmonic.psip_state == "Bad"
    args = (1e-5, 0.0, 0.0, 0.0)
    assert isfinite(harmonic.alpha_of_psi(*args))
    assert isfinite(harmonic.phase_of_psi(*args))
    assert isfinite(harmonic.h_of_psi(*args))
    assert isfinite(harmonic.dh_dpsi(*args))
    assert isfinite(harmonic.dh_of_psi_dtheta(*args))
    assert isfinite(harmonic.dh_of_psi_dzeta(*args))
    with pytest.raises(BaseException):
        harmonic.alpha_of_psip(*args)
    with pytest.raises(BaseException):
        harmonic.phase_of_psip(*args)
    with pytest.raises(BaseException):
        harmonic.h_of_psip(*args)
    with pytest.raises(BaseException):
        harmonic.dh_dpsip(*args)
    with pytest.raises(BaseException):
        harmonic.dh_of_psip_dtheta(*args)
    with pytest.raises(BaseException):
        harmonic.dh_of_psip_dzeta(*args)

    assert isfinite(harmonic.dh_of_psi_dt(*args))
    assert isfinite(harmonic.dh_of_psip_dt(*args))


def test_poloidal_nc_harmonic():
    harmonic = NcHarmonic(poloidal_netcdf_path, "Cubic", 3, 2, "Interpolation")
    assert harmonic.psi_state == "Bad"
    assert harmonic.psip_state == "Good"
    args = (1e-5, 0.0, 0.0, 0.0)
    assert isfinite(harmonic.alpha_of_psip(*args))
    assert isfinite(harmonic.phase_of_psip(*args))
    assert isfinite(harmonic.h_of_psip(*args))
    assert isfinite(harmonic.dh_dpsip(*args))
    assert isfinite(harmonic.dh_of_psip_dtheta(*args))
    assert isfinite(harmonic.dh_of_psip_dzeta(*args))
    with pytest.raises(BaseException):
        harmonic.alpha_of_psi(*args)
    with pytest.raises(BaseException):
        harmonic.phase_of_psi(*args)
    with pytest.raises(BaseException):
        harmonic.h_of_psi(*args)
    with pytest.raises(BaseException):
        harmonic.dh_dpsi(*args)
    with pytest.raises(BaseException):
        harmonic.dh_of_psi_dtheta(*args)
    with pytest.raises(BaseException):
        harmonic.dh_of_psi_dzeta(*args)

    assert isfinite(harmonic.dh_of_psi_dt(*args))
    assert isfinite(harmonic.dh_of_psip_dt(*args))


def _test_toroidal_harmonic_vectorized_evals(harmonic: Harmonic):
    methods = [
        harmonic.alpha_of_psi,
        harmonic.phase_of_psi,
        harmonic.h_of_psi,
        harmonic.dh_dpsi,
        harmonic.dh_of_psi_dtheta,
        harmonic.dh_of_psi_dzeta,
        harmonic.dh_of_psi_dt,
    ]
    _test_harmonic_vectorized_evals(methods)


def _test_poloidal_harmonic_vectorized_evals(harmonic: Harmonic):
    methods = [
        harmonic.alpha_of_psip,
        harmonic.phase_of_psip,
        harmonic.h_of_psip,
        harmonic.dh_dpsip,
        harmonic.dh_of_psip_dtheta,
        harmonic.dh_of_psip_dzeta,
        harmonic.dh_of_psip_dt,
    ]
    _test_harmonic_vectorized_evals(methods)


def _test_harmonic_vectorized_evals(methods: list[Callable]):
    # 0D evaluations
    flux = 1e-5
    theta = PI
    zeta = 2 * PI
    t = 10
    for method in methods:
        assert isfinite(method(flux, theta, zeta, t))
        assert isinstance(method(flux, theta, zeta, t), float)

    # 1D Evaluations
    fluxes = np.linspace(1e-5, 1e-4, 5)
    thetas = np.linspace(0, PI, 5)
    zetas = np.linspace(0, 2 * PI, 5)
    ts = np.linspace(0, 100, 5)
    for method in methods:
        assert method(fluxes, thetas, zetas, ts).ndim == 1
        assert isinstance(method(fluxes, thetas, zetas, ts), np.ndarray)

    # 4D Evaluations
    flux_grid = np.random.random([2] * 4) * 1e-5
    theta_grid = np.random.random([2] * 4) * PI
    zeta_grid = np.random.random([2] * 4) * 2 * PI
    t_grid = np.random.random([2] * 4) * 1e-5
    assert flux_grid.ndim == 4
    for method in methods:
        assert method(flux_grid, theta_grid, zeta_grid, t_grid).ndim == 4
        assert isinstance(method(flux_grid, theta_grid, zeta_grid, t_grid), np.ndarray)

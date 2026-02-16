import pytest

from dexter.equilibrium import _TEST_NETCDF_PATH as netcdf_path
from dexter.equilibrium import _TOROIDAL_TEST_NETCDF_PATH as toroidal_netcdf_path
from dexter.equilibrium import _POLOIDAL_TEST_NETCDF_PATH as poloidal_netcdf_path

from dexter import Harmonic, CosHarmonic, NcHarmonic


def test_cos_harmonic():
    harmonic = CosHarmonic(1e-3, 3, 2, 0)
    assert isinstance(harmonic.__str__(), str)
    assert isinstance(harmonic.__repr__(), str)
    assert harmonic.alpha == 1e-3
    assert harmonic.m == 3
    assert harmonic.n == 2
    assert harmonic.phase == 0
    _test_harmonic_evals(harmonic)


def test_nc_harmonic_getters(nc_harmonic: NcHarmonic):
    assert isinstance(nc_harmonic.__str__(), str)
    assert isinstance(nc_harmonic.__repr__(), str)
    assert nc_harmonic.equilibrium_type == "Numerical"
    assert isinstance(nc_harmonic.path, str)
    assert isinstance(nc_harmonic.netcdf_version, str)
    assert isinstance(nc_harmonic.equilibrium_type, str)
    assert isinstance(nc_harmonic.interp_type, str)
    assert isinstance(nc_harmonic.m, int)
    assert isinstance(nc_harmonic.n, int)
    assert nc_harmonic.phase_method == "Zero"
    assert isinstance(nc_harmonic.psi_wall, float)
    assert isinstance(nc_harmonic.psip_wall, float)
    assert nc_harmonic.psi_state == "Good"
    assert nc_harmonic.psip_state == "Good"
    assert nc_harmonic.psi_array.ndim == 1
    assert nc_harmonic.psip_array.ndim == 1
    assert nc_harmonic.alpha_array.ndim == 1
    assert nc_harmonic.phase_array.ndim == 1


def test_nc_harmonic_eval(nc_harmonic: NcHarmonic):
    _test_harmonic_evals(nc_harmonic)


def test_nc_harmonic_phase_methods():
    harmonic = NcHarmonic(netcdf_path, "Cubic", 3, 2, "Zero")
    assert harmonic.phase_method == "Zero"
    assert harmonic.phase_of_psi(0, 0, 0, 0) == 0
    _test_harmonic_evals(harmonic)

    harmonic = NcHarmonic(netcdf_path, "Cubic", 3, 2, "Average")
    assert harmonic.phase_method == "Average"
    _test_harmonic_evals(harmonic)

    harmonic = NcHarmonic(netcdf_path, "Cubic", 3, 2, "Interpolation")
    assert harmonic.phase_method == "Interpolation"
    _test_harmonic_evals(harmonic)

    # harmonic = NcHarmonic(poloidal_netcdf_path, "Cubic", 3, 2, "Resonance")
    # assert harmonic.phase_method == "Resonance"
    # # Not all evals work in this one of course

    harmonic = NcHarmonic(netcdf_path, "Cubic", 3, 2, ("Custom", 10))
    assert harmonic.phase_method == "Custom(10.0)"
    assert harmonic.phase_of_psi(0, 0, 0, 0) == 10
    _test_harmonic_evals(harmonic)

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
    args = (0.0, 0.0, 0.0, 0.0)
    assert isinstance(harmonic.alpha_of_psi(*args), float)
    assert isinstance(harmonic.phase_of_psi(*args), float)
    assert isinstance(harmonic.h_of_psi(*args), float)
    assert isinstance(harmonic.dh_dpsi(*args), float)
    assert isinstance(harmonic.dh_of_psi_dtheta(*args), float)
    assert isinstance(harmonic.dh_of_psi_dzeta(*args), float)
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

    assert isinstance(harmonic.dh_of_psi_dt(*args), float)
    assert isinstance(harmonic.dh_of_psip_dt(*args), float)


def test_poloidal_nc_harmonic():
    harmonic = NcHarmonic(poloidal_netcdf_path, "Cubic", 3, 2, "Interpolation")
    assert harmonic.psi_state == "Bad"
    assert harmonic.psip_state == "Good"
    args = (0.0, 0.0, 0.0, 0.0)
    assert isinstance(harmonic.alpha_of_psip(*args), float)
    assert isinstance(harmonic.phase_of_psip(*args), float)
    assert isinstance(harmonic.h_of_psip(*args), float)
    assert isinstance(harmonic.dh_dpsip(*args), float)
    assert isinstance(harmonic.dh_of_psip_dtheta(*args), float)
    assert isinstance(harmonic.dh_of_psip_dzeta(*args), float)
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

    assert isinstance(harmonic.dh_of_psi_dt(*args), float)
    assert isinstance(harmonic.dh_of_psip_dt(*args), float)


def _test_harmonic_evals(harmonic: Harmonic):
    psi = 0.01
    psip = 0.015
    theta = 3.1415
    zeta = 6.2831
    t = 0
    assert isinstance(harmonic.alpha_of_psi(psi, theta, zeta, t), float)
    assert isinstance(harmonic.alpha_of_psip(psip, theta, zeta, t), float)
    assert isinstance(harmonic.phase_of_psi(psi, theta, zeta, t), float)
    assert isinstance(harmonic.phase_of_psip(psip, theta, zeta, t), float)
    assert isinstance(harmonic.h_of_psi(psi, theta, zeta, t), float)
    assert isinstance(harmonic.h_of_psip(psip, theta, zeta, t), float)
    assert isinstance(harmonic.dh_dpsi(psi, theta, zeta, t), float)
    assert isinstance(harmonic.dh_dpsip(psip, theta, zeta, t), float)
    assert isinstance(harmonic.dh_of_psi_dtheta(psi, theta, zeta, t), float)
    assert isinstance(harmonic.dh_of_psip_dtheta(psip, theta, zeta, t), float)
    assert isinstance(harmonic.dh_of_psi_dzeta(psi, theta, zeta, t), float)
    assert isinstance(harmonic.dh_of_psip_dzeta(psip, theta, zeta, t), float)
    assert isinstance(harmonic.dh_of_psi_dt(psi, theta, zeta, t), float)
    assert isinstance(harmonic.dh_of_psip_dt(psip, theta, zeta, t), float)

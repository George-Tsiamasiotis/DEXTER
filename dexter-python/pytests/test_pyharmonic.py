from dexter import NcHarmonic
from dexter import _LAR_NETCDF_PATH as netcdf_path


def test_nc_harmonic_getters(nc_harmonic1: NcHarmonic):
    assert isinstance(nc_harmonic1.path, str)
    assert isinstance(nc_harmonic1.typ, str)
    assert isinstance(nc_harmonic1.phase_average, float | None)
    assert isinstance(nc_harmonic1.phase_resonance, float | None)
    assert isinstance(nc_harmonic1.phase_method, str | None)
    assert isinstance(nc_harmonic1.m, int)
    assert isinstance(nc_harmonic1.n, int)
    assert nc_harmonic1.psip_data.ndim == 1
    assert nc_harmonic1.a_data.ndim == 1
    assert nc_harmonic1.phase_data.ndim == 1
    assert len(nc_harmonic1) == len(nc_harmonic1.psip_data)
    str(nc_harmonic1)


def test_nc_harmonic_evals(nc_harmonic1: NcHarmonic):
    (psip, theta, zeta) = 0.015, 3.14, 0
    assert isinstance(nc_harmonic1.h(psip, theta, zeta), float)
    assert isinstance(nc_harmonic1.dh_dpsip(psip, theta, zeta), float)
    assert isinstance(nc_harmonic1.dh_dtheta(psip, theta, zeta), float)
    assert isinstance(nc_harmonic1.dh_dzeta(psip, theta, zeta), float)
    assert isinstance(nc_harmonic1.dh_dt(psip, theta, zeta), float)
    assert isinstance(nc_harmonic1.a(psip), float)
    assert isinstance(nc_harmonic1.da_dpsip(psip), float)
    assert isinstance(nc_harmonic1.phase(psip), float)


def test_nc_harmonic_phase_methods():
    NcHarmonic(netcdf_path, "Steffen", 2, 1, phase_method="Zero")
    NcHarmonic(netcdf_path, "Steffen", 2, 1, phase_method="Average")
    NcHarmonic(netcdf_path, "Steffen", 2, 1, phase_method="Resonance")
    NcHarmonic(netcdf_path, "Steffen", 2, 1, phase_method="Interpolation")

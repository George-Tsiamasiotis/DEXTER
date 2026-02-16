import pytest
from math import isfinite
from dexter.equilibrium import _TOROIDAL_TEST_NETCDF_PATH, _POLOIDAL_TEST_NETCDF_PATH
from dexter import Qfactor, UnityQfactor, ParabolicQfactor, NcQfactor
from dexter.types import FluxWall


def test_unity_qfactor():
    qfactor = UnityQfactor()
    assert isinstance(qfactor.__str__(), str)
    assert isinstance(qfactor.__repr__(), str)
    assert qfactor.equilibrium_type == "Analytical"
    _test_qfactor_evals(qfactor)


def test_parabolic_qfactor():
    flux_wall: FluxWall = ("Toroidal", 0.45)
    qfactor = ParabolicQfactor(1.1, 3.8, flux_wall)
    assert isinstance(qfactor.__str__(), str)
    assert isinstance(qfactor.__repr__(), str)
    assert qfactor.equilibrium_type == "Analytical"
    assert qfactor.qaxis == 1.1
    assert qfactor.qwall == 3.8
    assert qfactor.psi_wall == 0.45
    assert isfinite(qfactor.psip_wall)
    _test_qfactor_evals(qfactor)


def test_nc_qfactor_getters(nc_qfactor: NcQfactor):
    assert isinstance(nc_qfactor.__str__(), str)
    assert isinstance(nc_qfactor.__repr__(), str)
    assert nc_qfactor.equilibrium_type == "Numerical"
    assert isinstance(nc_qfactor.path, str)
    assert isinstance(nc_qfactor.netcdf_version, str)
    assert isinstance(nc_qfactor.equilibrium_type, str)
    assert isinstance(nc_qfactor.interp_type, str)
    assert isfinite(nc_qfactor.qaxis)
    assert isfinite(nc_qfactor.qwall)
    assert isfinite(nc_qfactor.psi_wall)
    assert isfinite(nc_qfactor.psip_wall)
    assert nc_qfactor.psi_state == "Good"
    assert nc_qfactor.psip_state == "Good"
    assert nc_qfactor.psi_array.ndim == 1
    assert nc_qfactor.psip_array.ndim == 1
    assert nc_qfactor.q_array.ndim == 1


def test_nc_qfactor_eval(nc_qfactor: NcQfactor):
    _test_qfactor_evals(nc_qfactor)


def test_toroidal_nc_qfactor():
    nc_qfactor = NcQfactor(_TOROIDAL_TEST_NETCDF_PATH, "Steffen")
    assert nc_qfactor.psi_state == "Good"
    assert nc_qfactor.psip_state == "Bad"
    assert isfinite(nc_qfactor.q_of_psi(0.01))
    assert isfinite(nc_qfactor.iota_of_psi(0.01))
    assert isfinite(nc_qfactor.psip_of_psi(0.01))
    assert isfinite(nc_qfactor.dpsip_dpsi(0.01))
    with pytest.raises(BaseException):
        nc_qfactor.q_of_psip(0.01)
    with pytest.raises(BaseException):
        nc_qfactor.iota_of_psip(0.01)
    with pytest.raises(BaseException):
        nc_qfactor.psi_of_psip(0.01)
    with pytest.raises(BaseException):
        nc_qfactor.dpsi_dpsip(0.01)


def test_poloidal_nc_qfactor():
    nc_qfactor = NcQfactor(_POLOIDAL_TEST_NETCDF_PATH, "Steffen")
    assert nc_qfactor.psi_state == "Bad"
    assert nc_qfactor.psip_state == "Good"
    assert isfinite(nc_qfactor.q_of_psip(0.01))
    assert isfinite(nc_qfactor.iota_of_psip(0.01))
    assert isfinite(nc_qfactor.psi_of_psip(0.01))
    assert isfinite(nc_qfactor.dpsi_dpsip(0.01))
    with pytest.raises(BaseException):
        nc_qfactor.q_of_psi(0.01)
    with pytest.raises(BaseException):
        nc_qfactor.iota_of_psi(0.01)
    with pytest.raises(BaseException):
        nc_qfactor.psip_of_psi(0.01)
    with pytest.raises(BaseException):
        nc_qfactor.dpsip_dpsi(0.01)


def _test_qfactor_evals(qfactor: Qfactor):
    psi = 0.01
    psip = 0.015
    assert isfinite(qfactor.q_of_psi(psi))
    assert isfinite(qfactor.q_of_psip(psip))
    assert isfinite(qfactor.psip_of_psi(psi))
    assert isfinite(qfactor.psi_of_psip(psip))
    assert isfinite(qfactor.iota_of_psi(psi))
    assert isfinite(qfactor.iota_of_psip(psip))
    assert isfinite(qfactor.dpsip_dpsi(psi))
    assert isfinite(qfactor.dpsi_dpsip(psip))

import pytest
from dexter.equilibrium import _TOROIDAL_TEST_NETCDF_PATH, _POLOIDAL_TEST_NETCDF_PATH
from dexter import Qfactor, UnityQfactor, ParabolicQfactor, NcQfactor


def test_unity_qfactor():
    qfactor = UnityQfactor()
    assert isinstance(qfactor.__str__(), str)
    assert isinstance(qfactor.__repr__(), str)
    assert qfactor.equilibrium_type == "Analytical"
    _test_qfactor_evals(qfactor)


def test_parabolic_qfactor():
    qfactor = ParabolicQfactor(1.1, 3.8, 0.45)
    assert isinstance(qfactor.__str__(), str)
    assert isinstance(qfactor.__repr__(), str)
    assert qfactor.equilibrium_type == "Analytical"
    assert qfactor.qaxis == 1.1
    assert qfactor.qwall == 3.8
    assert qfactor.psi_wall == 0.45
    assert isinstance(qfactor.psip_wall, float)
    _test_qfactor_evals(qfactor)


def test_nc_qfactor_getters(nc_qfactor: NcQfactor):
    assert isinstance(nc_qfactor.__str__(), str)
    assert isinstance(nc_qfactor.__repr__(), str)
    assert nc_qfactor.equilibrium_type == "Numerical"
    assert isinstance(nc_qfactor.path, str)
    assert isinstance(nc_qfactor.netcdf_version, str)
    assert isinstance(nc_qfactor.equilibrium_type, str)
    assert isinstance(nc_qfactor.interp_type, str)
    assert isinstance(nc_qfactor.qaxis, float)
    assert isinstance(nc_qfactor.qwall, float)
    assert isinstance(nc_qfactor.psi_wall, float)
    assert isinstance(nc_qfactor.psip_wall, float)
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
    assert isinstance(nc_qfactor.q_of_psi(0.01), float)
    assert isinstance(nc_qfactor.iota_of_psi(0.01), float)
    assert isinstance(nc_qfactor.psip_of_psi(0.01), float)
    assert isinstance(nc_qfactor.dpsip_dpsi(0.01), float)
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
    assert isinstance(nc_qfactor.q_of_psip(0.01), float)
    assert isinstance(nc_qfactor.iota_of_psip(0.01), float)
    assert isinstance(nc_qfactor.psi_of_psip(0.01), float)
    assert isinstance(nc_qfactor.dpsi_dpsip(0.01), float)
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
    assert isinstance(qfactor.q_of_psi(psi), float)
    assert isinstance(qfactor.q_of_psip(psip), float)
    assert isinstance(qfactor.psip_of_psi(psi), float)
    assert isinstance(qfactor.psi_of_psip(psip), float)
    assert isinstance(qfactor.iota_of_psi(psi), float)
    assert isinstance(qfactor.iota_of_psip(psip), float)
    assert isinstance(qfactor.dpsip_dpsi(psi), float)
    assert isinstance(qfactor.dpsi_dpsip(psip), float)

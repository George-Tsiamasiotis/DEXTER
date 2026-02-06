import pytest
from dexter import (
    _TOROIDAL_TEST_NETCDF_PATH,
    _POLOIDAL_TEST_NETCDF_PATH,
    UnityQfactor,
    ParabolicQfactor,
    NcQfactor,
)


def test_nc_qfactor_getters(nc_qfactor: NcQfactor):
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
    assert isinstance(nc_qfactor.__str__(), str)
    assert isinstance(nc_qfactor.__repr__(), str)


def test_nc_qfactor_eval(nc_qfactor: NcQfactor):
    psi = 0.01
    psip = 0.015
    assert isinstance(nc_qfactor.q_of_psi(psi), float)
    assert isinstance(nc_qfactor.q_of_psip(psip), float)
    assert isinstance(nc_qfactor.psip_of_psi(psi), float)
    assert isinstance(nc_qfactor.psi_of_psip(psip), float)
    assert isinstance(nc_qfactor.iota_of_psi(psi), float)
    assert isinstance(nc_qfactor.iota_of_psip(psip), float)
    assert isinstance(nc_qfactor.dpsip_dpsi(psi), float)
    assert isinstance(nc_qfactor.dpsi_dpsip(psip), float)


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


def test_unity_qfactor():
    qfactor = UnityQfactor()
    assert qfactor.equilibrium_type == "Analytical"
    p = 0.01
    assert qfactor.q_of_psi(p) == 1
    assert qfactor.q_of_psip(p) == 1
    assert qfactor.psip_of_psi(p) == p
    assert qfactor.psi_of_psip(p) == p
    assert qfactor.dpsip_dpsi(p) == 1
    assert qfactor.dpsi_dpsip(p) == 1
    assert qfactor.iota_of_psi(p) == 1
    assert qfactor.iota_of_psip(p) == 1


def test_parabolic_qfactor():
    qfactor = ParabolicQfactor(1.1, 3.8, 0.45)
    assert qfactor.equilibrium_type == "Analytical"
    assert qfactor.qaxis == 1.1
    assert qfactor.qwall == 3.8
    assert qfactor.psi_wall == 0.45
    assert isinstance(qfactor.psip_wall, float)
    assert isinstance(qfactor.q_of_psi(0.01), float)
    assert isinstance(qfactor.q_of_psip(0.01), float)
    assert isinstance(qfactor.psip_of_psi(0.01), float)
    assert isinstance(qfactor.psi_of_psip(0.01), float)
    assert isinstance(qfactor.dpsip_dpsi(0.01), float)
    assert isinstance(qfactor.dpsi_dpsip(0.01), float)
    assert isinstance(qfactor.iota_of_psi(0.01), float)
    assert isinstance(qfactor.iota_of_psip(0.01), float)

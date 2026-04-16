import pytest
import numpy as np
from math import isfinite
from dexter.equilibrium import _TOROIDAL_TEST_NETCDF_PATH, _POLOIDAL_TEST_NETCDF_PATH
from dexter import (
    Qfactor,
    UnityQfactor,
    ParabolicQfactor,
    NcQfactor,
    LastClosedFluxSurface,
)


def test_unity_qfactor():
    LCFS = LastClosedFluxSurface("Toroidal", 0.45)
    qfactor = UnityQfactor(LCFS)
    assert qfactor.equilibrium_type == "Analytical"
    assert qfactor.psi_state == "Good"
    assert qfactor.psip_state == "Good"
    assert qfactor.psi_last == 0.45
    assert qfactor.psip_last == 0.45
    assert qfactor.qlast == 1.0
    assert qfactor.qaxis == 1.0
    _test_qfactor_vectorized_evals(qfactor)
    assert isinstance(qfactor.__str__(), str)
    assert isinstance(qfactor.__repr__(), str)


def test_parabolic_qfactor_toroidal_lcfs():
    LCFS = LastClosedFluxSurface("Toroidal", 0.45)
    qfactor = ParabolicQfactor(1.1, 3.8, LCFS)
    assert qfactor.equilibrium_type == "Analytical"
    assert qfactor.psi_state == "Good"
    assert qfactor.psip_state == "Good"
    assert qfactor.qaxis == 1.1
    assert qfactor.qlast == 3.8
    assert qfactor.psi_last == 0.45
    assert isfinite(qfactor.psip_last)
    _test_qfactor_vectorized_evals(qfactor)
    _test_inverse_qfactor_vectorized_evals(qfactor)
    assert isinstance(qfactor.__str__(), str)
    assert isinstance(qfactor.__repr__(), str)


def test_parabolic_qfactor_poloidal_lcfs():
    LCFS = LastClosedFluxSurface("Poloidal", 0.45)
    qfactor = ParabolicQfactor(1.1, 3.8, LCFS)
    assert qfactor.equilibrium_type == "Analytical"
    assert qfactor.psi_state == "Good"
    assert qfactor.psip_state == "Good"
    assert qfactor.qaxis == 1.1
    assert qfactor.qlast == 3.8
    assert isfinite(qfactor.psip_last)
    assert qfactor.psip_last == 0.45
    _test_qfactor_vectorized_evals(qfactor)
    _test_inverse_qfactor_vectorized_evals(qfactor)
    assert isinstance(qfactor.__str__(), str)
    assert isinstance(qfactor.__repr__(), str)


def test_nc_qfactor_getters(nc_qfactor: NcQfactor):
    assert nc_qfactor.equilibrium_type == "Numerical"
    assert isinstance(nc_qfactor.path, str)
    assert isinstance(nc_qfactor.netcdf_version, str)
    assert isinstance(nc_qfactor.equilibrium_type, str)
    assert isinstance(nc_qfactor.interp_type, str)
    assert isfinite(nc_qfactor.qaxis)
    assert isfinite(nc_qfactor.qlast)
    assert isfinite(nc_qfactor.psi_last)
    assert isfinite(nc_qfactor.psip_last)
    assert nc_qfactor.psi_state == "Good"
    assert nc_qfactor.psip_state == "Good"
    assert nc_qfactor.psi_array.ndim == 1
    assert nc_qfactor.psip_array.ndim == 1
    assert nc_qfactor.q_array.ndim == 1
    assert isinstance(nc_qfactor.__str__(), str)
    assert isinstance(nc_qfactor.__repr__(), str)


def test_nc_qfactor_eval(nc_qfactor: NcQfactor):
    _test_qfactor_vectorized_evals(nc_qfactor)
    _test_inverse_qfactor_vectorized_evals(nc_qfactor)


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


def _test_qfactor_vectorized_evals(qfactor: Qfactor):
    methods = [
        qfactor.q_of_psi,
        qfactor.q_of_psip,
        qfactor.psip_of_psi,
        qfactor.psip_of_psi,
        qfactor.iota_of_psi,
        qfactor.iota_of_psip,
        qfactor.dpsip_dpsi,
        qfactor.dpsi_dpsip,
    ]

    # 0D evaluations
    flux = 1e-5
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


def _test_inverse_qfactor_vectorized_evals(qfactor: Qfactor):
    methods = [
        qfactor.psi_of_q,
        qfactor.psip_of_q,
    ]

    # 0D evaluations
    q = 1.1
    for method in methods:
        assert isfinite(method(q))
        assert isinstance(method(q), float)

    # 1D Evaluations
    fluxes = np.linspace(1.1, 1.4, 5)
    for method in methods:
        assert method(fluxes).ndim == 1
        assert isinstance(method(fluxes), np.ndarray)

    # 4D Evaluations
    grid = 1.1 + np.random.random([2] * 4)
    assert grid.ndim == 4
    for method in methods:
        assert method(grid).ndim == 4
        assert isinstance(method(grid), np.ndarray)

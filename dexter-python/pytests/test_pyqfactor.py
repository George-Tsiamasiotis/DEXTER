import pytest
from dexter import Qfactor


def test_pyqfactor_derived_fields(qfactor: Qfactor):
    assert isinstance(qfactor.path, str)
    assert isinstance(qfactor.typ, str)
    assert isinstance(qfactor.psip_wall, float)
    assert isinstance(qfactor.psi_wall, float)


def test_pyqfactor_eval(qfactor: Qfactor):
    psip = 0.015
    assert isinstance(qfactor.q(psip), float)
    assert isinstance(qfactor.r(psip), float)
    assert isinstance(qfactor.psi(psip), float)


def test_data_extraction(qfactor: Qfactor):
    assert qfactor.psip_data.ndim == 1
    assert qfactor.q_data.ndim == 1
    assert qfactor.r_data.ndim == 1
    assert qfactor.psi_data.ndim == 1
    assert qfactor.q_data_derived.ndim == 1


def test_immutability(qfactor: Qfactor):
    with pytest.raises(AttributeError):
        qfactor.psip_data *= 2
    with pytest.raises(AttributeError):
        qfactor.q_data *= 2
    with pytest.raises(AttributeError):
        qfactor.r_data *= 2
    with pytest.raises(AttributeError):
        qfactor.psi_data *= 2
    with pytest.raises(AttributeError):
        qfactor.q_data_derived *= 2
    with pytest.raises(AttributeError):
        qfactor.psip_wall += 1
    with pytest.raises(AttributeError):
        qfactor.psi_wall += 1
    with pytest.raises(AttributeError):
        qfactor.path = ""
    with pytest.raises(AttributeError):
        qfactor.typ = ""


def test_repr(qfactor: Qfactor):
    str(qfactor)

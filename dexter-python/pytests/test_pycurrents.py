import pytest
from dexter import Currents


def test_derived_fields(currents: Currents):
    assert isinstance(currents.path, str)
    assert isinstance(currents.typ, str)
    assert isinstance(currents.psip_wall, float)


def test_eval(currents: Currents):
    psip = 0.015
    assert isinstance(currents.g(psip), float)
    assert isinstance(currents.i(psip), float)
    assert isinstance(currents.dg_dpsip(psip), float)
    assert isinstance(currents.di_dpsip(psip), float)


def test_data_extraction(currents: Currents):
    assert currents.psip_data.ndim == 1
    assert currents.g_data.ndim == 1
    assert currents.i_data.ndim == 1


def test_immutability(currents: Currents):
    """Tests that current fields are immutable."""
    with pytest.raises(AttributeError):
        currents.psip_data *= 2
    with pytest.raises(AttributeError):
        currents.g_data *= 2
    with pytest.raises(AttributeError):
        currents.i_data *= 2
    with pytest.raises(AttributeError):
        currents.psip_wall += 1
    with pytest.raises(AttributeError):
        currents.path = ""
    with pytest.raises(AttributeError):
        currents.typ = ""


def test_repr(currents: Currents):
    str(currents)

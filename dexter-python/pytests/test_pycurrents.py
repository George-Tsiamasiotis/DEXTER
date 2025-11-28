from dexter import Currents


def test_fields(currents: Currents):
    assert isinstance(currents.path, str)
    assert isinstance(currents.typ, str)
    assert isinstance(currents.psip_wall, float)
    assert currents.psip_data.ndim == 1
    assert currents.g_data.ndim == 1
    assert currents.i_data.ndim == 1


def test_eval(currents: Currents):
    psip = 0.015
    assert isinstance(currents.g(psip), float)
    assert isinstance(currents.i(psip), float)
    assert isinstance(currents.dg_dpsip(psip), float)
    assert isinstance(currents.di_dpsip(psip), float)


def test_magic(currents: Currents):
    assert len(currents) == len(currents.psip_data)
    assert isinstance(currents.__str__(), str)
    assert isinstance(currents.__repr__(), str)

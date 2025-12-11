from dexter import LarCurrent, NcCurrent


def test_nc_current_getters(nc_current: NcCurrent):
    assert isinstance(nc_current.path, str)
    assert isinstance(nc_current.typ, str)
    assert nc_current.psip_data.ndim == 1
    assert nc_current.g_data.ndim == 1
    assert nc_current.i_data.ndim == 1
    assert len(nc_current) == len(nc_current.psip_data)
    assert isinstance(nc_current.__str__(), str)
    assert isinstance(nc_current.__repr__(), str)


def test_nc_current_eval(nc_current: NcCurrent):
    psip = 0.015
    assert isinstance(nc_current.g(psip), float)
    assert isinstance(nc_current.i(psip), float)
    assert isinstance(nc_current.dg_dpsip(psip), float)
    assert isinstance(nc_current.di_dpsip(psip), float)


def test_lar_current(lar_current: LarCurrent):
    psip = 0.015
    assert lar_current.g(psip) == 1
    assert lar_current.i(psip) == 0
    assert lar_current.dg_dpsip(psip) == 0
    assert lar_current.di_dpsip(psip) == 0
    str(lar_current)

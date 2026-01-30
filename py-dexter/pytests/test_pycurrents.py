from dexter import NcCurrent


def test_nc_current_getters(nc_current: NcCurrent):
    assert isinstance(nc_current.path, str)
    assert isinstance(nc_current.typ, str)
    assert isinstance(nc_current.psi_wall, float)
    assert isinstance(nc_current.psip_wall, float)
    assert nc_current.psi_array.ndim == 1
    assert nc_current.psip_array.ndim == 1
    assert nc_current.g_array.ndim == 1
    assert nc_current.i_array.ndim == 1
    assert len(nc_current) == len(nc_current.psi_array)
    assert isinstance(nc_current.__str__(), str)
    assert isinstance(nc_current.__repr__(), str)


def test_nc_current_eval(nc_current: NcCurrent):
    psi = 0.01
    psip = 0.015
    assert isinstance(nc_current.g_of_psi(psi), float)
    assert isinstance(nc_current.g_of_psip(psip), float)
    assert isinstance(nc_current.i_of_psi(psi), float)
    assert isinstance(nc_current.i_of_psip(psip), float)
    assert isinstance(nc_current.dg_dpsi(psi), float)
    assert isinstance(nc_current.dg_dpsip(psip), float)
    assert isinstance(nc_current.di_dpsi(psi), float)
    assert isinstance(nc_current.di_dpsip(psip), float)

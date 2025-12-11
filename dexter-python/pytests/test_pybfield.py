from dexter import NcBfield


def test_nc_bfield_getters(nc_bfield: NcBfield):
    assert isinstance(nc_bfield.path, str)
    assert isinstance(nc_bfield.typ, str)
    assert nc_bfield.psip_data.ndim == 1
    assert nc_bfield.theta_data.ndim == 1
    assert nc_bfield.b_data.ndim == 2
    assert nc_bfield.db_dpsip_data.ndim == 2
    assert nc_bfield.db_dtheta_data.ndim == 2
    assert nc_bfield.shape == (len(nc_bfield.psip_data), len(nc_bfield.theta_data))
    assert nc_bfield.shape == nc_bfield.b_data.shape
    str(nc_bfield)


def test_nc_bfield_eval(nc_bfield: NcBfield):
    psip = 0.015
    theta = 1
    assert isinstance(nc_bfield.b(psip, theta), float)
    assert isinstance(nc_bfield.db_dpsip(psip, theta), float)
    assert isinstance(nc_bfield.db_dtheta(psip, theta), float)

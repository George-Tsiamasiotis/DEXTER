from dexter import Bfield


def test_fields(bfield: Bfield):
    assert isinstance(bfield.path, str)
    assert isinstance(bfield.typ, str)
    assert isinstance(bfield.psip_wall, float)

    assert bfield.psip_data.ndim == 1
    assert bfield.theta_data.ndim == 1
    assert bfield.b_data.ndim == 2
    assert bfield.db_dpsip_data.ndim == 2
    assert bfield.db_dtheta_data.ndim == 2

    assert bfield.shape == (len(bfield.psip_data), len(bfield.theta_data))
    assert bfield.shape == bfield.b_data.shape


def test_eval(bfield: Bfield):
    psip = 0.015
    theta = 1
    assert isinstance(bfield.b(psip, theta), float)
    assert isinstance(bfield.db_dpsip(psip, theta), float)
    assert isinstance(bfield.db_dtheta(psip, theta), float)
    assert isinstance(bfield.d2b_dpsip2(psip, theta), float)
    assert isinstance(bfield.d2b_dtheta2(psip, theta), float)
    assert isinstance(bfield.d2b_dpsip_dtheta(psip, theta), float)


def test_magic(bfield: Bfield):
    assert isinstance(bfield.__str__(), str)
    assert isinstance(bfield.__repr__(), str)

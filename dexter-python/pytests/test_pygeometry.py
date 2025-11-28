from dexter import Geometry


def test_fields(geometry: Geometry):
    assert isinstance(geometry.path, str)
    assert isinstance(geometry.typ1d, str)
    assert isinstance(geometry.typ2d, str)
    assert isinstance(geometry.baxis, float)
    assert isinstance(geometry.raxis, float)
    assert isinstance(geometry.psip_wall, float)
    assert isinstance(geometry.psi_wall, float)
    assert isinstance(geometry.r_wall, float)
    assert geometry.theta_data.ndim == 1
    assert geometry.psip_data.ndim == 1
    assert geometry.psi_data.ndim == 1
    assert geometry.r_data.ndim == 1
    assert geometry.rlab_data.ndim == 2
    assert geometry.zlab_data.ndim == 2

    assert geometry.shape == (len(geometry.psip_data), len(geometry.theta_data))
    assert geometry.shape == geometry.rlab_data.shape


def test_eval(geometry: Geometry):
    r = 0.1
    psip = 0.015
    theta = 3.14
    assert isinstance(geometry.r(psip), float)
    assert isinstance(geometry.psip(r), float)
    assert isinstance(geometry.rlab(psip, theta), float)
    assert isinstance(geometry.zlab(psip, theta), float)


def test_repr(geometry: Geometry):
    str(geometry)

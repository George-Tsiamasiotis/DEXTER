from dexter import NcGeometry


def test_nc_geometry_getters(nc_geometry: NcGeometry):
    assert isinstance(nc_geometry.path, str)
    assert isinstance(nc_geometry.typ1d, str)
    assert isinstance(nc_geometry.typ2d, str)
    assert isinstance(nc_geometry.baxis, float)
    assert isinstance(nc_geometry.raxis, float)
    assert isinstance(nc_geometry.zaxis, float)
    assert isinstance(nc_geometry.rgeo, float)
    assert isinstance(nc_geometry.psip_wall, float)
    assert isinstance(nc_geometry.psi_wall, float)
    assert isinstance(nc_geometry.r_wall, float)
    assert nc_geometry.theta_data.ndim == 1
    assert nc_geometry.psip_data.ndim == 1
    assert nc_geometry.psi_data.ndim == 1
    assert nc_geometry.r_data.ndim == 1
    assert nc_geometry.rlab_data.ndim == 2
    assert nc_geometry.zlab_data.ndim == 2
    assert nc_geometry.jacobian_data.ndim == 2

    assert nc_geometry.shape == (
        len(nc_geometry.psip_data),
        len(nc_geometry.theta_data),
    )
    assert nc_geometry.shape == nc_geometry.rlab_data.shape
    str(nc_geometry)


def test_nc_geometry_evals(nc_geometry: NcGeometry):
    r = 0.1
    psip = 0.015
    theta = 3.14
    assert isinstance(nc_geometry.r(psip), float)
    assert isinstance(nc_geometry.psi(psip), float)
    assert isinstance(nc_geometry.psip(r), float)
    assert isinstance(nc_geometry.rlab(psip, theta), float)
    assert isinstance(nc_geometry.zlab(psip, theta), float)
    assert isinstance(nc_geometry.jacobian(psip, theta), float)

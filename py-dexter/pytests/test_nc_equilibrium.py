from dexter import NcEquilibrium
from dexter import NcGeometry, NcQfactor, NcCurrent, NcBfield


def test_objects(
    toroidal_nc_equilibrium: NcEquilibrium,
    poloidal_nc_equilibrium: NcEquilibrium,
):
    geometry, qfactor, current, bfield = toroidal_nc_equilibrium.objects()
    assert isinstance(geometry, NcGeometry)
    assert isinstance(qfactor, NcQfactor)
    assert isinstance(current, NcCurrent)
    assert isinstance(bfield, NcBfield)
    geometry, qfactor, current, bfield = poloidal_nc_equilibrium.objects()
    assert isinstance(geometry, NcGeometry)
    assert isinstance(qfactor, NcQfactor)
    assert isinstance(current, NcCurrent)
    assert isinstance(bfield, NcBfield)


def test_plots(
    toroidal_nc_equilibrium: NcEquilibrium,
    poloidal_nc_equilibrium: NcEquilibrium,
):
    _test_plots(toroidal_nc_equilibrium)
    _test_plots(poloidal_nc_equilibrium)


def _test_plots(nc_equilibrium: NcEquilibrium):
    nc_equilibrium.plot_b()
    nc_equilibrium.plot_db()

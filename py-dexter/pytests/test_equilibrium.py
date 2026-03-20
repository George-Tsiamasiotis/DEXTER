import pytest
from dexter import Equilibrium
from dexter import (
    NcGeometry,
    NcQfactor,
    NcCurrent,
    NcBfield,
    Perturbation,
    LarGeometry,
    UnityQfactor,
    LarCurrent,
    LarBfield,
    CosHarmonic,
)


def test_numerical_equilibrium(
    toroidal_nc_equilibrium: Equilibrium,
    poloidal_nc_equilibrium: Equilibrium,
):
    eq = toroidal_nc_equilibrium
    assert isinstance(eq.geometry, NcGeometry)
    assert isinstance(eq.qfactor, NcQfactor)
    assert isinstance(eq.current, NcCurrent)
    assert isinstance(eq.bfield, NcBfield)
    assert isinstance(eq.psi_wall, float)
    assert isinstance(eq.psip_wall, float)

    eq = poloidal_nc_equilibrium
    assert isinstance(eq.geometry, NcGeometry)
    assert isinstance(eq.qfactor, NcQfactor)
    assert isinstance(eq.current, NcCurrent)
    assert isinstance(eq.bfield, NcBfield)
    assert isinstance(eq.psi_wall, float)
    assert isinstance(eq.psip_wall, float)


def test_minimum_equilibrium():
    equilibrium = Equilibrium(
        qfactor=UnityQfactor(),
        current=LarCurrent(),
        bfield=LarBfield(),
    )
    with pytest.raises(AttributeError):
        equilibrium.psi_wall
    with pytest.raises(AttributeError):
        equilibrium.psip_wall
    with pytest.raises(TypeError):
        equilibrium.plot_b()
    with pytest.raises(TypeError):
        equilibrium.plot_db()


def test_full_equilibrium():
    equilibrium = Equilibrium(
        geometry=LarGeometry(2, 1.75, 0.5),
        qfactor=UnityQfactor(),
        current=LarCurrent(),
        bfield=LarBfield(),
        perturbation=Perturbation(
            [
                CosHarmonic(8e-4, 3, 1, 0),
                CosHarmonic(8e-4, 5, 3, 0),
            ]
        ),
    )
    with pytest.raises(AttributeError):
        equilibrium.psi_wall
    with pytest.raises(AttributeError):
        equilibrium.psip_wall


def test_plots(
    toroidal_nc_equilibrium: Equilibrium,
    poloidal_nc_equilibrium: Equilibrium,
):
    _test_plots(toroidal_nc_equilibrium)
    _test_plots(poloidal_nc_equilibrium)


# FIXME: add analytical after vectorizing eval methods
def _test_plots(nc_equilibrium: Equilibrium):
    nc_equilibrium.plot_b()
    nc_equilibrium.plot_db()
    nc_equilibrium.plot_flux_surfaces()
    nc_equilibrium.plot_boozer_theta()

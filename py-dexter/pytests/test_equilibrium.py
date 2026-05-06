import pytest
from dexter import Equilibrium
from dexter import (
    LastClosedFluxSurface,
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
    assert isinstance(eq.psi_last, float)
    assert isinstance(eq.psip_last, float)
    assert isinstance(eq.rlast, float)
    assert isinstance(eq.raxis, float)
    assert isinstance(eq.baxis, float)

    eq = poloidal_nc_equilibrium
    assert isinstance(eq.geometry, NcGeometry)
    assert isinstance(eq.qfactor, NcQfactor)
    assert isinstance(eq.current, NcCurrent)
    assert isinstance(eq.bfield, NcBfield)
    assert isinstance(eq.psi_last, float)
    assert isinstance(eq.psip_last, float)
    assert isinstance(eq.rlast, float)
    assert isinstance(eq.raxis, float)
    assert isinstance(eq.baxis, float)


def test_minimum_equilibrium():
    LCFS = LastClosedFluxSurface("Toroidal", 0.1)
    equilibrium = Equilibrium(
        qfactor=UnityQfactor(LCFS),
        current=LarCurrent(),
        bfield=LarBfield(),
    )
    assert equilibrium.psi_last == 0.1
    assert equilibrium.psip_last == 0.1
    with pytest.raises(AttributeError):
        equilibrium.plot_b()
    with pytest.raises(AttributeError):
        equilibrium.plot_db()


def test_full_equilibrium():
    lcfs = LastClosedFluxSurface("Toroidal", 0.45)
    equilibrium = Equilibrium(
        geometry=LarGeometry(2, 1.75, 0.5),
        qfactor=UnityQfactor(lcfs),
        current=LarCurrent(),
        bfield=LarBfield(),
        perturbation=Perturbation(
            [
                CosHarmonic(8e-4, lcfs, 3, 1, 0),
                CosHarmonic(8e-4, lcfs, 5, 3, 0),
            ]
        ),
    )
    assert equilibrium.psi_last == 0.45
    assert equilibrium.psip_last == 0.45


def test_plots(
    toroidal_nc_equilibrium: Equilibrium,
    poloidal_nc_equilibrium: Equilibrium,
):
    _test_plots(toroidal_nc_equilibrium)
    _test_plots(poloidal_nc_equilibrium)


# FIXME: add analytical after vectorizing eval methods
def _test_plots(nc_equilibrium: Equilibrium):
    nc_equilibrium.plot_b(units="SI")
    nc_equilibrium.plot_b(units="NU")
    nc_equilibrium.plot_db()
    nc_equilibrium.plot_flux_surfaces()
    nc_equilibrium.plot_boozer_theta()
    nc_equilibrium.plot_midplane()

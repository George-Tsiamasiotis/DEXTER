from math import isclose

from dexter import (
    Equilibrium,
    LarCurrent,
    UnityQfactor,
    LarBfield,
    LarGeometry,
    LastClosedFluxSurface,
)


def test_gcmotion_units():
    geometry = LarGeometry(3.5, 1.75, 0.5)
    LCFS = LastClosedFluxSurface("Toroidal", geometry.psi_last)
    equilibrium = Equilibrium(
        geometry=geometry,
        qfactor=UnityQfactor(LCFS),
        current=LarCurrent(),
        bfield=LarBfield(),
        species="Deuterium",
    )

    assert isclose(equilibrium._frequency_unit, 167629580.2981652)
    assert isclose(equilibrium._energy_unit, 1.4393791168013255e-10)

    Q = equilibrium.quantity
    assert isclose(Q(1e-5, "energy_units").to("keV").m, 8.983897819104792)
    assert isclose(Q(0.003, "frequency_units").to("kilohertz").m, 502.88874089449564)

    assert isclose(Q(geometry.baxis, "Tesla").to("bfield_units").m, 1)
    assert isclose(Q(1, "bfield_units").to("Tesla").m, geometry.baxis)

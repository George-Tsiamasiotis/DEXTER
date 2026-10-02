import dexter as dex

geometry = dex.LarGeometry(baxis=2, raxis=1.75, rlast=0.5)
LCFS = geometry.psi_last
machine = dex.Machine(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.9, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

plane = dex.EnergyPzetaPlane(machine, 3e-5)
plane.show()

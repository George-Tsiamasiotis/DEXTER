import numpy as np
import dexter as dex

geometry = dex.LarGeometry(baxis=2, raxis=1.75, rlast=0.5)
LCFS = geometry.psi_last
machine = dex.Machine(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.9, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

mu = 8e-5

plane = dex.EnergyPzetaPlane(machine, mu)

pzeta = np.linspace(-1.9, 0.6, 1000) * machine.psip_last.value
plane.calculate_critical_points(pzeta)
plane.calculate_boundary_extrema(pzeta)
print(plane)
plane.show(xlim=(-1.9, 0.6), ylim=(0.4, 1.8), curves=True)

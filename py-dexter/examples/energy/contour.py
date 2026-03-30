import numpy as np
import dexter as dex
from math import pi as PI

geometry = dex.LarGeometry(4, 1.75, 0.5)
LCFS = dex.LastClosedFluxSurface("Toroidal", geometry.psi_last)
equilibrium = dex.Equilibrium(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(1.1, 3.2, LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

coms = dex.COMs(
    pzeta=-0.014,
    mu=2e-5,
)

psi_array = np.linspace(0, equilibrium.psi_last, 200)
theta_array = np.linspace(-PI, PI, 200)

energy_array = coms.energy_of_psi_grid(equilibrium, psi_array, theta_array)

dex.energy_contour(psi_array, theta_array, energy_array, levels=40)

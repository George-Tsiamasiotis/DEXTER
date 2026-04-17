"""Classification of all particles in a LAR equilibrium.

This script is a combination of all previous tests.
"""

import numpy as np
import dexter as dex
from math import sqrt
import matplotlib.pyplot as plt

LCFS = dex.LastClosedFluxSurface("Toroidal", 0.03)
raxis = 1.75
rlast = sqrt(2 * LCFS.value) * raxis  # `rlast` must be in [m]
equilibrium = dex.Equilibrium(
    geometry=dex.LarGeometry(baxis=1, raxis=raxis, rlast=rlast),
    qfactor=dex.ParabolicQfactor(1.1, 3.9, LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

mu = 6e-5

pzeta0 = (
    np.asarray([-0.8, -0.1, -1.3, -0.8, -0.5, -0.0448, 0, -0.4]) * equilibrium.psip_last
)
psi0 = np.asarray([0.001, 0.015, 0.015, 0.01, 0.018, 0.0045, 0.0014, 0.0014])
particle_count = len(pzeta0)
initial_conditions = dex.QueueInitialConditions.mixed(
    t0=np.zeros(particle_count),
    flux0=dex.InitialFluxArray("Toroidal", psi0),
    theta0=np.ones(particle_count),
    zeta0=np.zeros(particle_count),
    pzeta0=pzeta0,
    mu0=np.full(particle_count, mu),
)

queue = dex.Queue(initial_conditions)
queue.classify(equilibrium)

dex.plot_parabolas(equilibrium, mu, particles=queue.particles, xlim=(-2, 0.9))

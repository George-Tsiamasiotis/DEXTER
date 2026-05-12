"""Poincare map calculation on a 'θ=const' surface in an analytical equilibrium.

Also see './analytical_const_zeta.py', which tracks the same orbits but for a different cross section.
"""

import numpy as np
import dexter as dex
from math import pi as PI

geometry = dex.LarGeometry(4, 1.75, 0.5)
LCFS = dex.LastClosedFluxSurface("Toroidal", geometry.psi_last)
equilibrium = dex.Equilibrium(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(1.1, 4.2, LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
    perturbation=dex.Perturbation(
        [
            dex.CosHarmonic(8e-5, LCFS, 3, 1, PI),
            dex.CosHarmonic(3e-4, LCFS, 5, 4, PI),
        ]
    ),
)

particle_count = 50
rlast = np.sqrt(equilibrium.psi_last * 2)
r0s = np.linspace(0, rlast, particle_count)
psi0s = equilibrium.geometry.psi_of_r(r0s)

initial_conditions = dex.QueueInitialConditions.boozer(
    t0=np.zeros(particle_count),
    flux0=dex.InitialFluxArray("Toroidal", psi0s),
    theta0=np.zeros(particle_count),
    zeta0=np.zeros(particle_count),
    rho0=1e-7 * np.ones(particle_count),
    mu0=0 * np.ones(particle_count),
)

intersect_params = dex.IntersectParams("ConstTheta", 0.0, 500)
queue = dex.Queue(initial_conditions)

queue.intersect(
    equilibrium=equilibrium,
    intersect_params=intersect_params,
    max_steps=2_000_000,
    energy_rel_tol=1e-11,
)
queue.plot_const_theta_cartesian_poincare()

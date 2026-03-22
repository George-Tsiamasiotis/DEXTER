"""Poincare map calculation on a 'θ=const' surface in an analytical equilibrium.

Also see './analytical_const_zeta.py', which tracks the same orbits but for a different cross section.
"""

import numpy as np
import dexter as dex
from math import pi as PI

geometry = dex.LarGeometry(4, 1.75, 0.5)
equilibrium = dex.Equilibrium(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(1.1, 4.2, ("Toroidal", geometry.psi_wall)),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
    perturbation=dex.Perturbation(
        [
            dex.CosHarmonic(8e-5, 3, 1, PI),
            dex.CosHarmonic(8e-5, 5, 3, PI),
        ]
    ),
)

particle_count = 80
psips = np.linspace(0, 1, particle_count) * equilibrium.psip_wall
psis = np.asarray([equilibrium.qfactor.psi_of_psip(psip) for psip in psips])

initial_conditions = dex.QueueInitialConditions.boozer(
    t0=np.zeros(particle_count),
    flux0=dex.InitialFluxArray("Toroidal", psis),
    theta0=np.zeros(particle_count),
    zeta0=np.zeros(particle_count),
    rho0=1e-7 * np.ones(particle_count),
    mu0=0 * np.ones(particle_count),
)

intersect_params = dex.IntersectParams("ConstTheta", 0.0, 300)
queue = dex.Queue(initial_conditions)

queue.intersect(
    equilibrium=equilibrium,
    intersect_params=intersect_params,
    max_steps=2_000_000,
    energy_rel_tol=1e-11,
)
print(queue)
queue.plot_const_theta_cartesian_poincare()

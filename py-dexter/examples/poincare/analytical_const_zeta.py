"""Poincare map calculation on a 'ζ=const' surface in an analytical equilibrium.

Also see './analytical_const_theta.py', which tracks the same orbits but for a different cross section.
"""

import numpy as np
import dexter as dex
from math import pi as PI

flux_wall = ("Toroidal", 0.5)
qfactor = dex.ParabolicQfactor(1.1, 4.2, flux_wall)
current = dex.LarCurrent()
bfield = dex.LarBfield()
perturbation = dex.Perturbation(
    [
        dex.CosHarmonic(8e-4, 3, 1, PI),
        dex.CosHarmonic(8e-4, 5, 3, PI),
    ]
)

particle_count = 80
initial_conditions = dex.QueueInitialConditions(
    t0=np.zeros(particle_count),
    flux0=dex.InitialFluxArray1(
        "Toroidal",
        np.linspace(0, 1, particle_count) * qfactor.psi_wall,
    ),
    theta0=np.zeros(particle_count),
    zeta0=np.zeros(particle_count),
    rho0=1e-7 * np.ones(particle_count),
    mu0=0 * np.ones(particle_count),
)

intersect_params = dex.IntersectParams("ConstZeta", 0.0, 300)
queue = dex.Queue(initial_conditions)

queue.intersect(
    qfactor=qfactor,
    current=current,
    bfield=bfield,
    perturbation=perturbation,
    intersect_params=intersect_params,
    max_steps=1_000_000,
    energy_rel_tol=1e-11,
)
print(queue)
queue.plot_const_zeta_cartesian_poincare(initial=True)

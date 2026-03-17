"""Poincare map calculation in an analytical equilibrium."""

import numpy as np
import dexter as dex

flux_wall = ("Toroidal", 0.45)
qfactor = dex.ParabolicQfactor(1.1, 3.8, flux_wall)
current = dex.LarCurrent()
bfield = dex.LarBfield()
perturbation = dex.Perturbation(
    [
        dex.CosHarmonic(3e-3, 3, 1, 0),
        dex.CosHarmonic(2e-3, 7, 2, 0),
        dex.CosHarmonic(1e-3, 15, 4, 0),
    ]
)

particle_count = 100
initial_conditions = dex.QueueInitialConditions(
    t0=np.zeros(particle_count),
    flux0=dex.InitialFluxArray1(
        "Toroidal",
        np.linspace(0.01, 0.99, particle_count) * qfactor.psi_wall,
    ),
    theta0=np.zeros(particle_count),
    zeta0=np.zeros(particle_count),
    rho0=1e-5 * np.ones(particle_count),
    mu0=1e-6 * np.ones(particle_count),
)

intersect_params = dex.IntersectParams("ConstZeta", 0.0, 1000)
queue = dex.Queue(initial_conditions)

queue.intersect(
    qfactor=qfactor,
    current=current,
    bfield=bfield,
    perturbation=perturbation,
    intersect_params=intersect_params,
)
print(queue)

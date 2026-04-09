"""Calculation of ζ=0 intersections of a particle in an analytical equilibrium."""

import dexter as dex
from math import sqrt

LCFS = dex.LastClosedFluxSurface("Toroidal", 0.45)
equilibrium = dex.Equilibrium(
    geometry=dex.LarGeometry(baxis=2, raxis=1.75, rlast=sqrt(2 * LCFS.value)),
    qfactor=dex.ParabolicQfactor(1.1, 3.8, LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
    perturbation=dex.Perturbation(
        [
            dex.CosHarmonic(5e-3, LCFS, 3, 1, 0),
            dex.CosHarmonic(6e-3, LCFS, 7, 2, 0),
            dex.CosHarmonic(7e-3, LCFS, 15, 4, 0),
        ]
    ),
)

initial_conditions = dex.InitialConditions.boozer(
    t0=0,
    flux0=dex.InitialFlux("Toroidal", 0.03),
    theta0=1.0,
    zeta0=0.0,
    rho0=8e-3,
    mu0=1e-6,
)

particle = dex.Particle(initial_conditions)
intersect_params = dex.IntersectParams("ConstZeta", 0.0, 1000)

particle.intersect(
    equilibrium=equilibrium,
    intersect_params=intersect_params,
    stepping_method=("FixedStep", 0.4),
    max_steps=10_000_000,
)
print(particle)

particle.plot_evolution()
dex.plot_particle_poloidal_drift(particle, equilibrium, levels=100)

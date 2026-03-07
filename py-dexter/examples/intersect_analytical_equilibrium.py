"""Calculation of ζ=0 intersections of a particle in an analytical equilibrium."""

import dexter as dex

flux_wall = ("Toroidal", 0.45)
qfactor = dex.ParabolicQfactor(1.1, 3.8, flux_wall)
current = dex.LarCurrent()
bfield = dex.LarBfield()
perturbation = dex.perturbation(
    [
        dex.CosHarmonic(3e-3, 3, 1, 0),
        dex.CosHarmonic(2e-3, 7, 2, 0),
        dex.CosHarmonic(1e-3, 15, 4, 0),
    ]
)

initial_conditions = dex.InitialConditions(
    t0=0,
    flux0=("Toroidal", 0.02),
    theta0=0.0,
    zeta0=0.0,
    rho0=8e-3,
    mu0=1e-6,
)

particle = dex.Particle(initial_conditions)
intersect_params = dex.IntersectParams("ConstZeta", 0.0, 1000)

particle.intersect(
    qfactor=qfactor,
    current=current,
    bfield=bfield,
    perturbation=perturbation,
    intersect_params=intersect_params,
    stepping_method=("FixedStep", 0.4),
    max_steps=10_000_000,
)
print(particle)

particle.plot_poloidal_drift()

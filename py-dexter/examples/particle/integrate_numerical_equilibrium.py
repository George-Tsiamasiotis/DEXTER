"""Integration of a particle in a numerical equilibrium."""

import dexter as dex
from dexter.equilibrium import _TEST_NETCDF_PATH

path = _TEST_NETCDF_PATH

equilibrium = dex.numerical_equilibrium(path, "Steffen", "Bicubic")
equilibrium.perturbation = dex.Perturbation([])

initial_conditions = dex.InitialConditions.mixed(
    t0=0,
    flux0=dex.InitialFlux("Toroidal", 0.77 * equilibrium.psi_last),
    theta0=0.0,
    zeta0=0.0,
    pzeta0=-0.18,
    mu0=5e-3,
)

particle = dex.Particle(initial_conditions)

particle.integrate(
    equilibrium=equilibrium,
    teval=(0, 8e3),
)
print(particle)

particle.plot_evolution()
dex.plot_particle_poloidal_drift(particle, equilibrium)

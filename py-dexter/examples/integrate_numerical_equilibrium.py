"""Integration of a particle in a numerical equilibrium."""

import dexter as dex
from dexter.equilibrium import _TEST_NETCDF_PATH

path = _TEST_NETCDF_PATH

equilibrium = dex.numerical_equilibrium(path, "Steffen", "Bicubic")
equilibrium.perturbation = dex.Perturbation(
    [
        dex.NcHarmonic(
            path, interp_type="Steffen", m=2, n=1, phase_method="Interpolation"
        ),
        dex.NcHarmonic(
            path, interp_type="Steffen", m=3, n=2, phase_method="Interpolation"
        ),
    ]
)

initial_conditions = dex.InitialConditions(
    t0=0,
    flux0=dex.InitialFlux("Poloidal", 0.8 * equilibrium.psip_wall),
    theta0=1.0,
    zeta0=0.0,
    rho0=1e-3,
    mu0=0,
)

particle = dex.Particle(initial_conditions)

particle.integrate(
    equilibrium=equilibrium,
    teval=(0, 7e5),
)
print(particle)

particle.plot_evolution()

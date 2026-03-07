"""Integration of a particle in a numerical equilibrium."""

import dexter as dex

path = "./data.nc"

geometry = dex.NcGeometry(path, interp1d_type="Steffen", interp2d_type="Bicubic")
qfactor = dex.NcQfactor(path, interp_type="Steffen")
current = dex.NcCurrent(path, interp_type="Steffen")
bfield = dex.NcBfield(path, interp_type="Bicubic")
perturbation = dex.perturbation(
    [
        dex.NcHarmonic(path, interp_type="Steffen", m=m, n=n, phase_method="Zero")
        for n in range(1, 24)
        for m in range(-1, 3)
    ]
)

initial_conditions = dex.InitialConditions(
    t0=0,
    flux0=("Poloidal", 0.1 * qfactor.psip_wall),
    theta0=1.0,
    zeta0=0.0,
    rho0=1e-3,
    mu0=0,
)

particle = dex.Particle(initial_conditions)

particle.integrate(
    qfactor=qfactor,
    current=current,
    bfield=bfield,
    perturbation=perturbation,
    teval=(0, 6e3),
)
print(particle)

particle.plot_evolution()
particle.plot_poloidal_drift()

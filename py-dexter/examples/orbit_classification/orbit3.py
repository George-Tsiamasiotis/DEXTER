"""Classification of Cu-Passing particle in a LAR equilibrium.

CuPassing-Lost inside the right wall parabola, outside the left wall parabola and left from the
`Pζ/ψp_last = -1` line.

This script accompanies Rust's `orbit_classification` tests.
"""

import dexter as dex
from math import sqrt

LCFS = dex.LastClosedFluxSurface("Toroidal", 0.03)
raxis = 1.75
rlast = sqrt(2 * LCFS.value) * raxis  # `rlast` must be in [m]
equilibrium = dex.Equilibrium(
    geometry=dex.LarGeometry(baxis=1, raxis=raxis, rlast=rlast),
    qfactor=dex.ParabolicQfactor(1.1, 3.9, LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

pzeta = -equilibrium.psip_last * 1.3
mu = 6e-5

initial_conditions = dex.InitialConditions.mixed(
    t0=0,
    flux0=dex.InitialFlux("Toroidal", 0.015),
    theta0=1.0,
    zeta0=0.0,
    pzeta0=pzeta,
    mu0=mu,
)

particle = dex.Particle(initial_conditions)
particle.close(equilibrium=equilibrium)
particle.classify(equilibrium=equilibrium)
assert particle.orbit_type == "CuPassingLost"
energy_pzeta_point = (
    pzeta / equilibrium.psip_last,
    particle.initial_energy / mu,
)
print(particle)

# =========================

fig, ax = dex.plot_parabolas(equilibrium, mu, particles=[particle], show=False)
ax.scatter(*energy_pzeta_point, c="r", label=rf"${particle.orbit_type}\ Particle$")
ax.legend()

# particle.plot_evolution()
dex.plot_particle_poloidal_drift(particle, equilibrium, locator="Log")

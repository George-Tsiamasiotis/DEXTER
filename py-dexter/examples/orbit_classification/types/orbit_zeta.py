"""Classification of a particle in a LAR equilibrium.

CoPassing-Lost inside the right wall parabola and inside axis parabola.

This script accompanies Rust's `orbit_classification` tests.
"""

import dexter as dex
from math import sqrt, pi

LCFS = dex.LastClosedFluxSurface("Toroidal", 0.03)
raxis = 1.75
rlast = sqrt(2 * LCFS.value) * raxis  # `rlast` must be in [m]
equilibrium = dex.Equilibrium(
    geometry=dex.LarGeometry(baxis=1, raxis=raxis, rlast=rlast),
    qfactor=dex.ParabolicQfactor(1.1, 3.9, LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

pzeta = -equilibrium.psip_last * 0.4
mu = 6e-5

initial_conditions = dex.InitialConditions.mixed(
    t0=0,
    flux0=dex.InitialFlux("Toroidal", 0.025),
    theta0=pi,
    zeta0=0.0,
    pzeta0=pzeta,
    mu0=mu,
)

particle = dex.Particle(initial_conditions)
particle.close(equilibrium=equilibrium)
particle.classify(equilibrium=equilibrium)
assert particle.energy_pzeta_position == "Zeta"
assert particle.orbit_type == "CoPassingLost"
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

"""Classification of a particle in a LAR machine.

CuPassing-Confined inside the left wall parabola.

This script accompanies Rust's `orbit_classification` tests.
"""

import dexter as dex
from math import sqrt, pi

LCFS = dex.MagneticFlux.Toroidal(0.03)
raxis = 1.75
rlast = sqrt(2 * LCFS.value) * raxis  # `rlast` must be in [m]
machine = dex.Machine(
    geometry=dex.LarGeometry(baxis=1, raxis=raxis, rlast=rlast),
    qfactor=dex.ParabolicQfactor(1.1, 3.9, LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

pzeta = -machine.psip_last * 0.4
mu = 6e-5

initial_conditions = dex.InitialConditions.Mixed(
    t0=0,
    flux0=dex.MagneticFlux.Toroidal(0.025),
    theta0=pi,
    zeta0=0.0,
    pzeta0=pzeta,
    mu0=mu,
)

particle = dex.Particle(initial_conditions)
particle.close(machine=machine)
particle.classify(machine=machine)
assert particle.energy_pzeta_position == "Zeta"
assert particle.orbit_type == "CoPassingLost"
print(particle)

# =========================

dex.plot_poloidal_drift(machine, particle, show=False)

plane = dex.EnergyPzetaPlane(machine, mu)
plane.add_particles(particle)
plane.show(particles=True)

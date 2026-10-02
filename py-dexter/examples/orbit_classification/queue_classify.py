"""Orbit classification for particles with `ψ=const` and constant `ρ` and `μ` in an analytical machine."""

import numpy as np
import dexter as dex
from math import pi as PI

RNG = np.random.default_rng(42)

# Equilibrium setup
geometry = dex.LarGeometry(1, 1.75, 0.5)
LCFS = geometry.psi_last
machine = dex.Machine(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.9, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
    perturbation=dex.Perturbation([]),
)

mu = 7e-6

# Initial Conditions setup
num = 50000
psi0s = dex.MagneticFluxArray.Toroidal(RNG.random(num) * LCFS.value)

initial_conditions = dex.QueueInitialConditions.Mixed(
    t0=np.zeros(num),
    flux0=psi0s,
    theta0=2 * PI * RNG.random(num),
    zeta0=np.zeros(num),
    pzeta0=np.linspace(-1.4, 0.2, num) * machine.psip_last.value,
    mu0=np.full(num, mu),
)

# Queue setup
queue = dex.Queue(initial_conditions)

# Run
queue.classify(machine)
queue.retain_energy(0, 3 * mu)

# =========================

# Plot orbits on the E-Pζ space
plane = dex.EnergyPzetaPlane(machine, mu)
plane.add_particles(queue)
plane.show()
plane.clear_particles()

# =========================

confined_particles = queue.retain_orbit_types(
    [
        "CoPassingConfined",
        "CuPassingConfined",
        "TrappedConfined",
        "Stagnated",
        "Potato",
    ]
)

plane.add_particles(queue)
plane.show()
plane.clear_particles()

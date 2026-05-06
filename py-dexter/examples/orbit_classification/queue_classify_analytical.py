"""Orbit classification for particles with `ψ=const` and constant `ρ` and `μ` in an analytical equilibrium."""

import numpy as np
import dexter as dex
from math import pi as PI
import matplotlib.pyplot as plt

# Equilibrium setup
geometry = dex.LarGeometry(1, 1.75, 0.5)
LCFS = dex.LastClosedFluxSurface("Toroidal", geometry.psi_last)
equilibrium = dex.Equilibrium(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.9, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
    perturbation=dex.Perturbation([]),
)

mu = 7e-6

# Initial Conditions setup
num = 100000
psi0s = dex.InitialFluxArray("Toroidal", np.random.random(num) * LCFS.value)

initial_conditions = dex.QueueInitialConditions.mixed(
    t0=np.zeros(num),
    flux0=psi0s,
    theta0=2 * PI * np.random.random(num),
    zeta0=np.zeros(num),
    pzeta0=np.linspace(-1.4, 0.2, num) * equilibrium.psip_last,
    mu0=np.full(num, mu),
)

# Queue setup
queue = dex.Queue(initial_conditions)

# Run
queue.classify(equilibrium)

# =========================

# Plot orbits on the E-Pζ space
dex.plot_parabolas(equilibrium, mu, particles=queue.particles, ymax=3)

# =========================

trapped_particles = [
    particle
    for particle in queue.particles
    if particle.orbit_type in ["TrappedConfined", "TrappedLost"]
]
_, ax = dex.plot_parabolas(equilibrium, mu, particles=trapped_particles, show=False)
ax.set_title("Trapped only")
plt.show()

# =========================

cupassing_confined_particles = [
    particle for particle in queue.particles if "Passing" in particle.orbit_type
]
_, ax = dex.plot_parabolas(
    equilibrium, mu, particles=cupassing_confined_particles, show=False
)
ax.set_title("All passing particles")
plt.show()

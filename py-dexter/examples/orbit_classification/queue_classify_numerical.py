"""Orbit classification for particles with `ψ=const` and constant `ρ` and `μ` in an analytical equilibrium."""

import numpy as np
import dexter as dex
from math import pi as PI
import matplotlib.pyplot as plt

# Equilibrium setup
equilibrium = dex.numerical_equilibrium("./data.nc", "Steffen", "Bicubic")

mu = 1e-3

# Initial Conditions setup
num = 10000
psip0s = dex.InitialFluxArray("Poloidal", np.random.random(num) * equilibrium.psip_last)

initial_conditions = dex.QueueInitialConditions.mixed(
    t0=np.zeros(num),
    flux0=psip0s,
    theta0=2 * PI * np.random.random(num),
    zeta0=np.zeros(num),
    pzeta0=np.linspace(-1.05, 0.5, num) * equilibrium.psip_last,
    mu0=np.full(num, mu),
)

# Queue setup
queue = dex.Queue(initial_conditions)

# Run
queue.classify(equilibrium)

# =========================

# Plot orbits on the E-Pζ space
xlim = (-1.1, 0.5)
dex.plot_parabolas(equilibrium, mu, particles=queue.particles, ymax=5, xlim=xlim)

# =========================

trapped_particles = [
    particle
    for particle in queue.particles
    if particle.orbit_type in ["TrappedConfined", "TrappedLost"]
]
_, ax = dex.plot_parabolas(
    equilibrium, mu, particles=trapped_particles, show=False, xlim=xlim
)
ax.set_title("Trapped only")
plt.show()

# =========================

unclassified_particles = [
    particle for particle in queue.particles if particle.orbit_type == "Unclassified"
]
_, ax = dex.plot_parabolas(
    equilibrium, mu, particles=unclassified_particles, show=False, xlim=xlim
)
ax.set_title("Unclassified only")
plt.show()

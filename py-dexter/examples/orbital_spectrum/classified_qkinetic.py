"""Selects a specific family of particles in the COM space and calculates their frequencies."""

import numpy as np
import dexter as dex
import matplotlib.pyplot as plt

# Equilibrium setup
LCFS = dex.LastClosedFluxSurface(kind="Toroidal", value=0.05)
equilibrium = dex.Equilibrium(
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=4.9, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
    perturbation=dex.Perturbation([]),
)

mu = 7e-6

# Initial Conditions setup
num = 50000
psi0s = dex.InitialFluxArray("Toroidal", np.random.random(num) * LCFS.value)

initial_conditions = dex.QueueInitialConditions.mixed(
    t0=np.zeros(num),
    flux0=psi0s,
    theta0=np.zeros(num),
    zeta0=np.zeros(num),
    pzeta0=np.linspace(-1.35, 0.3, num) * equilibrium.psip_last,
    mu0=np.full(num, mu),
)

# Queue setup
queue = dex.Queue(initial_conditions)
queue.classify(equilibrium)

# =========================

# Plot all orbits on the E-Pζ space and calculate all frequencies
_, ax = dex.plot_parabolas(equilibrium, mu, particles=queue.particles, show=False)
ax.set_title("All discovered particles")
plt.show()

# Filter out lost particles that are going to escape
confined_queue = dex.Queue.from_particles(
    [particle for particle in queue.particles if "Lost" not in particle.orbit_type]
)
_, ax = dex.plot_parabolas(
    equilibrium, mu, particles=confined_queue.particles, show=False
)
ax.set_title(r"All non-Lost particles")
plt.show()
confined_queue.close(equilibrium, max_steps=100000)
confined_queue.plot_qkinetic_pzeta_sweep(psip_last=equilibrium.psip_last)
confined_queue.plot_qkinetic_energy_sweep()

# =========================

# Calculate the frequencies of the trapped particles only
trapped_particles = [
    particle
    for particle in queue.particles
    if particle.orbit_type in ["TrappedConfined", "Stagnated", "Potato"]
]
_, ax = dex.plot_parabolas(equilibrium, mu, particles=trapped_particles, show=False)
ax.set_title("Trapped particles only")
plt.show()

trapped_queue = dex.Queue.from_particles(trapped_particles)
trapped_queue.close(equilibrium)

trapped_queue.plot_qkinetic_pzeta_sweep(psip_last=equilibrium.psip_last)
trapped_queue.plot_qkinetic_energy_sweep()

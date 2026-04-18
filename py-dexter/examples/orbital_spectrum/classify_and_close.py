"""Selects a specific family of particles in the COM space and calculates their frequencies."""

import numpy as np
import dexter as dex
import matplotlib.pyplot as plt

# Equilibrium setup
LCFS = dex.LastClosedFluxSurface(kind="Toroidal", value=0.05)
equilibrium = dex.Equilibrium(
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.9, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
    perturbation=dex.Perturbation([]),
)

mu = 7e-6

# Initial Conditions setup
num = 20000
psi0s = dex.InitialFluxArray("Toroidal", np.random.random(num) * LCFS.value)

initial_conditions = dex.QueueInitialConditions.mixed(
    t0=np.zeros(num),
    flux0=psi0s,
    theta0=np.zeros(num),
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
_, ax = dex.plot_parabolas(equilibrium, mu, particles=queue.particles, show=False)
ax.set_title("All discovered particles")
plt.show()

queue.close(equilibrium, max_steps=100000)
queue.plot_qkinetic_pzeta_sweep()

# =========================

trapped_particles = [
    particle
    for particle in queue.particles
    if particle.orbit_type in ["TrappedConfined", "TrappedLost"]
]
_, ax = dex.plot_parabolas(equilibrium, mu, particles=trapped_particles, show=False)
ax.set_title("Trapped particles only")
plt.show()

close_queue = dex.Queue.from_particles(trapped_particles)
close_queue.close(equilibrium)

close_queue.plot_qkinetic_pzeta_sweep()

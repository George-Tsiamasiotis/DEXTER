"""Filters orbits in COM space and calculates their frequencies."""

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

mu = 3e-5

# Initial Conditions setup
num = 50000
psi0s = dex.InitialFluxArray("Toroidal", np.random.random(num) * LCFS.value)

initial_conditions = dex.QueueInitialConditions.mixed(
    t0=np.zeros(num),
    flux0=psi0s,
    theta0=np.zeros(num),
    zeta0=np.zeros(num),
    pzeta0=np.linspace(-1.4, 0.2, num) * equilibrium.psip_last,
    mu0=np.full(num, mu),
)

# Queue setup and classification
queue = dex.Queue(initial_conditions)
queue.classify(equilibrium)


max_energy = 1.8
filtered_particles = [
    p
    for p in queue.particles
    if "Lost" not in p.orbit_type and p.initial_energy / mu < max_energy
]
filtered_queue = dex.Queue.from_particles(filtered_particles)

dex.plot_parabolas(
    equilibrium,
    mu,
    particles=filtered_queue.particles,
    xlim=(-1.4, 0.2),
    ymax=max_energy,
)

filtered_queue.close(equilibrium, max_steps=100000)
_, ax = filtered_queue.plot_qkinetic_pzeta_energy3d(
    psip_last=equilibrium.psip_last,
    show=False,
)
plt.show()

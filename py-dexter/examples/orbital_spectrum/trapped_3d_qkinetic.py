"""Selects TrappedConfined, Potato and Stagnated orbits in COM space and calculates their frequencies."""

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

mu = 1e-4

# Initial Conditions setup
num = 4_000_000
psi0s = dex.InitialFluxArray("Toroidal", np.random.random(num) * LCFS.value)

initial_conditions = dex.QueueInitialConditions.mixed(
    t0=np.zeros(num),
    flux0=psi0s,
    theta0=np.zeros(num),
    zeta0=np.zeros(num),
    pzeta0=np.linspace(-1.0, 0.2, num) * equilibrium.psip_last,
    mu0=np.full(num, mu),
)

# Queue setup and classification
queue = dex.Queue(initial_conditions)
queue.classify(equilibrium)

# =========================

# Filter particles and try to only keep a similar amount of each family, to make the plot prettier
particles = queue.particles
trapped_particles = [p for p in particles if p.orbit_type == "TrappedConfined"][::50]
potato_particles = [p for p in particles if p.orbit_type == "Potato"]
stagnated_particles = [p for p in particles if p.orbit_type == "Stagnated"]
filtered_particles = potato_particles + trapped_particles + stagnated_particles[::20]

filtered_queue = dex.Queue.from_particles(filtered_particles)
_, ax = dex.plot_parabolas(
    equilibrium,
    mu,
    particles=filtered_queue.particles,
    xlim=(-1.1, 0.2),
    show=False,
)
ax.set_title(r"All Potato and Stagnated particles")
ax.set_ybound(lower=0.6, upper=2)
plt.show()

# Frequency calculation
filtered_queue.close(equilibrium, max_steps=100000)

filtered_queue.plot_qkinetic_pzeta_sweep(psip_last=equilibrium.psip_last)
filtered_queue.plot_qkinetic_energy_sweep()

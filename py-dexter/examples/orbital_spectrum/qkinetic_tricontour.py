"""Selects a specific family of particles in the COM space and calculates their frequencies."""

import numpy as np
import dexter as dex

from dexter.types import EnergyPzetaPosition

RNG = np.random.default_rng(42)

# Equilibrium setup
geometry = dex.LarGeometry(1, 1.75, 0.5)
LCFS = dex.LastClosedFluxSurface("Toroidal", geometry.psi_last)
equilibrium = dex.Equilibrium(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=4.9, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
    perturbation=dex.Perturbation([]),
)

mu = 7e-6
ymax = 2.51

# Initial Conditions setup
num = 100000
psi0s = dex.InitialFluxArray("Toroidal", RNG.random(num) * LCFS.value)

initial_conditions = dex.QueueInitialConditions.mixed(
    t0=np.zeros(num),
    flux0=psi0s,
    theta0=2 * np.pi * RNG.random(num),
    zeta0=np.zeros(num),
    pzeta0=np.linspace(-1.18, 0.25, num) * equilibrium.psip_last,
    mu0=np.full(num, mu),
)

# Queue setup
queue = dex.Queue(initial_conditions)
queue.classify(equilibrium)

# =========================

dex.plot_parabolas(equilibrium, mu, queue.particles)

# Keep only specific E-Pζ regions
filter: list[EnergyPzetaPosition] = ["Alpha", "Epsilon", "Kappa", "Iota"]
filtered_queue = dex.Queue.from_particles(
    [
        particle
        for particle in queue.particles
        if particle.initial_energy / mu < ymax + 0.1
        and particle.energy_pzeta_position in filter
    ]
)

dex.plot_parabolas(equilibrium, mu, filtered_queue.particles, ymax=ymax)


filtered_queue.close(
    equilibrium,
    max_steps=100000,
    energy_rel_tol=1e-13,
    energy_abs_tol=1e-14,
)

dex.plot_qkinetic_tricontour(
    equilibrium,
    queue=filtered_queue,
    ymax=2.5,
)

dex.plot_qkinetic_tricontour(
    equilibrium,
    queue=filtered_queue,
    ymax=2.5,
    levels=np.arange(-4, 4, 0.5),
)

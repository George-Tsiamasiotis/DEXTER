"""Selects a specific family of particles in the COM space and calculates their frequencies."""

import numpy as np
import dexter as dex

from dexter.types import EnergyPzetaPosition

RNG = np.random.default_rng(42)

# Equilibrium setup
LCFS = dex.LastClosedFluxSurface(kind="Toroidal", value=0.04)
equilibrium = dex.Equilibrium(
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=4.9, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
    perturbation=dex.Perturbation([]),
)

mu = 7e-6

# Initial Conditions setup
num = 100000
psi0s = dex.InitialFluxArray("Toroidal", RNG.random(num) * LCFS.value)

initial_conditions = dex.QueueInitialConditions.mixed(
    t0=np.zeros(num),
    flux0=psi0s,
    theta0=2 * np.pi * RNG.random(num),
    zeta0=np.zeros(num),
    pzeta0=np.linspace(-1.16, -0.03, num) * equilibrium.psip_last,
    mu0=np.full(num, mu),
)

# Queue setup
queue = dex.Queue(initial_conditions)
queue.classify(equilibrium)

# =========================

# Keep only specific E-Pζ regions
filter: list[EnergyPzetaPosition] = ["Alpha", "Zeta", "Mu"]
filtered_queue = dex.Queue.from_particles(
    [
        particle
        for particle in queue.particles
        if particle.initial_energy / mu < 2.51
        and particle.energy_pzeta_position in filter
    ]
)

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
    clabel=False,
)

dex.plot_qkinetic_tricontour(
    equilibrium,
    queue=filtered_queue,
    ymax=2.5,
    levels=[-3, -2.5, -2, -1.5, -0.5, 0, 0.5],
)

"""qkinetic calculation for particles with `ψ=0...ψlast` and constant `ρ` and `μ` in a numerical equilibrium."""

import numpy as np
import dexter as dex

# Equilibrium setup
equilibrium = dex.numerical_equilibrium("data.nc", "Steffen", "Bicubic")

# Initial Conditions setup
num = 1000
psips = dex.InitialFluxArray("Poloidal", np.full(num, 0.2 * equilibrium.psip_last))

initial_conditions = dex.QueueInitialConditions.mixed(
    t0=np.zeros(num),
    flux0=psips,
    theta0=np.zeros(num),
    zeta0=np.zeros(num),
    pzeta0=np.linspace(-0.7, 0.5, num) * equilibrium.psip_last,
    mu0=np.full(num, 2e-4),
)

# Queue setup
queue = dex.Queue(initial_conditions)

# Run
queue.close(
    equilibrium=equilibrium,
    max_steps=100_000,
    energy_rel_tol=1e-11,
)
queue.classify(equilibrium)
queue.plot_qkinetic_pzeta_sweep(psip_last=equilibrium.psip_last)

"""qkinetic calculation for particles with `ψ=const` and constant `ρ` and `μ` in a numerical equilibrium."""

import numpy as np
import dexter as dex

# Equilibrium setup
equilibrium = dex.numerical_equilibrium("data.nc", "Steffen", "Bicubic")

# Initial Conditions setup
num = 10000
psip0s = dex.InitialFluxArray("Poloidal", np.full(num, 0.25 * equilibrium.psip_last))

initial_conditions = dex.QueueInitialConditions.mixed(
    t0=np.zeros(num),
    flux0=psip0s,
    theta0=np.zeros(num),
    zeta0=np.zeros(num),
    pzeta0=np.linspace(-0.3, 0.25, num),
    mu0=np.full(num, 5e-5),
)

# Queue setup
queue = dex.Queue(initial_conditions)

# Run
queue.close(equilibrium=equilibrium, max_steps=50_000)
queue.plot_qkinetic_energy_sweep()

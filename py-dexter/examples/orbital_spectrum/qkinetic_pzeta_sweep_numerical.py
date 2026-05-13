"""qkinetic calculation for particles with `ψ=0...ψlast` and constant `ρ` and `μ` in a numerical equilibrium."""

import numpy as np
import dexter as dex
import matplotlib.pyplot as plt

# Equilibrium setup
equilibrium = dex.numerical_equilibrium("data.nc", "Steffen", "Bicubic")

# Initial Conditions setup
num = 10000
psip0s = dex.InitialFluxArray("Poloidal", np.full(num, 0.25 * equilibrium.psip_last))

initial_conditions = dex.QueueInitialConditions.boozer(
    t0=np.zeros(num),
    flux0=psip0s,
    theta0=np.zeros(num),
    zeta0=np.zeros(num),
    rho0=np.linspace(1e-4, 1e-2, num),
    mu0=np.full(num, 5e-5),
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
_, ax = queue.plot_qkinetic_pzeta_sweep(psip_last=equilibrium.psip_last, show=False)
plt.show()

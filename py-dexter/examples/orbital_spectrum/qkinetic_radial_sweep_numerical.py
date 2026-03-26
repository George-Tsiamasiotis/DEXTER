"""qkinetic calculation for particles with `ψ=0...ψlast` and constant `ρ` and `μ` in a numerical equilibrium"""

import numpy as np
import dexter as dex
import matplotlib.pyplot as plt

# Equilibrium setup
equilibrium = dex.numerical_equilibrium("data.nc", "Steffen", "Bicubic")

# Initial Conditions setup
num = 4000
psips = dex.InitialFluxArray("Poloidal", np.linspace(0, equilibrium.psip_last, num))

initial_conditions = dex.QueueInitialConditions.boozer(
    t0=np.zeros(num),
    flux0=psips,
    theta0=np.zeros(num),
    zeta0=np.zeros(num),
    rho0=np.full(num, 9e-3),
    mu0=np.full(num, 1e-4),
)

# Queue setup
queue = dex.Queue(initial_conditions)

# Run
queue.close(equilibrium=equilibrium)
queue.plot_qkinetic_radial_sweep(show=False)

# Plot reversal
ax = queue.ax
ax.axvline(x=0.091818499, c="r", linewidth=2)

plt.show()

"""qkinetic calculation for particles with `ψ=0...ψlast` and constant `ρ` and `μ` in an analytical equilibrium"""

import numpy as np
import dexter as dex

# Equilibrium setup
LCFS = dex.LastClosedFluxSurface(kind="Toroidal", value=0.05)
equilibrium = dex.Equilibrium(
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.5, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
    perturbation=dex.Perturbation([]),
)

# Initial Conditions setup
num = 4000
psi0s = dex.InitialFluxArray("Toroidal", np.linspace(0, LCFS.value, num))

initial_conditions = dex.QueueInitialConditions.boozer(
    t0=np.zeros(num),
    flux0=psi0s,
    theta0=np.zeros(num),
    zeta0=np.zeros(num),
    rho0=np.full(num, 1e-2),
    mu0=np.full(num, 1e-4),
)

# Queue setup
queue = dex.Queue(initial_conditions)

# Run
queue.close(
    equilibrium=equilibrium,
    max_steps=100_000,
)
queue.classify(equilibrium)
queue.plot_qkinetic_radial_sweep(lcfs_value=LCFS.value)

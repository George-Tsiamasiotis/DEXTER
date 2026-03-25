"""qkinetic calculation for particles with `ψ=0...ψlast` and constant `ρ` and `μ` in an analytical equilibrium."""

import numpy as np
import dexter as dex

# Equilibrium setup
LCFS = dex.LastClosedFluxSurface(kind="Toroidal", value=0.02)
equilibrium = dex.Equilibrium(
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.5, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
    perturbation=dex.Perturbation([]),
)

# Initial Conditions setup
num = 4000
psis = dex.InitialFluxArray("Toroidal", np.full(num, 0.5 * LCFS.value))

initial_conditions = dex.QueueInitialConditions.mixed(
    t0=np.zeros(num),
    flux0=psis,
    theta0=np.zeros(num),
    zeta0=np.zeros(num),
    pzeta0=np.linspace(-0.014, -0.002, num),
    mu0=np.full(num, 7e-6),
)

# Queue setup
queue = dex.Queue(initial_conditions)

# Run
queue.close(equilibrium=equilibrium)
queue.plot_qkinetic_pzeta_sweep()

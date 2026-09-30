import numpy as np
import dexter as dex
import matplotlib.pyplot as plt

geometry = dex.LarGeometry(1, 1.75, 0.5)
LCFS = geometry.psi_last
machine = dex.Machine(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=6.5, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
    perturbation=dex.Perturbation(
        [
            dex.FluteMode(7e-5, LCFS, 5, 4, 0),
        ]
    ),
)

num = 30
mu = 3e-5
energy = machine.quantity(20, "keV").to("NormJoule").m

rs = np.linspace(0.03, geometry.rlast * 0.999, num)
psis = (rs / geometry.raxis) ** 2 / 2
theta0s = np.full(num, np.pi)

bs = machine.bfield.eval_b(psi=psis, theta=theta0s)
rho0s = -np.sqrt(2 * energy - 2 * mu * bs) / bs


initial_conditions = dex.QueueInitialConditions.Boozer(
    t0=np.zeros(num),
    flux0=dex.MagneticFluxArray.Toroidal(psis),
    theta0=theta0s,
    zeta0=np.zeros(num),
    rho0=rho0s,
    mu0=np.full(num, 1e-5),
)
queue = dex.Queue(initial_conditions)


intersect_params = dex.IntersectParams("ConstZeta", 0, 2000, "Initial")
solver_params = dex.SolverParams(
    method=dex.SteppingMethod.EnergyAdaptiveStep(1e-8, 1e-10),
    max_steps=10_000_000,
)
queue.intersect(
    machine,
    intersect_params,
    solver_params,
)
print(f"Min energy variance: {np.nanmin(queue.energy_rsd_array)}")
print(f"Max energy variance: {np.nanmax(queue.energy_rsd_array)}")

dex.plot_pzeta_poincare(machine, queue, intersect_params, initial=True)
dex.plot_rz_poincare(machine, queue, intersect_params, initial=True)

# =========================

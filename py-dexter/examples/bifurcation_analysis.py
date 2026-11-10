import numpy as np
import dexter as dex

geometry = dex.LarGeometry(2, 1.75, 0.5)
LCFS = geometry.psit_last
machine = dex.Machine(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=25, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

strength = machine.quantity(500, "kiloVolt / meter").to("NormVolt / NormMeter").m
rmean = 0.7 * geometry.rlast
mean = dex.MagneticFlux.Toroidal(float(geometry.eval_psit_of_r(rmean)))
waist = LCFS.value * 0.1
efield = dex.GaussianEfield(strength, mean, waist)
machine.efield = efield

mu = 7e-6

# ==========================

num1 = 10
num2 = 10
num3 = 9
r0s = np.concat(
    (
        np.linspace(0.31, 0.49, num1),
        np.linspace(0.04, 0.28, num2),
        np.linspace(0.29, 0.49, num3),
    ),
)
theta0s = np.concat(
    (
        np.full(num1, 0),
        np.full(num2, np.pi),
        np.full(num3, np.pi),
    ),
)
psit0s = geometry.eval_psit_of_r(r0s)
num = len(psit0s)

initial_conditions = dex.QueueInitialConditions.Mixed(
    t0=np.zeros(num),
    flux0=dex.MagneticFluxArray.Toroidal(psit0s),
    theta0=theta0s,
    zeta0=np.zeros(num),
    pzeta0=np.full(num, -1.15 * machine.psip_last.value),
    mu0=np.full(num, mu),
)

queue = dex.Queue(initial_conditions)
solver_params = dex.SolverParams(max_steps=100_000)
queue.close(machine, solver_params=solver_params, discard_arrays=False)


# ====================

plane = dex.EnergyPzetaPlane(machine, mu)
pzeta_lim = (-1.5, 0.4)
pzeta_array = np.linspace(*pzeta_lim, 3000) * machine.psip_last.value
plane.bifurcation_analysis(pzeta_array)
plane.add_particles(queue)
plane.print_durations()
fig, ax = plane.show(
    xlim=pzeta_lim, ylim=(0.0, 5), colors="Cycle", curves=True, show=False
)
fig.set_figheight(4)
fig.set_figwidth(6)

dex.plot_poloidal_drift(machine, queue)

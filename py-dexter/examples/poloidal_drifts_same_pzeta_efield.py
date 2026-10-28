import numpy as np
import dexter as dex
import matplotlib.pyplot as plt

geometry = dex.LarGeometry(1, 1.65, 0.5)
LCFS = geometry.psit_last
machine = dex.Machine(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.5, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

strength = machine.quantity(750, "kiloVolt / meter").to("NormVolt / NormMeter").m
rmean = 0.7 * geometry.rlast
mean = dex.MagneticFlux.Toroidal(float(geometry.eval_psit_of_r(rmean)))
waist = LCFS.value * 0.1
efield = dex.GaussianEfield(strength, mean, waist)
machine.efield = efield

# ==========================

mu = 2e-4

num1 = 11
num2 = 17
num3 = 9
r0s = np.concat(
    (
        np.linspace(0.22, 0.49, num1),
        np.linspace(0.04, 0.34, num2),
        np.linspace(0.35, 0.49, num3),
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
    pzeta0=np.full(num, -0.85 * machine.psip_last.value),
    mu0=np.full(num, mu),
)

queue = dex.Queue(initial_conditions)


queue.close(machine, discard_arrays=False)
dex.plot_poloidal_drift(machine, queue, colors="Cycle", show=False)

plane = dex.EnergyPzetaPlane(machine, mu)
pzeta_array = np.linspace(-2.2, 0.8, 1000) * machine.psip_last.value
plane.calculate_critical_points(pzeta_array)
plane.calculate_boundary_extrema(pzeta_array)
plane.add_particles(queue)
fig, ax = plane.show(
    xlim=(-2.2, 0.8), ylim=(0.6, 1.8), colors="Cycle", curves=True, show=False
)
fig.set_figheight(4)
fig.set_figwidth(6)
ax.legend().remove()

plt.show()

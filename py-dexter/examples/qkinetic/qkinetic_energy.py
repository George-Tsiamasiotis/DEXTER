import numpy as np
import dexter as dex

geometry = dex.LarGeometry(1, 1.75, 0.5)
LCFS = geometry.psit_last
machine = dex.Machine(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.5, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

# =========================

mu = 6e-5
pzeta = -0.3 * machine.psip_last.value

num = 1000
psit0_values = np.linspace(1e-6, LCFS.value, num)
psit0s = dex.MagneticFluxArray.Toroidal(psit0_values)
theta0s = np.zeros(num)

initial_conditions = dex.QueueInitialConditions.Mixed(
    t0=np.zeros(num),
    flux0=psit0s,
    theta0=theta0s,
    zeta0=np.zeros(num),
    pzeta0=np.full(num, pzeta),
    mu0=np.full(num, mu),
)

queue = dex.Queue(initial_conditions)

# =========================

queue.classify(machine)
queue.retain_energy(0, machine.quantity(60, "keV").to("NormJoule").m)

plane = dex.EnergyPzetaPlane(machine, mu)
pzeta_array = np.linspace(-1.3, 0.3, 2000) * machine.psip_last.value
plane.calculate_critical_points(pzeta_array)
plane.add_particles(queue)
fig, ax = plane.show(xlim=(-1.1, 0.3), ylim=(0.6, 2.8), curves=True, show=False)

queue.close(machine)
queue.plot_qkinetic_energy(machine)

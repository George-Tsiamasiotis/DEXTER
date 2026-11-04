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

mu = 3e-5

num = 1000
psit0_values = np.linspace(1e-6, LCFS.value, num)
psit0s = dex.MagneticFluxArray.Toroidal(psit0_values)
theta0s = np.zeros(num)

energy = machine.quantity(10.2, "keV").to("NormJoule").m
rho0s = dex.GCHamiltonian.eval_rho_array(machine, theta0s, psit0s, energy, mu)
rho0s[::2] *= -1  # alternate ρ signs to include counter-passing orbits


initial_conditions = dex.QueueInitialConditions.Boozer(
    t0=np.zeros(num),
    flux0=psit0s,
    theta0=theta0s,
    zeta0=np.zeros(num),
    rho0=rho0s,
    mu0=np.full(num, mu),
)

queue = dex.Queue(initial_conditions)

# =========================

queue.classify(machine)

plane = dex.EnergyPzetaPlane(machine, mu)
pzeta_array = np.linspace(-1.3, 0.3, 2000) * machine.psip_last.value
plane.calculate_critical_points(pzeta_array)
plane.add_particles(queue)
fig, ax = plane.show(xlim=(-1.1, 0.3), ylim=(0.6, 1.8), curves=True)

# Discard lost particles before integrating
plane.clear_particles()
queue.retain_orbit_types(
    [
        "Stagnated",
        "Potato",
        "TrappedConfined",
        "CoPassingConfined",
        "CuPassingConfined",
        "Undefined",
    ]
)
plane.add_particles(queue)
fig, ax = plane.show(xlim=(-1.1, 0.3), ylim=(0.6, 1.8), curves=True, show=False)

queue.close(machine)
dex.plot_qkinetic_pzeta(machine, queue)

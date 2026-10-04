import numpy as np
import dexter as dex

geometry = dex.LarGeometry(1, 1.75, 0.5)
LCFS = geometry.psi_last
machine = dex.Machine(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.5, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)


mu = 7e-5

num1 = 10
num2 = 7
num3 = 12
r0s = np.concat(
    (
        np.linspace(0.056, 0.205, num1),
        np.linspace(0.04, 0.15, num2),
        np.linspace(0.19, 0.4, num3),
    ),
)
theta0s = np.concat(
    (
        np.full(num1, 0),
        np.full(num2, np.pi),
        np.full(num3, np.pi),
    ),
)
psi0s = geometry.eval_psi_of_r(r0s)
num = len(psi0s)

initial_conditions = dex.QueueInitialConditions.Mixed(
    t0=np.zeros(num),
    flux0=dex.MagneticFluxArray.Toroidal(psi0s),
    theta0=theta0s,
    zeta0=np.zeros(num),
    pzeta0=np.full(num, -0.2 * machine.psip_last.value),
    mu0=np.full(num, mu),
)

queue = dex.Queue(initial_conditions)

queue.classify(machine)
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
plane = dex.EnergyPzetaPlane(machine, mu)
plane.add_particles(queue)
fig, ax = plane.show(xlim=(-0.4, 0.1), ylim=(0.8, 1.4), show=False)
fig.set_figheight(4)
fig.set_figwidth(3)
ax.legend().remove()


queue.close(machine, discard_arrays=False)
dex.plot_poloidal_drift(machine, queue)

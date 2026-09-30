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
            dex.FluteMode(2e-4, LCFS, 4, 3, np.pi),
            dex.FluteMode(4e-4, LCFS, 3, 1, np.pi),
        ]
    ),
)

num = 42
rs = np.linspace(0.01, geometry.rlast * 0.99, num)
num = len(rs)
psis = (rs / geometry.raxis) ** 2 / 2

initial_conditions = dex.QueueInitialConditions.Boozer(
    t0=np.zeros(num),
    flux0=dex.MagneticFluxArray.Toroidal(psis),
    theta0=np.zeros(num),
    zeta0=np.zeros(num),
    rho0=np.full(num, 1e-8),
    mu0=np.zeros(num),
)
queue = dex.Queue(initial_conditions)

intersect_params = dex.IntersectParams("ConstZeta", 0, 2000)
solver_params = dex.SolverParams(
    method=dex.SteppingMethod.EnergyAdaptiveStep(1e-7, 1e-8),
    max_steps=10_000_000,
)
queue.intersect(
    machine,
    intersect_params,
    solver_params,
)
print(f"Min energy variance: {queue.energy_rsd_array.min()}")
print(f"Max energy variance: {queue.energy_rsd_array.max()}")

# =========================


fig = plt.figure(layout="constrained", dpi=200, figsize=(3.2, 3))
ax = fig.subplots()

rlab_last = geometry.rlab_last
zlab_last = geometry.zlab_last

for p in queue.particles():
    if p.steps_stored == 0:
        continue
    r = geometry.eval_rlab(psi=p.psi_array, theta=p.theta_array)
    z = geometry.eval_zlab(psi=p.psi_array, theta=p.theta_array)
    ax.plot(
        r,
        z,
        c="b",
        marker=".",
        markersize=1.2,
        markeredgewidth=0,
        alpha=0.7,
        linestyle="",
    )


ax.plot(rlab_last, zlab_last, c="k", linewidth=2, zorder=2)
ax.scatter(1.75, 0, marker="+", s=20, c="k")

ax.margins(0.005, 0.005)
ax.set_aspect("equal")
ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.set_xticks((1.25, 1.5, 1.75, 2, 2.25))
ax.set_yticks((-0.5, -0.25, 0, 0.25, 0.5))
plt.show()

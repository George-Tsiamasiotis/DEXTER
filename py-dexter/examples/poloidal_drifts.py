"""Classification of all particles in a LAR machine."""

import numpy as np
import dexter as dex
import matplotlib.pyplot as plt
from dexter.simulate.plot import orbit_color

geometry = dex.LarGeometry(1, 1.75, 0.5)
LCFS = geometry.psi_last
machine = dex.Machine(
    geometry=geometry,
    qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.5, lcfs=LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

mu = 6e-5

points = (
    (-0.1, 0.026, 1),  # CoPassing-Confined
    (-0.8, 0.001, 1),  # CuPassing-Confined
    (-0.5, 0.018, 1),  # Trapped-Confined
    (-0.8, 0.02, 1),  # Trapped-Lost
    (-0.0448, 0.0045, 1),  # Potato
    (-0.0, 0.0014, 1),  # Stagnated
    (-0.3, 0.029, 2),  # CoPassing-Lost
    (-1.2, 0.027, 0),  # CuPassing-Lost
)

points_array = np.asarray(points).T
pzeta0 = points_array[0] * machine.psip_last.value
psi0 = points_array[1]
theta0 = points_array[2]
particle_count = len(pzeta0)

initial_conditions = dex.QueueInitialConditions.Mixed(
    t0=np.zeros(particle_count),
    flux0=dex.MagneticFluxArray.Toroidal(psi0),
    theta0=theta0,
    zeta0=np.zeros(particle_count),
    pzeta0=pzeta0,
    mu0=np.full(particle_count, mu),
)
queue = dex.Queue(initial_conditions)

queue.classify(machine)
queue.close(machine, discard_arrays=False)


# ======================

fig, ax = dex.plot_poloidal_drift(machine, queue, show=False)
fig.set_figheight(2.3)
fig.set_figwidth(4.5)

for i, p in enumerate(queue.particles()):
    color = orbit_color(p.orbit_type)
    label = rf"${p.orbit_type}$"
    ax.plot([], [], c=color, linewidth=1.5, label=label, zorder=-4)

ax.scatter(1.75, 0, c="k", marker="+", s=40, zorder=3, label=r"$R_{axis}$")
ax.set_xticks((1.25, 1.5, 1.75, 2, 2.25))
ax.set_yticks((-0.5, -0.25, 0, 0.25, 0.5))
ax.legend(bbox_to_anchor=(1.1, 1.05))

plt.show()

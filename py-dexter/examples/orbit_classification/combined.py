"""Classification of all particles in a LAR machine.

This script is a combination of all previous tests.
"""

import numpy as np
import dexter as dex
from math import sqrt, pi
import matplotlib.pyplot as plt
from matplotlib import patheffects

from dexter.simulate.plot import orbit_color

LCFS = dex.MagneticFlux.Toroidal(0.03)
raxis = 1.75
rlast = sqrt(2 * LCFS.value) * raxis  # `rlast` must be in [m]
machine = dex.Machine(
    geometry=dex.LarGeometry(baxis=1, raxis=raxis, rlast=rlast),
    qfactor=dex.ParabolicQfactor(1.1, 3.9, LCFS),
    current=dex.LarCurrent(),
    bfield=dex.LarBfield(),
)

mu = 6e-5

points = (
    # Pζ, ψ, θ
    (-0.8, 0.001, 0),  # Alpha
    (-1.5, 0.02, 1),  # Beta
    (-0.8, 0.01, 1),  # Gamma
    (-0.6, 0.018, pi),  # Delta
    (-0.6, 0.003, pi),  # Epsilon
    (-0.4, 0.025, pi),  # Zeta
    (-0.1, 0.015, 1),  # Eta
    (-0.0448, 0.0045, 1),  # Theta
    (-0.6, 0.016, 1),  # Iota
    (-0.36, 0.025, 0),  # Kappa
    (-0.36, 0.001, 0),  # Lambda
    (-0.4, 0.01, 0),  # Mu - Trapped
    (0.1, 0.0013, 0),  # Mu - CoPassing
    (-0.3, 0.0005, pi),  # Mu - CuPassing
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

plane = dex.EnergyPzetaPlane(machine, mu)
plane.add_particles(queue)
fig, ax = plane.show(particles=True, ylim=(0.4, 2.5), show=False)
fig.set_figwidth(9)

for p in queue.particles():
    xy = (
        p.initial_conditions.pzeta0 / machine.psip_last.value,
        p.initial_energy / mu,
    )
    color = orbit_color(p.orbit_type)
    ax.annotate(
        rf"${p.energy_pzeta_position}-{p.orbit_type}$",
        xy=xy,
        xytext=(xy[0], xy[1] + 0.2),
        zorder=100,
        fontsize=8,
        arrowprops=dict(arrowstyle="->", connectionstyle="angle3", lw=2, color=color),
        path_effects=[patheffects.withStroke(linewidth=3, foreground="w")],
    )

plt.show()

"""Classification of all particles in a LAR equilibrium.

This script is a combination of all previous tests.
"""

import numpy as np
import dexter as dex
from math import sqrt, pi
import matplotlib.pyplot as plt
from matplotlib import patheffects

from dexter.simulate.colors import orbit_color

LCFS = dex.LastClosedFluxSurface("Toroidal", 0.03)
raxis = 1.75
rlast = sqrt(2 * LCFS.value) * raxis  # `rlast` must be in [m]
equilibrium = dex.Equilibrium(
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
    (-0.0, 0.0014, 1),  # Mu
)

points_array = np.asarray(points).T
pzeta0 = points_array[0] * equilibrium.psip_last
psi0 = points_array[1]
theta0 = points_array[2]
particle_count = len(pzeta0)
initial_conditions = dex.QueueInitialConditions.mixed(
    t0=np.zeros(particle_count),
    flux0=dex.InitialFluxArray("Toroidal", psi0),
    theta0=theta0,
    zeta0=np.zeros(particle_count),
    pzeta0=pzeta0,
    mu0=np.full(particle_count, mu),
)

queue = dex.Queue(initial_conditions)
queue.classify(equilibrium)

fig, ax = dex.plot_parabolas(equilibrium, mu, show=False)
for p in queue.particles:
    xy = (
        p.initial_conditions.pzeta0 / equilibrium.psip_last,
        p.initial_energy / mu,
    )
    color = orbit_color(p.orbit_type)
    ax.scatter(*xy, s=20, c=color, zorder=5)
    ax.annotate(
        rf"${p.energy_pzeta_position}-{p.orbit_type}$",
        xy=xy,
        xytext=(xy[0], xy[1] + 0.2),
        zorder=10,
        arrowprops=dict(arrowstyle="->", connectionstyle="angle3", lw=2, color=color),
        path_effects=[patheffects.withStroke(linewidth=3, foreground="w")],
    )

plt.show()

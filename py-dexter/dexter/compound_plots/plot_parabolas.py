"""Plots the orbit classification parabolas on the (E, Pζ, μ=const) space."""

from collections import Counter
from math import isfinite
import numpy as np
import matplotlib.pyplot as plt

from dexter import EnergyPzetaPlane, Equilibrium, COMs, OrbitType, Particle, Canvas
from dexter.simulate._plotters import orbit_color, orbit_color_legend_handles

PARABOLAS_FIG_KW = {"figsize": (9, 6), "layout": "constrained", "dpi": 120}
MAGNETIC_AXIS_KW = {
    "color": "xkcd:cobalt blue",
    "label": r"$Magnetic\ Axis$",
    "linewidth": 2,
    "linestyle": "--",
}
LEFT_WALL_KW = {
    "color": "xkcd:coral",
    "label": r"$Left\ Wall$",
    "linewidth": 2,
    "linestyle": "--",
}
RIGHT_WALL_KW = {
    "color": "xkcd:forrest green",
    "label": r"$Right\ Wall$",
    "linewidth": 2,
    "linestyle": "--",
}
TP_BOUNDARY_KW = {
    "color": "xkcd:bright pink",
    "linewidth": 2,
    "linestyle": "--",
}


def plot_parabolas(
    equilibrium: Equilibrium,
    mu: float,
    particles: list[Particle] | None = None,
    xlim: tuple[float, float] = (-1.6, 0.5),
    ymax: float = 3,
    density: int = 1000,
    show: bool = True,
) -> Canvas:
    r"""Plots the orbit classification parabolas on the $(E, P_\zeta, \mu=const)$ space.

    Parameters
    ----------
    equilibrium
        The equilibrium in which to construct the COM space.
    mu
        The magnetic moment $\mu$.

    Other Parameters
    ----------------
    particles
        A list of particles to project on the E-Pζ plane.
    xlim
        The xaxis limits, normalized to $\psi_{p,wall}$. Defaults to (-1.6, 0.5).
    ymax
        The yaxis upper limit. If not provided, the axes are autoscaled by the
        parabolas. Note that by definition, the minimum of the magnetic axis parabola is
        always the point $(0, 1)$. Defaults to 3.
    density
        The number of points to evaluate each parabola on. Defaults to 1000.
    show
        Whether or not to call `plt.show()`. Defaults to True.

    Returns
    -------
    Canvas
        The produced `Figure` and `Ax`.
    """

    fig = plt.figure(**PARABOLAS_FIG_KW)
    ax = fig.add_subplot()
    fig.suptitle(
        rf"Orbit classification parabolas on the $(E, P_\zeta, \mu={mu})$ space."
    )

    psip_last = equilibrium.psip_last
    xaxis_span = np.linspace(xlim[0], xlim[1], density)
    pzeta_span = xaxis_span * psip_last
    energy_norm = mu

    coms = COMs(mu=mu)
    plane: EnergyPzetaPlane = coms.build_energy_pzeta_plane(equilibrium)

    axis = plane.axis_parabola.eval_array1(pzeta_span) / energy_norm
    left_wall = plane.left_wall_parabola.eval_array1(pzeta_span) / energy_norm
    right_wall = plane.right_wall_parabola.eval_array1(pzeta_span) / energy_norm
    tp_upper = plane.tp_upper / energy_norm
    tp_lower = plane.tp_lower / energy_norm
    tpx = plane.tp_pzeta_interval / psip_last  # not always a linspace

    ax.plot(xaxis_span, axis, **MAGNETIC_AXIS_KW)
    ax.plot(xaxis_span, left_wall, **LEFT_WALL_KW)
    ax.plot(xaxis_span, right_wall, **RIGHT_WALL_KW)

    ax.plot(tpx, tp_upper, **TP_BOUNDARY_KW)
    ax.plot(tpx, tp_lower, **TP_BOUNDARY_KW)
    ax.plot(
        [tpx[0]] * 2,
        [tp_lower[0], tp_upper[0]],
        **(TP_BOUNDARY_KW | dict(label=r"$Trapped-Passing\ Boundary$")),
    )  # close

    if particles is not None:
        pzetas = []
        energies = []
        colors = []
        for particle in particles:
            try:
                pzeta = particle.initial_conditions.pzeta0
                energy = particle.initial_energy
                if isfinite(pzeta) and isfinite(energy):
                    pzetas.append(pzeta / equilibrium.psip_last)
                    energies.append(particle.initial_energy / energy_norm)
                    colors.append(orbit_color(particle.orbit_type))
            except AttributeError:
                continue
        ax.scatter(pzetas, energies, c=colors, s=1)
        found_orbit_types = Counter([particle.orbit_type for particle in particles])
        ax.legend(handles=orbit_color_legend_handles(found_orbit_types))
    else:
        ax.legend()

    ax.margins(0)
    ax.set_ylim(0, ymax)
    ax.set_xlabel(r"$P_\zeta/\psi_{p,wall}$")
    ax.set_ylabel(r"$E/\mu$")
    ax.grid(True)

    if show:
        plt.show()
        plt.close()

    return (fig, ax)

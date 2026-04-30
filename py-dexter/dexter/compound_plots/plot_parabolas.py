"""Plots the orbit classification parabolas on the (E, Pζ, μ=const) space."""

from collections import Counter
from math import isfinite
from typing import Literal
from matplotlib.tri import Triangulation
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from alpha_shapes import Alpha_Shaper

from dexter import EnergyPzetaPlane, Equilibrium, COMs, Particle, Queue
from dexter.simulate._plotters import orbit_color, orbit_color_legend_handles
from dexter.types import Array1, Canvas

PARABOLAS_FIG_KW = {"figsize": (9, 6), "layout": "constrained", "dpi": 120}
MAGNETIC_AXIS_KW = {
    "color": "xkcd:cobalt blue",
    "label": r"$Magnetic\ Axis$",
    "linewidth": 3,
    "linestyle": "-",
}
LEFT_WALL_KW = {
    "color": "xkcd:coral",
    "label": r"$Left\ Wall$",
    "linewidth": 3,
    "linestyle": "-",
}
RIGHT_WALL_KW = {
    "color": "xkcd:forrest green",
    "label": r"$Right\ Wall$",
    "linewidth": 3,
    "linestyle": "-",
}
TP_BOUNDARY_KW = {
    "color": "xkcd:bright pink",
    "linewidth": 3,
    "linestyle": "-",
    "zorder": 5,
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
    )

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


def plot_qkinetic_tricontour(
    equilibrium: Equilibrium,
    queue: Queue,
    levels: Array1 | list[float] | None = None,
    clabel: bool = True,
    qmax: float = np.inf,
    alpha: float = 30,
    xlim: tuple[float, float] = (-1.6, 0.5),
    ymax: float = 3,
    density: int = 1000,
    show: bool = True,
) -> Canvas:
    r"""Plots a tricontour of the $q_{kinetic}$ of a scatter grid of particles on the $E-P_\zeta$ plane.

    Note that the particles must have the same $\mu$ for this plot to make sense.

    Parameters
    ----------
    equilibrium
        The Equilibrium in which the particles where integrated.
    queue
        The Queue with the already integrated particles.

    Other Parameters
    ----------------
    levels
        The levels of the contour lines. If an array is passed then its values are plotter. If None, the
        levels are $[q_{min}, ..., q_{max}]$, rounded to the closest half-integer. Defaults to None.
    clabel
        Whether or not to add a label on the contour lines. Defaults
        to True.
    qmax
        An upper limit to the $q_kinetic$ values. Useful when particle close to the separatrix appear.
        Defaults to np.inf.
    alpha
        The $\alpha$ parameter passed to `alpha_shapes`. Defaults to 30.
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

    """

    mu = queue[0].initial_conditions.mu0
    if not np.all(queue.initial_conditions.mu_array == mu):
        raise ValueError("All particles must have the same 'μ'")

    particles = queue.particles
    fig, ax = plot_parabolas(
        equilibrium,
        mu,
        xlim=xlim,
        ymax=ymax,
        density=density,
        show=False,
    )
    fig.suptitle(r"$q_{qkinetic}$ on the " rf"$(E, P_\zeta, \mu={mu})$ space.")
    ax.grid(False)

    energy_norm = mu

    trapped = [
        particle
        for particle in particles
        if particle.energy_pzeta_position in ["Iota", "Kappa", "Nu"]
    ]
    copassing = [
        particle
        for particle in particles
        if particle.energy_pzeta_position in ["Theta", "Lambda"]
    ]
    cupassing = [
        particle
        for particle in particles
        if particle.energy_pzeta_position in ["Alpha", "Zeta", "Epsilon", "Mu"]
    ]

    if levels is None:
        qkinetic_array = queue.qkinetic_array
        qmin = np.nanmin(qkinetic_array)
        qmax = np.nanmax(qkinetic_array)
    else:
        qmin = min(levels)
        qmax = max(levels)
    levels = np.arange(np.floor(qmin), np.ceil(qmax), 0.5)

    for family in [trapped, copassing, cupassing]:
        if len(family) == 0:
            continue
        pzetas = []
        energies = []
        qkinetics = []
        colors = []
        for particle in family:
            try:
                pzeta = particle.initial_conditions.pzeta0
                energy = particle.initial_energy
                qkinetic = particle.qkinetic
            except AttributeError:
                continue
            if qkinetic > qmax:
                continue
            if not isfinite(qkinetic):
                continue
            pzetas.append(pzeta / equilibrium.psip_last)
            energies.append(energy / energy_norm)
            qkinetics.append(qkinetic)
            colors.append(orbit_color(particle.orbit_type))

        if len(qkinetics) < 3:  # Problems with Triangulation
            continue

        points = np.asarray((pzetas, energies))
        if len(levels) > 4:  # requirement for alpha_shapes
            triang = Alpha_Shaper(points.T, normalize=False)
            triang.set_mask_at_alpha(alpha)
        else:
            triang = Triangulation(points[0], points[1])

        tri = ax.tricontour(
            triang,
            qkinetics,
            levels=levels,
            vmin=qmin,
            vmax=qmax,
            cmap="cool",
        )

        if clabel:
            ax.clabel(
                tri,
                levels=levels,
                fmt=contour_clabel_fmt,
                fontsize=10,
                manual=False,
            )

    if not clabel:
        cmap = plt.get_cmap("cool")
        norm = BoundaryNorm(levels, cmap.N)
        cm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
        fig.colorbar(cm, ticks=levels, ax=ax, label="$q_{kinetic}$")

    ax.margins(0)
    ax.set_ylim(0, ymax)
    ax.set_xlabel(r"$P_\zeta/\psi_{p,wall}$")
    ax.set_ylabel(r"$E/\mu$")

    if show:
        plt.show()
        plt.close()

    return (fig, ax)


# ================================================================================================


def contour_clabel_fmt(x: float) -> str:
    """Formats the contour line labels."""
    ratio = x.as_integer_ratio()
    numerator = ratio[0]
    denominator = ratio[1]

    _numerator_str = f"{numerator:.1f}"
    _denominator_str = f"{denominator:.1f}"
    if _numerator_str.endswith("0"):
        _numerator_str = f"{numerator:.0f}"
    if _denominator_str.endswith("0"):
        _denominator_str = f"{denominator:.0f}"

    return f"{_numerator_str}/{_denominator_str}"

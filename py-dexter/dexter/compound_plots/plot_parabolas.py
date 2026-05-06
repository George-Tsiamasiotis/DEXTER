"""Plots the orbit classification parabolas on the (E, Pζ, μ=const) space."""

from collections import Counter
from fractions import Fraction
from math import isfinite
from matplotlib.ticker import Formatter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from alpha_shapes import Alpha_Shaper

from dexter import EnergyPzetaPlane, Equilibrium, COMs, Particle, Queue
from dexter.simulate.colors import orbit_color, _orbit_color_legend_handles
from dexter.types import Array1, Canvas, EnergyPzetaPosition

PARABOLAS_FIG_KW = {"figsize": (7, 5), "layout": "constrained", "dpi": 160}
MAGNETIC_AXIS_KW = {
    "color": "xkcd:cobalt blue",
    "label": r"$Magnetic\ Axis$",
    "linewidth": 3,
    "linestyle": "-",
    "zorder": 5,
}
LEFT_WALL_KW = {
    "color": "xkcd:coral",
    "label": r"$Left\ Wall$",
    "linewidth": 3,
    "linestyle": "-",
    "zorder": 5,
}
RIGHT_WALL_KW = {
    "color": "xkcd:forrest green",
    "label": r"$Right\ Wall$",
    "linewidth": 3,
    "linestyle": "-",
    "zorder": 5,
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
        **(TP_BOUNDARY_KW | dict(label=r"$Trapped-Passing\ Boundary$")),  # type: ignore
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
        ax.legend(handles=_orbit_color_legend_handles(found_orbit_types))
    else:
        ax.legend()

    if equilibrium._has_pint:
        si_ax = ax.twinx()
        max_energy_nu = ymax * mu
        max_energy_si = equilibrium.quantity(max_energy_nu, "energy_units").to("keV")
        si_ax.set_ybound(0, max_energy_si.m)
        si_ax.set_ylabel(r"$Energy\ [keV]$")

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
    levels: Array1 | list[float] | int = 20,
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
        The levels of the contour lines. If an array is passed then its values are plotted. If an int $n$ is passed,
        the levels are a linear space from $q_{min}$ to $q_{max}$ of length $n$. Defaults to 20.
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

    all_particles = queue.particles
    fig, ax = plot_parabolas(
        equilibrium,
        mu,
        xlim=xlim,
        ymax=ymax,
        density=density,
        show=False,
    )
    fig.suptitle(r"$q_{kinetic}$ on the " rf"$(E, P_\zeta, \mu={mu})$ space.")
    ax.grid(False)

    energy_norm = mu

    all_qkinetics = queue.qkinetic_array
    qmin = np.nanmin(all_qkinetics)
    qmax = np.nanmax(all_qkinetics)

    if isinstance(levels, int):
        levels = np.linspace(qmin, qmax, levels)
        formatter = None
    else:
        levels = np.sort(levels)
        formatter = FractionFormatter()

    families: tuple[EnergyPzetaPosition] = EnergyPzetaPosition.__args__

    for family in families:
        pzetas = []
        energies = []
        qkinetics = []
        particles = [p for p in all_particles if p.energy_pzeta_position == family]
        if len(particles) == 0:
            continue
        for particle in particles:
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

        if len(qkinetics) < 3:  # Problems with Triangulation
            continue

        points = np.asarray((pzetas, energies)).T

        triang = Alpha_Shaper(points, normalize=True)
        triang.set_mask_at_alpha(alpha)

        ax.tricontour(
            triang,
            qkinetics,
            levels=levels,
            vmin=qmin,
            vmax=qmax,
            cmap="cool",
        )

    cmap = plt.get_cmap("cool")
    norm = BoundaryNorm(levels, cmap.N)
    cm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    fig.colorbar(
        cm,
        ticks=levels,
        ax=ax,
        format=formatter,
        drawedges=False,
        label="$q_{kinetic}$",
    )

    ax.margins(0)
    ax.set_ylim(0, ymax)
    ax.set_xlabel(r"$P_\zeta/\psi_{p,wall}$")
    ax.set_ylabel(r"$E/\mu$")

    if show:
        plt.show()
        plt.close()

    return (fig, ax)


# ================================================================================================


class FractionFormatter(Formatter):
    r"""Formats values as integer ratios."""

    denominator_limit: int = 25

    def format_data(self, value: float) -> str:
        return str(Fraction(value).limit_denominator(self.denominator_limit))

    def __call__(self, x: float, _: int | None = None) -> str:
        return self.format_data(x)

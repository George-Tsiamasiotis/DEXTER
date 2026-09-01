"""Plots the orbit classification parabolas on the (E, Pζ, μ=const) space."""

from collections import Counter
from fractions import Fraction
from math import isfinite
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.contour import ContourSet
from matplotlib.ticker import Formatter
from alpha_shapes import Alpha_Shaper

from dexter import EnergyPzetaPlane, Equilibrium, COMs, Particle, Queue
from dexter.simulate.colors import orbit_color, _orbit_color_legend_handles
from dexter.types import Array1, Canvas, EnergyPzetaPosition

PSIP_REV = 0.09157457

PARABOLAS_FIG_KW = {"figsize": (7, 5), "layout": "constrained", "dpi": 160}
PARABOLAS_DENSITY: int = 2000
MAGNETIC_AXIS_KW = {
    "color": "xkcd:cobalt blue",
    # "color": "gray",
    "alpha": 1,
    "label": r"$Magnetic\ Axis$",
    "linewidth": 2,
    "linestyle": "-",
    "zorder": 5,
}
LEFT_WALL_KW = {
    "color": "xkcd:coral",
    # "color": "gray",
    "alpha": 1,
    "label": r"$Left\ Wall$",
    "linewidth": 2,
    "linestyle": "-",
    "zorder": 5,
}
RIGHT_WALL_KW = {
    "color": "xkcd:forrest green",
    # "color": "gray",
    "alpha": 1,
    "label": r"$Right\ Wall$",
    "linewidth": 2,
    "linestyle": "-",
    "zorder": 5,
}
TP_BOUNDARY_KW = {
    "color": "xkcd:bright pink",
    # "color": "gray",
    "alpha": 1,
    "linewidth": 2,
    "linestyle": "-",
    "zorder": 5,
}


def plot_parabolas(
    equilibrium: Equilibrium,
    mu: float,
    particles: list[Particle] | None = None,
    ymax: float = 3,
    emin: bool = True,
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
    ymax
        The yaxis upper limit. If not provided, the axes are autoscaled by the
        parabolas. Note that by definition, the minimum of the magnetic axis parabola is
        always the point $(0, 1)$. Defaults to 3.
    show
        Whether or not to call `plt.show()`. Defaults to True.

    Returns
    -------
    Canvas
        The produced `Figure` and `Ax`.
    """

    fig = plt.figure(**PARABOLAS_FIG_KW)
    ax = fig.add_subplot()
    # fig.suptitle(
    #     rf"Orbit classification parabolas on the $(E, P_\zeta, \mu={mu})$ space."
    # )

    coms = COMs(mu=mu)
    energy_norm = mu
    psip_last = equilibrium.psip_last
    plane: EnergyPzetaPlane = coms.build_energy_pzeta_plane(equilibrium)

    # Creates x-spans from the intersept points
    # Note that parabolas are defined as `E(Pζ)`
    intercept = ymax * mu
    ax_ips = plane.axis_parabola.horizontal_line_intercepts(intercept)
    lw_ips = plane.left_wall_parabola.horizontal_line_intercepts(intercept)
    rw_ips = plane.right_wall_parabola.horizontal_line_intercepts(intercept)
    ax_xaxis_span = np.linspace(ax_ips[0][0], ax_ips[1][0], PARABOLAS_DENSITY)
    lw_xaxis_span = np.linspace(lw_ips[0][0], lw_ips[1][0], PARABOLAS_DENSITY)
    rw_xaxis_span = np.linspace(rw_ips[0][0], rw_ips[1][0], PARABOLAS_DENSITY)
    tp_upper = plane.tp_upper / energy_norm
    tp_lower = plane.tp_lower / energy_norm

    axis = plane.axis_parabola.eval_array1(ax_xaxis_span) / energy_norm
    left_wall = plane.left_wall_parabola.eval_array1(lw_xaxis_span) / energy_norm
    right_wall = plane.right_wall_parabola.eval_array1(rw_xaxis_span) / energy_norm
    tpx = plane.tp_pzeta_interval / psip_last  # not always a linspace

    ax.plot(ax_xaxis_span / psip_last, axis, **MAGNETIC_AXIS_KW)
    ax.plot(lw_xaxis_span / psip_last, left_wall, **LEFT_WALL_KW)
    ax.plot(rw_xaxis_span / psip_last, right_wall, **RIGHT_WALL_KW)

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
        ax.scatter(
            pzetas, energies, c=colors, s=3, zorder=2, edgecolor=None, linewidths=0
        )
        found_orbit_types = Counter([particle.orbit_type for particle in particles])
        handles = _orbit_color_legend_handles(found_orbit_types)
        ax.legend(
            handles=handles,
            bbox_to_anchor=(1.2, 0.44),
        )
    else:
        ax.legend(loc="lower right")

    if emin:
        pzeta_span = (rw_ips[0][0], ax_ips[1][0])
        pzetas_min, emins, emaxs = calculate_energy_min(equilibrium, mu, pzeta_span)
        ax.plot(
            pzetas_min / equilibrium.psip_last,
            emins / mu,
            c="cyan",
            linewidth=2,
            linestyle="--",
            zorder=10,
            label=r"$E_{min}$",
        )
        ax.plot(
            pzetas_min / equilibrium.psip_last,
            emaxs / mu,
            c="xkcd:light green",
            linewidth=2,
            linestyle="--",
            zorder=10,
            label=r"$E_{max}$",
        )

    if equilibrium._has_pint:
        si_ax = ax.twinx()
        max_energy_nu = ymax * mu
        max_energy_si = equilibrium.quantity(max_energy_nu, "energy_units").to("keV")
        si_ax.set_ybound(0, max_energy_si.m)
        si_ax.set_ylabel(r"$Energy\ [keV]$")

    ax.set_ylim(0, ymax)
    ax.set_xlabel(r"$P_\zeta/\psi_{p,LCFS}$")
    ax.set_ylabel(r"$E/\mu$")
    ax.grid(True)
    # ax.legend()

    if show:
        plt.show()
        plt.close()

    return (fig, ax)


def calculate_energy_min(
    equilibrium: Equilibrium,
    mu: float,
    pzeta_span: tuple[float, float],
) -> tuple[Array1, Array1, Array1]:
    theta = np.concat(
        (
            np.full(1000, 0),
            np.full(1000, np.pi),
        )
    )

    if equilibrium.geometry.psi_state == "Good":
        psi = np.concat(
            (
                np.linspace(0, equilibrium.psi_last, 1000),
                np.linspace(0, equilibrium.psi_last, 1000),
            )
        )
        psip = equilibrium.qfactor.psip_of_psi(psi)
        b = equilibrium.bfield.b_of_psi(psi, theta)
        g = equilibrium.current.g_of_psi(psi)
    else:
        psip = np.concat(
            (
                np.linspace(0, equilibrium.psip_last, 1000),
                np.linspace(0, equilibrium.psip_last, 1000),
            )
        )
        b = equilibrium.bfield.b_of_psip(psip, theta)
        g = equilibrium.current.g_of_psip(psip)

    def energy(pzeta: float) -> np.ndarray:
        return ((pzeta + psip) * b) ** 2 / (2 * g**2) + mu * b

    emins = []
    emaxs = []
    pzeta = np.linspace(*pzeta_span, 1000)
    for i in range(0, len(pzeta)):
        energy_array = energy(pzeta[i])
        emin = energy_array.min()
        emax = energy_array.max()
        emins.append(emin)
        emaxs.append(emax)

    return np.asarray(pzeta), np.asarray(emins), np.asarray(emaxs)


def plot_qkinetic_tricontour(
    equilibrium: Equilibrium,
    queue: Queue,
    levels: Array1 | list[float] | int = 20,
    qmax: float = np.inf,
    alpha: float = 30,
    ymax: float = 3,
    emin: bool = True,
    clabel: bool = True,
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
    ymax
        The yaxis upper limit. If not provided, the axes are autoscaled by the
        parabolas. Note that by definition, the minimum of the magnetic axis parabola is
        always the point $(0, 1)$. Defaults to 3.
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
        ymax=ymax,
        emin=emin,
        show=False,
    )
    scale = 0.6
    fig.set_figheight(5 * scale)
    fig.set_figwidth(8 * scale)
    # fig.suptitle(
    #     r"$q_{kinetic}\ on\ the\ $" + rf"$(E, P_\zeta, \mu={mu})\ $" + "$space.$"
    # )
    ax.grid(False)

    # pitch = equilibrium.quantity(112, "keV").to("energy_unit").m / mu
    # print(pitch)
    # ax.set_xbound(-1.05, 0.05)
    # ax.axhline(
    #     y=pitch, c="k", alpha=0.7, linestyle="--", label=r"$E=112keV$", zorder=-2
    # )

    energy_norm = mu

    all_qkinetics = queue.qkinetic_array
    vmin = max(np.nanmin(all_qkinetics), -qmax)
    vmax = min(np.nanmax(all_qkinetics), qmax)

    if isinstance(levels, int):
        levels = np.sort(np.linspace(vmin, vmax, levels))
        formatter = None
    else:
        # sort, remove duplicates and flatten
        vmin = max(np.nanmin(all_qkinetics), -qmax, np.min(levels))
        vmax = min(np.nanmax(all_qkinetics), qmax, np.max(levels))
        levels = np.sort(list(set(levels)))
        formatter = FractionFormatter()

    # family_groups: list[list[EnergyPzetaPosition]] = [
    #     ["Kappa", "Eta"],
    #     ["Theta", "Iota", "Mu"],
    #     ["Alpha", "Epsilon", "Lambda"],
    # ]
    family_groups: list[list[EnergyPzetaPosition]] = [
        [
            "Lambda",
        ],
        [
            "Iota",
            "Theta",
            "Mu",
        ],
        [
            "Kappa",
            "Eta",
        ],
    ]
    family_groups = [EnergyPzetaPosition.__args__]
    allsegs = []
    allkinds = []
    dup_levels = []
    for families in family_groups:
        pzetas = []
        energies = []
        qkinetics = []
        particles = [p for p in all_particles if p.energy_pzeta_position in families]
        if len(particles) == 0:
            continue
        for particle in particles:
            try:
                pzeta = particle.initial_conditions.pzeta0
                energy = particle.initial_energy
                qkinetic = particle.qkinetic
            except AttributeError:
                continue
            # if abs(qkinetic) > qmax:
            #     continue
            # if abs(qkinetic) < 0.05 or abs(qkinetic) > 1.4:
            #     continue
            if not isfinite(qkinetic):
                continue
            pzetas.append(pzeta)
            energies.append(energy)
            qkinetics.append(qkinetic)

        if len(qkinetics) < 20:  # Problems with Triangulation
            continue

        pzetas = np.asarray(pzetas) / equilibrium.psip_last
        energies = np.asarray(energies) / energy_norm
        points = np.asarray((pzetas, energies)).T

        triang = Alpha_Shaper(points, normalize=True)
        triang.set_mask_at_alpha(alpha)

        tri = ax.tricontour(
            triang,
            qkinetics,
            levels=levels,
            vmin=vmin,
            vmax=vmax,
            cmap="cool",
            linewidths=1,
        )
        allsegs += tri.allsegs
        allkinds += tri.allkinds
        dup_levels += list(levels)

    # =========================

    # =========================

    if clabel:
        # Need to combine previous tricontours to use manual clabels
        cs = ContourSet(
            ax, dup_levels, allsegs, linewidths=0, cmap="cool", vmin=vmin, vmax=vmax
        )
        ax.clabel(
            cs,
            fmt=FractionFormatter(),
            fontsize=8,
            zorder=10,
            manual=True,
        )
    # else:
    # cmap = plt.get_cmap("cool")
    # norm = BoundaryNorm(levels, cmap.N)
    # cm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    # fig.colorbar(
    #     cm,
    #     ticks=levels,
    #     ax=ax,
    #     format=formatter,
    #     drawedges=False,
    #     label="$q_{kinetic}$",
    # )

    ax.margins(0)
    ax.set_ylim(0, ymax)
    ax.set_xlabel(r"$P_\zeta/\psi_{p,wall}$")
    ax.set_ylabel(r"$E/\mu$")

    if show:
        raise Exception
        plt.show()
        plt.close()

    return (fig, ax)


# ================================================================================================


class FractionFormatter(Formatter):
    r"""Formats values as integer ratios."""

    denominator_limit: int = 250

    def format_data(self, value: float) -> str:
        return str(Fraction(value).limit_denominator(self.denominator_limit))

    def format_ticks(self, values: list[float]) -> list[str]:
        return [self.format_data(value) for value in values]

    def __call__(self, x: float, _: int | None = None) -> str:
        return self.format_data(x)

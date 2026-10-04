r"""Plotting functions for simulated objects.

Functions
---------
plot_evolution
    Plots the time evolution of a particle's dynamical variables.
plot_poloidal_drift
    Plots a particle's drift on the $R-Z$ plane, overlaid on a contour plot of the Hamiltonian.
"""

import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.ticker import LogLocator, MaxNLocator
from cycler import cycler
from math import floor, log10
from typing import cast
from pint.facets.plain import PlainQuantity

from dexter.machine.machine import Machine
from dexter.machine.geometries import GeometryObject
from dexter.simulate.params import IntersectParams
from dexter.simulate.particle import Particle
from dexter.simulate.queue import Queue
from dexter.simulate.energy import (
    create_poloidal_grid,
    energy_of_psi_grid,
    energy_of_psip_grid,
)
from dexter._utils import tex_unit
from dexter.machine.plot import _resolve_magnetic_flux_kind
from dexter.types import Array2, ArrayShape, Locator, Array1, OrbitType

TAU = 2 * np.pi
PI = np.pi


def plot_evolution(
    machine: Machine,
    particle: Particle,
    downsample: bool = True,
    show: bool = True,
) -> tuple[Figure, tuple[Axes]]:
    r"""Plots the time evolution of a particle's dynamical variables.

    Parameters
    ----------
    machine
        The machine in which the particle was integrated. It is used to convert specific
        quantities to SI units.
    particle
        The integrated particle.
    downsample
        Whether or not to downsample the evolution arrays. This can be
        really helpful when plotting arrays with a lot of points, since it
        drastically improves both figure creation and interaction
        performance. Downsampling is done by increasing the [::step] just
        enough so that the final number of points is more than 50.000.
    show
        Whether or not to call `plt.show()`.

    Raises
    ------
    RuntimeError
        If the particle has not been integrated (`#!python particle.steps_taken == 0`).

    """
    if particle.steps_taken == 0:
        raise RuntimeError("Particle has not been integrated")

    fig = plt.figure(figsize=(9, 5), dpi=140)
    axes = fig.subplots(4, 2, sharex=True)
    axpsi = axes[0, 0]
    axtheta = axes[1, 0]
    axptheta = axes[2, 0]
    axrho = axes[3, 0]
    axpsip = axes[0, 1]
    axzeta = axes[1, 1]
    axpzeta = axes[2, 1]
    axenergy = axes[3, 1]

    DOWNSAMPLE_TARGET = 50_000

    t_array = particle.t_array

    step = 1
    if downsample:
        length = len(t_array)
        oom = floor(log10(length))  # order of magnitude of number of points
        target_oom = floor(log10(DOWNSAMPLE_TARGET))
        if oom > target_oom:
            step = 10 ** (oom - target_oom)

    if downsample and not len(t_array) < DOWNSAMPLE_TARGET * 10:
        warnings.warn("Downsampling did not work..")

    t = t_array[::step]
    psi = particle.psi_array[::step] / machine.psi_last.value
    psip = particle.psip_array[::step] / machine.psip_last.value
    theta = particle.theta_array[::step]
    zeta = particle.zeta_array[::step]
    rho = particle.rho_array[::step]
    ptheta = particle.ptheta_array[::step] / machine.psi_last.value
    pzeta = particle.pzeta_array[::step] / particle.initial_conditions.pzeta0 - 1
    energy = particle.energy_array[::step] / particle.initial_energy - 1

    FLUX_COLOR = "xkcd:forrest green"
    ANGLE_COLOR = "xkcd:royal blue"
    RHO_COLOR = "xkcd:crimson"
    PTHETA_COLOR = "xkcd:cerulean"
    PZETA_COLOR = "xkcd:cerulean"
    ENERGY_COLOR = "xkcd:dark orange"

    PLOT_KW = dict(linewidth=0, marker=".", markersize=2)

    axpsi.plot(t, psi, c=FLUX_COLOR, **PLOT_KW)
    axpsip.plot(t, psip, c=FLUX_COLOR, **PLOT_KW)
    axtheta.plot(t, theta, c=ANGLE_COLOR, **PLOT_KW)
    axzeta.plot(t, zeta, c=ANGLE_COLOR, **PLOT_KW)
    axrho.plot(t, rho, c=RHO_COLOR, **PLOT_KW)
    axptheta.plot(t, ptheta, c=PTHETA_COLOR, **PLOT_KW)
    axpzeta.plot(t, pzeta, c=PZETA_COLOR, **PLOT_KW)
    axenergy.plot(t, energy, c=ENERGY_COLOR, **PLOT_KW)

    MARGINS = (0, 0.03)
    for ax in axes.flatten():
        ax.margins(*MARGINS)

    for ax in axes[:, 0]:
        ax.yaxis.set_ticks_position("right")
        ax.yaxis.set_label_position("left")
    for ax in axes[:, 1]:
        ax.yaxis.set_ticks_position("right")
        ax.yaxis.set_label_position("left")

    LEFT_LABEL_KW = dict(fontsize=12)
    RIGHT_LABEL_KW = dict(fontsize=12)

    if machine.geometry is not None:
        # Possible divide by zero at t=0 during conversion
        with warnings.catch_warnings(action="ignore"):
            # FIXME: this does not rescale the time units as expected
            tu = cast(
                PlainQuantity[Array1],
                machine.quantity(t, "NormSecond").to("second").to_compact(),
            )
            t = tu.m
        tunits = tex_unit(tu)
        axrho.set_xlabel(rf"$t\ [{tunits}]$")
        axenergy.set_xlabel(rf"$t\ [{tunits}]$")
    else:
        axrho.set_xlabel(rf"$t\ [Normalized]$")
        axenergy.set_xlabel(rf"$t\ [Normalized]$")

    axpsi.set_ylabel(r"$\psi(t)/\psi_{LCFS}$", **LEFT_LABEL_KW)
    axtheta.set_ylabel(r"$\theta(t)\ [rad]$", **LEFT_LABEL_KW)
    axptheta.set_ylabel(r"$P_\theta(t)/\psi_{LCFS}$", **LEFT_LABEL_KW)
    axrho.set_ylabel(r"$\rho_{||}(t)\ [Normalized]$", **LEFT_LABEL_KW)
    axpsip.set_ylabel(r"$\psi_p(t)/\psi_{p,LCFS}$", **RIGHT_LABEL_KW)
    axzeta.set_ylabel(r"$\zeta(t)\ [rad]$", **RIGHT_LABEL_KW)
    axpzeta.set_ylabel(r"$\Delta P_\zeta(t)/P_{\zeta,0}$", **RIGHT_LABEL_KW)
    axenergy.set_ylabel(r"$\Delta E(t)/ E_0$", **RIGHT_LABEL_KW)

    if show:
        plt.show()

    return fig, axes


def plot_poloidal_drift(
    machine: Machine,
    obj: Particle | Queue,
    array_shape: ArrayShape = (300, 300),
    levels: int = 50,
    locator: Locator = "MaxN",
    show: bool = True,
) -> tuple[Figure, Axes]:
    r"""Plots a particle's drift on the $R-Z$ plane, overlaid on a contour plot of the Hamiltonian.

    Parameters
    ----------
    machine
        The machine in which the particle was integrated. It is used to convert specific
        quantities to SI units.
    obj
        The Particle or Queue containing the particles.
    levels
        The number of contour levels.
    show
        Whether or not to call `plt.show()`.

    Other parameters
    ----------------
    locator
        The tick locator to use to locate contour levels.

    Raises
    ------
    AttributeError
        If `machine` has not defined a `geometry`.

    """
    if isinstance(obj, Particle):
        particles = [obj]
    else:
        particles = obj.particles()

    geometry: GeometryObject = getattr(machine, "geometry")

    fig = plt.figure(figsize=(4, 3))
    ax = fig.add_subplot()

    flux = _resolve_magnetic_flux_kind(machine.bfield, "Toroidal")

    if flux == "Toroidal":
        flux_arg_name = "psi"
    else:
        flux_arg_name = "psip"

    CMAP = "plasma"
    LOG_LOCATOR_BASE = 1 + 1e-10
    LAST_COLOR = "k"
    MARGINS = (0.001, 0.001)

    for particle in particles:
        if particle.steps_stored == 0:
            continue
        if flux == "Toroidal":
            flux_array = particle.psi_array
        else:
            flux_array = particle.psip_array
        particle_eval_arg = {
            flux_arg_name: flux_array,
            "theta": particle.theta_array % TAU,
        }
        prlab = geometry.eval_rlab(**particle_eval_arg)
        pzlab = geometry.eval_zlab(**particle_eval_arg)
        color = orbit_color(particle.orbit_type)
        ax.plot(prlab, pzlab, c=color, linewidth=0, marker=".", markersize=1, zorder=10)

    # Plot energy contour only if all particles have the same Pζ and μ
    pzetas = np.asarray([particle.initial_conditions.pzeta0 for particle in particles])
    mus = np.asarray([particle.initial_conditions.mu0 for particle in particles])
    if np.all(pzetas == pzetas[0]) and np.all(mus == mus[0]):

        if flux == "Toroidal":
            energy_function = energy_of_psi_grid
            eval_flux_of_r_function = geometry.eval_psi_of_r
        else:
            energy_function = energy_of_psip_grid
            eval_flux_of_r_function = geometry.eval_psip_of_r

        r_array = np.linspace(0, machine.rlast, array_shape[1]) * 0.99999
        flux_array: Array1 = eval_flux_of_r_function(r_array)  # pyright: ignore
        theta_grid, flux_grid = create_poloidal_grid(
            theta_array=np.linspace(0, TAU, array_shape[0]),
            flux_array=flux_array,
        )
        grid_eval_arg = {"theta": theta_grid, flux_arg_name: flux_grid}
        rlab_grid = geometry.eval_rlab(**grid_eval_arg)
        zlab_grid = geometry.eval_zlab(**grid_eval_arg)
        energy_grid = cast(
            Array2,
            machine.quantity(
                energy_function(
                    machine,
                    pzetas[0],
                    mus[0],
                    theta_grid,
                    flux_grid,
                ),
                "NormJoule",
            )
            .to("keV")
            .m,
        )

        _locator = locator.lower()
        _locator = (
            LogLocator(base=LOG_LOCATOR_BASE, numticks=levels)
            if _locator == "log"
            else MaxNLocator(nbins=levels)
        )

        contourf = ax.contourf(
            rlab_grid,
            zlab_grid,
            energy_grid,
            levels=levels,
            locator=_locator,
            cmap=CMAP,
        )
        ax.contour(
            contourf,
            linewidths=0.1,
            colors="k",
        )
        fig.colorbar(contourf, label=r"$Energy\ [keV]$")

    ax.plot(geometry.rlab_last, geometry.zlab_last, color=LAST_COLOR)
    ax.scatter(
        geometry.raxis,
        geometry.zaxis,
        marker="+",
        c="k",
        s=40,
        zorder=10,
    )

    ax.set_xlabel(r"$R\ [m]$")
    ax.set_ylabel(r"$Z\ [m]$")
    ax.set_aspect("equal")
    ax.margins(*MARGINS)
    ax.grid(False)

    if show:
        plt.show()

    return fig, ax


def plot_pzeta_poincare(
    machine: Machine,
    queue: Queue,
    intersect_params: IntersectParams,
    *,
    color: bool = True,
    initial: bool = False,
    show: bool = True,
) -> tuple[Figure, Axes]:
    r"""Plots a $P_\zeta-\theta$ or $P_\zeta-\zeta$ Poincare plot.

    The kind of plot depends on the `intersection` field of `intersect_params`.

    Parameters
    ----------
    machine
        The machine in which the `queue` was integrated.
    queue
        The `Queue` containing the integrated particles.
    intersect_params
        The queue's intersection parameters.
    color
        Whether or not to color each orbit with a different color.
    initial
        Whether or not to plot each particle's initial point.
    show
        Whether or not to call `plt.show()`.
    """

    fig = plt.figure(layout="constrained", figsize=(5, 3))
    ax = fig.add_subplot()

    if color:
        ax.set_prop_cycle(cycler(color="brcmk"))
    else:
        ax.set_prop_cycle(cycler(color=["blue"]))

    if intersect_params.intersection == "ConstTheta":
        xlabel = r"$\zeta\ [rads]$"
        array_name = "zeta_array"
    else:
        xlabel = r"$\theta\ [rads]$"
        array_name = "theta_array"

    for p in queue.particles():
        if p.steps_stored == 0:
            continue
        angle_array = _pi_mod(p._r.get_array(array_name))
        pzeta_array = p.pzeta_array / machine.psip_last.value
        ax.plot(
            angle_array,
            pzeta_array,
            linewidth=0,
            marker=".",
            markersize=1.5,
            markeredgewidth=0,
        )
        if initial:
            angle0 = angle_array[0]
            pzeta0 = pzeta_array[0]
            ax.scatter(
                angle0,
                pzeta0,
                s=10,
                c="k",
                marker="x",
            )

    ax.margins(0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"$P_\zeta/\psi_{p,LCFS}$")

    if show:
        plt.show()

    return fig, ax


def plot_rz_poincare(
    machine: Machine,
    queue: Queue,
    intersect_params: IntersectParams,
    *,
    color: bool = True,
    initial: bool = False,
    show: bool = True,
) -> tuple[Figure, Axes]:
    r"""Plots an $R-Z$ Poincare plot.

    The kind of plot depends on the `intersection` field of `intersect_params`.

    Parameters
    ----------
    machine
        The machine in which the `queue` was integrated.
    queue
        The `Queue` containing the integrated particles.
    intersect_params
        The queue's intersection parameters.
    color
        Whether or not to color each orbit with a different color.
    initial
        Whether or not to plot each particle's initial point.
    show
        Whether or not to call `plt.show()`.

    Raises
    ------
    RuntimeError
        If `queue` was run with a 'ConstTheta' intersection.
    RuntimeError
        If `machine.geometry` is not defined.
    """

    if intersect_params.intersection == "ConstTheta":
        raise RuntimeError(
            "Cannot plot R-Z Poincare map with a 'ConstTheta' intersection"
        )
    if machine.geometry is None:
        raise RuntimeError("Geometry must be defined")

    fig = plt.figure(layout="constrained", figsize=(3.5, 3))
    ax = fig.add_subplot()

    if color:
        ax.set_prop_cycle(cycler(color="brcmk"))
    else:
        ax.set_prop_cycle(cycler(color=["blue"]))

    for p in queue.particles():
        if p.steps_stored == 0:
            continue
        theta_array = p.theta_array % TAU
        if machine.qfactor.psi_state == "Good":
            eval_arg = {"psi": p.psi_array, "theta": theta_array}
        else:
            eval_arg = {"psip": p.psip_array, "theta": theta_array}
        rlab_array = machine.geometry.eval_rlab(**eval_arg)
        zlab_array = machine.geometry.eval_zlab(**eval_arg)
        ax.plot(
            rlab_array,
            zlab_array,
            linewidth=0,
            marker=".",
            markersize=1.3,
            markeredgewidth=0,
        )
        if initial:
            rlab0 = rlab_array[0]
            zlab0 = zlab_array[0]
            ax.scatter(
                rlab0,
                zlab0,
                s=5,
                c="k",
                marker="x",
            )

    LAST_COLOR = "k"
    MARGINS = (0.01, 0.01)

    axis_point = (machine.geometry.raxis, machine.geometry.zaxis)
    ax.plot(machine.geometry.rlab_last, machine.geometry.zlab_last, color=LAST_COLOR)
    ax.scatter(*axis_point, marker="+", c="k", s=40, zorder=10)

    ax.set_aspect("equal")
    ax.margins(*MARGINS)
    ax.set_xlabel(r"$R\ [m]$")
    ax.set_ylabel(r"$Z\ [m]$")

    if show:
        plt.show()

    return fig, ax


def orbit_color(orbit_type: OrbitType) -> str:
    r"""Returns each orbit type's color string.

    Helps with being consistent with the coloring on different plots.
    """
    # In general, Confined = 'bright', Lost = 'deep'
    match orbit_type:
        case "Undefined":
            return "xkcd:coral"
        case "TrappedConfined":
            return "xkcd:bright red"
        case "TrappedLost":
            return "xkcd:deep red"
        case "CoPassingConfined":
            return "xkcd:bright blue"
        case "CoPassingLost":
            return "xkcd:deep blue"
        case "CuPassingConfined":
            return "xkcd:bright green"
        case "CuPassingLost":
            return "xkcd:deep green"
        case "Potato":
            return "xkcd:tan"
        case "Stagnated":
            return "xkcd:sky"
        case "Unclassified":
            return "xkcd:bright purple"
        case _:
            return "xkcd:indigo"


def _pi_mod(arr: Array1) -> Array1:
    """Mods an angle time series in the interval [-π, π]."""
    a: Array1 = np.mod(arr, 2 * np.pi)
    a = a - 2 * np.pi * (a > np.pi)
    return a

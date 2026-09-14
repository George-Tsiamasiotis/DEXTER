r"""Plotting functions for machine objects.

Methods are also available as command line project scripts.

Functions
---------
plot_qfactor
    Plots a [`QfactorObject`][dexter.QfactorObject]'s $q(\psi)$, $q(\psi_p)$, $\psi_p(\psi)$ and $\psi(\psi_p)$.
plot_current
    Plots a [`CurrentObject`][dexter.CurrentObject]'s $g$, $I$ and their derivatives, with respect to
    $\psi$ and $\psi_p$.
plot_bfield
    Plots a [`BfieldObject`][dexter.BfieldObject]'s $B$ and its derivatives on the $R-Z$ plane.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from dexter.machine.machine import Machine
from dexter.machine.base import MachineObject
from dexter.types import MagneticFluxKind, Interpolation1dType, Array

plt.rcParams["figure.dpi"] = 180
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["figure.autolayout"] = False
plt.rcParams["figure.constrained_layout.use"] = True

TAU = 2 * np.pi
PI = np.pi


def plot_qfactor(
    machine: Machine,
    points: int = 1000,
    data: bool = False,
    show: bool = True,
) -> tuple[Figure, tuple[Axes, Axes, Axes]]:
    r"""Plots a [`QfactorObject`][dexter.QfactorObject]'s $q(\psi)$, $q(\psi_p)$, $\psi_p(\psi)$ and $\psi(\psi_p)$.

    The derivatives $d\psi_p/d\psi$ and $(d\psi/d\psi_p)^{-1}$ are calculated independently from
    $q(\psi)$ and $q(\psi_p)$. It is a good sanity check to ensure the two curves coincide.

    Parameters
    ----------
    machine
        The machine containing the current object.

    Other parameters
    ----------------
    points
        The number of points on which to evaluate the flux in each plot.
    data
        Whether or not to plot the data points, if `#!python machine.current.machine_type == "Numerical"`.
    show
        Whether or not to call `plt.show()`.

    Example
    -------

    ``` py title="Qfactor plot"
    >>> machine = dex.Machine.FromNetcdf(path, "Akima", "Bicubic")
    >>> fig, ax = dex.plot_qfactor(machine, data=True)

    ```

    ``` sh title="From command line"
    dexter-plot-qfactor ./netcdf.nc -i Cubic

    ```

    """
    qfactor = machine.qfactor

    fig = plt.figure(figsize=(5, 4))
    axes = fig.subplots(2, 2)
    axqt: Axes = axes[0, 0]
    axqp: Axes = axes[0, 1]
    axpt: Axes = axes[1, 0]
    axtp: Axes = axes[1, 1]

    EMPTY_PLOT_COLOR = "#9E9E9E"
    QFACTOR_COLOR = "r"
    DERIV_COLOR = "b"
    DERIV_STYLE = (0, (3, 2))
    MARGINS = (0, 0.01)

    DATA_COLOR = "k"
    DATA_MARKER = "+"
    DATA_SIZE = 15
    DATA_LABEL = r"$data\ points$"
    DATA_ZORDER = 10

    # ======================================================= q(ψ)
    axqt.set_xlabel(r"$\psi/\psi_{LCFS}$")
    axqt.set_ylabel(r"$q(\psi)$")
    if qfactor.psi_state == "Good":
        psis = np.linspace(0, qfactor.psi_last.value, points) * 0.99999
        psis_norm = psis / qfactor.psi_last.value
        q_values = qfactor.eval_q(psi=psis)
        axqt.plot(
            psis_norm,
            q_values,
            c=QFACTOR_COLOR,
            label=rf"$q(\psi)$",
        )
        if qfactor.psip_state == "Good":
            psips = qfactor.eval_other(psi=psis)
            d_values = qfactor.eval_deriv_of_other(psip=psips)
            axqt.plot(
                psis_norm,
                d_values,
                c=DERIV_COLOR,
                linestyle=DERIV_STYLE,
                label=rf"$d\psi/d\psi_p$",
            )
        if data and qfactor.machine_type == "Numerical":
            # Arrays always exist if `machine_type == "Numerical"`
            q_data = qfactor.q_array  # pyright: ignore
            psi_data = qfactor.psi_array  # pyright: ignore
            psi_data_norm = psi_data / qfactor.psi_last.value
            axqt.scatter(
                psi_data_norm,
                q_data,
                c=DATA_COLOR,
                marker=DATA_MARKER,
                s=DATA_SIZE,
                zorder=DATA_ZORDER,
                label=DATA_LABEL,
            )
        axqt.margins(*MARGINS)
        axqt.legend()
        axqt.grid()
    else:
        axqt.set_facecolor(EMPTY_PLOT_COLOR)
        axqt.text(0.5, 0.5, ha="center", va="center", s=r"$q(\psi)\ is\ not\ defined$")

    # ======================================================= q(ψp)
    axqp.set_xlabel(r"$\psi_p/\psi_{p,LCFS}$")
    axqp.set_ylabel(r"$q(\psi_p)$")
    if qfactor.psip_state == "Good":
        psips = np.linspace(0, qfactor.psip_last.value, points) * 0.99999
        psips_norm = psips / qfactor.psip_last.value
        q_values = qfactor.eval_q(psip=psips)
        axqp.plot(
            psips_norm,
            q_values,
            c=QFACTOR_COLOR,
            label=rf"$q(\psi_p)$",
        )
        if qfactor.psi_state == "Good":
            psis = qfactor.eval_other(psip=psips)
            d_values = 1 / qfactor.eval_deriv_of_other(psi=psis)
            axqp.plot(
                psips_norm,
                d_values,
                c=DERIV_COLOR,
                linestyle=DERIV_STYLE,
                label=r"$(d\psi_p/d\psi)^{-1}$",
            )
        if data and qfactor.machine_type == "Numerical":
            # Arrays always exist if `machine_type == "Numerical"`
            q_data = qfactor.q_array  # pyright: ignore
            psip_data = qfactor.psip_array  # pyright: ignore
            psip_data_norm = psip_data / qfactor.psip_last.value
            axqp.scatter(
                psip_data_norm,
                q_data,
                c=DATA_COLOR,
                marker=DATA_MARKER,
                s=DATA_SIZE,
                zorder=DATA_ZORDER,
                label=DATA_LABEL,
            )
        axqp.margins(*MARGINS)
        axqp.legend()
        axqp.grid()
    else:
        axqp.set_facecolor(EMPTY_PLOT_COLOR)
        axqp.text(
            0.5, 0.5, ha="center", va="center", s=r"$q(\psi_p)\ is\ not\ defined$"
        )

    # ======================================================= ψp(ψ)
    axpt.set_xlabel(r"$\psi/\psi_{LCFS}$")
    axpt.set_ylabel(r"$\psi_p(\psi)/\psi_{p,LCFS}$")
    if qfactor.psi_state == "Good":
        fluxes = np.linspace(0, qfactor.psi_last.value, points) * 0.99999
        fluxes_norm = fluxes / qfactor.psi_last.value
        other_values = qfactor.eval_other(psi=fluxes)
        other_values_norm = other_values / qfactor.psip_last.value
        axpt.plot(
            fluxes_norm,
            other_values_norm,
            c=QFACTOR_COLOR,
            label=rf"$\psi_p(\psi)$",
        )
        if data and qfactor.machine_type == "Numerical":
            # Arrays always exist if `machine_type == "Numerical"`
            psi_data = qfactor.psi_array  # pyright: ignore
            psip_data = qfactor.psip_array  # pyright: ignore
            psi_data_norm = psi_data / qfactor.psi_last.value
            psip_data_norm = psip_data / qfactor.psip_last.value
            axpt.scatter(
                psi_data_norm,
                psip_data_norm,
                c=DATA_COLOR,
                marker=DATA_MARKER,
                s=DATA_SIZE,
                zorder=DATA_ZORDER,
                label=DATA_LABEL,
            )
        axpt.margins(*MARGINS)
        axpt.legend()
        axpt.grid()
    else:
        axpt.set_facecolor(EMPTY_PLOT_COLOR)
        axpt.text(
            0.5, 0.5, ha="center", va="center", s=r"$\psi_p(\psi)\ is\ not\ defined$"
        )

    # ======================================================= ψ(ψp)
    axtp.set_xlabel(r"$\psi_p/\psi_{p,LCFS}$")
    axtp.set_ylabel(r"$\psi(\psi_p)/\psi_{LCFS}$")
    if qfactor.psip_state == "Good":
        fluxes = np.linspace(0, qfactor.psip_last.value, points) * 0.99999
        fluxes_norm = fluxes / qfactor.psip_last.value
        other_values = qfactor.eval_other(psip=fluxes)
        other_values_norm = other_values / qfactor.psi_last.value
        axtp.plot(
            fluxes_norm,
            other_values_norm,
            c=QFACTOR_COLOR,
            label=rf"$\psi(\psi_p)$",
        )
        if data and qfactor.machine_type == "Numerical":
            # Arrays always exist if `machine_type == "Numerical"`
            psi_data = qfactor.psi_array  # pyright: ignore
            psip_data = qfactor.psip_array  # pyright: ignore
            psi_data_norm = psi_data / qfactor.psi_last.value
            psip_data_norm = psip_data / qfactor.psip_last.value
            axtp.scatter(
                psip_data_norm,
                psi_data_norm,
                c=DATA_COLOR,
                marker=DATA_MARKER,
                s=DATA_SIZE,
                zorder=DATA_ZORDER,
                label=DATA_LABEL,
            )
        axtp.margins(*MARGINS)
        axtp.legend()
        axtp.grid()
    else:
        axtp.set_facecolor(EMPTY_PLOT_COLOR)
        axtp.text(
            0.5, 0.5, ha="center", va="center", s=r"$\psi(\psi_p)\ is\ not\ defined$"
        )

    if show:
        plt.show()

    return fig, axes


def plot_current(
    machine: Machine,
    flux: MagneticFluxKind = "Toroidal",
    points: int = 1000,
    data: bool = False,
    show: bool = True,
) -> tuple[Figure, tuple[Axes, Axes]]:
    r"""Plots a [`CurrentObject`][dexter.CurrentObject]'s $g$, $I$ and their derivatives, with respect to
    $\psi$ and $\psi_p$.

    Parameters
    ----------
    machine
        The machine containing the current object.

    Other parameters
    ----------------
    flux
        The kind of magnetic flux with respect to which to plot. If the toroidal flux is not a
        good coordinate, the poloidal flux is attempted.
    points
        The number of points on which to evaluate the flux in each plot.
    data
        Whether or not to plot the data points, if `#!python machine.current.machine_type == "Numerical"`.
    show
        Whether or not to call `plt.show()`.

    Example
    -------

    ``` py title="Current plot"
    >>> machine = dex.Machine.FromNetcdf(path, "Akima", "Bicubic")
    >>> fig, ax = dex.plot_current(machine, data=True)

    ```

    ``` sh title="From command line"
    dexter-plot-current ./netcdf.nc -d -f Poloidal

    ```

    """
    current = machine.current

    fig = plt.figure(figsize=(7, 2.5))
    axes = fig.subplots(1, 2)
    axg: Axes = axes[0]
    axi: Axes = axes[1]
    twg = axg.twinx()
    twi = axi.twinx()

    flux = _resolve_magnetic_flux_kind(current, flux)

    if flux == "Toroidal":
        lcfs = machine.psi_last
        flux_arg_name = "psi"
        flux_tex = r"\psi"
        flux_last_tex = r"\psi_{LCFS}"
    else:
        lcfs = machine.psip_last
        flux_arg_name = "psip"
        flux_tex = r"\psi_p"
        flux_last_tex = r"\psi_{p,LCFS}"

    fluxes = np.linspace(0, lcfs.value, points) * 0.99999
    fluxes_norm = fluxes / lcfs.value
    eval_arg = {flux_arg_name: fluxes}

    g_values = current.eval_g(**eval_arg)
    i_values = current.eval_i(**eval_arg)
    g_deriv_values = current.eval_g_deriv(**eval_arg)
    i_deriv_values = current.eval_i_deriv(**eval_arg)

    if data and current.machine_type == "Numerical":
        # Arrays always exist if `machine_type == "Numerical"`
        g_data = current.g_array  # pyright: ignore
        i_data = current.i_array  # pyright: ignore
        if flux == "Toroidal":
            flux_data: Array1 = current.psi_array  # pyright: ignore
            flux_data_norm = flux_data / machine.psi_last.value
        else:
            flux_data: Array1 = current.psip_array  # pyright: ignore
            flux_data_norm = flux_data / machine.psip_last.value
        DATA_COLOR = "k"
        DATA_MARKER = "+"
        DATA_SIZE = 15
        DATA_LABEL = r"$data\ points$"
        DATA_ZORDER = 10
        axg.scatter(
            flux_data_norm,
            g_data,
            c=DATA_COLOR,
            marker=DATA_MARKER,
            s=DATA_SIZE,
            zorder=DATA_ZORDER,
            label=DATA_LABEL,
        )
        axi.scatter(
            flux_data_norm,
            i_data,
            c=DATA_COLOR,
            marker=DATA_MARKER,
            s=DATA_SIZE,
            zorder=DATA_ZORDER,
            label=DATA_LABEL,
        )

    CURRENT_COLOR = "r"
    CURRENT_DERIV_COLOR = "b"
    MARGINS = (0, 0.01)

    axg.plot(
        fluxes_norm,
        g_values,
        c=CURRENT_COLOR,
        label=rf"$g({flux_tex})$",
    )
    twg.plot(  # legend
        [],
        [],
        c=CURRENT_COLOR,
        label=rf"$g({flux_tex})$",
    )
    twg.plot(
        fluxes_norm,
        g_deriv_values,
        c=CURRENT_DERIV_COLOR,
        label=rf"$dg({flux_tex})/d{flux_tex}$",
    )

    axi.plot(
        fluxes_norm,
        i_values,
        c=CURRENT_COLOR,
        label=rf"$I({flux_tex})$",
    )
    twi.plot(  # legend
        [],
        [],
        c=CURRENT_COLOR,
        label=rf"$I({flux_tex})$",
    )
    twi.plot(
        fluxes_norm,
        i_deriv_values,
        c=CURRENT_DERIV_COLOR,
        label=rf"$dI({flux_tex})/d{flux_tex}$",
    )

    axg.grid()
    axi.grid()

    twg.legend()
    twi.legend()

    axg.margins(*MARGINS)
    axi.margins(*MARGINS)

    twg.margins(*MARGINS)
    twi.margins(*MARGINS)

    twg.spines["left"].set_color(CURRENT_COLOR)
    twi.spines["left"].set_color(CURRENT_COLOR)

    twg.spines["right"].set_color(CURRENT_DERIV_COLOR)
    twi.spines["right"].set_color(CURRENT_DERIV_COLOR)

    axg.yaxis.label.set_color(CURRENT_COLOR)
    axi.yaxis.label.set_color(CURRENT_COLOR)

    twg.yaxis.label.set_color(CURRENT_DERIV_COLOR)
    twi.yaxis.label.set_color(CURRENT_DERIV_COLOR)

    axg.tick_params("y", which="both", colors=CURRENT_COLOR)
    axi.tick_params("y", which="both", colors=CURRENT_COLOR)

    twg.tick_params("y", which="both", colors=CURRENT_DERIV_COLOR)
    twi.tick_params("y", which="both", colors=CURRENT_DERIV_COLOR)

    axg.set_xlabel(rf"${flux_tex}/{flux_last_tex}$")
    axi.set_xlabel(rf"${flux_tex}/{flux_last_tex}$")

    axg.set_ylabel(rf"$g({flux_tex})\ [Normalized]$")
    axi.set_ylabel(rf"$I({flux_tex})\ [Normalized]$")

    twg.set_ylabel(rf"$dg({flux_tex})/d{flux_tex}\ [Normalized]$")
    twi.set_ylabel(rf"$dI({flux_tex})/d{flux_tex}\ [Normalized]$")

    if show:
        plt.show()

    return fig, (axg, axi)


def plot_bfield(
    machine: Machine,
    levels: int = 15,
    show: bool = True,
) -> tuple[Figure, tuple[Axes, Axes, Axes]]:
    """Plots a [`BfieldObject`][dexter.BfieldObject]'s $B$ and its derivatives on the $R-Z$ plane.

    Parameters
    ----------
    machine
        The machine containing the current object.

    Other parameters
    ----------------
    levels
        The number of contour levels.
    show
        Whether or not to call `plt.show()`.

    Example
    -------

    ``` py title="Bfield plot"
    >>> machine = dex.Machine.FromNetcdf(path, "Akima", "Bicubic")
    >>> fig, ax = dex.plot_bfield(machine, levels=20)

    ```

    ``` sh title="From command line"
    dexter-plot-bfield ./netcdf.nc -l 30

    ```

    """
    bfield = machine.bfield
    geometry = machine.geometry
    if geometry is None:
        raise RuntimeError("'Geometry' must be defined")

    fig = plt.figure(figsize=(6.1, 4))
    axes = fig.subplots(2, 2)
    axbb: Axes = axes[0, 0]
    axbf: Axes = axes[0, 1]
    axbt: Axes = axes[1, 0]
    axsc: Axes = axes[1, 1]

    flux_points, theta_points = getattr(
        bfield, "shape", (300, 300)
    )  # for analytical bfields
    thetas = np.linspace(0, TAU, theta_points)
    if bfield.psi_state == "Good":
        psi_last = machine.psi_last.value
        fluxes = np.linspace(1e-10, psi_last, flux_points) * 0.99999
        flux_arg_name = "psi"
        flux_tex = "psi"
    else:
        psip_last = machine.psip_last.value
        fluxes = np.linspace(1e-10, psip_last, flux_points) * 0.99999
        flux_arg_name = "psip"
        flux_tex = "psi_p"

    flux_grid, theta_grid = np.meshgrid(fluxes, thetas)
    lab_eval_args = {flux_arg_name: flux_grid, "theta": theta_grid}
    rlabs = geometry.eval_rlab(**lab_eval_args)
    zlabs = geometry.eval_zlab(**lab_eval_args)

    # The B arrays must be "rotated" to account for θ padding
    theta_offset = abs(getattr(bfield, "padding_theta", 0))
    bfield_eval_args = lab_eval_args | {"theta": (theta_grid + theta_offset) % TAU}
    bb = bfield.eval_b(**bfield_eval_args)
    bf = bfield.eval_deriv_flux(**bfield_eval_args)
    bt = bfield.eval_deriv_theta(**bfield_eval_args)
    bb = machine.quantity(bb, "NormTesla").to("Tesla")

    CMAP = "plasma"
    LEVEL_LINE_COLOR = "k"
    LEVEL_LINE_WIDTH = 0.5
    SC_COLOR = "xkcd:bright blue"
    LAST_COLOR = "k"
    MARGINS = (0.01, 0.01)

    csbb = axbb.contourf(rlabs, zlabs, bb.m, cmap=CMAP, levels=levels)
    csbf = axbf.contourf(rlabs, zlabs, bf, cmap=CMAP, levels=levels)
    csbt = axbt.contourf(rlabs, zlabs, bt, cmap=CMAP, levels=levels)
    cssc = axsc.contourf(rlabs, zlabs, bt, cmap=CMAP, levels=levels)
    axbb.contour(csbb, colors="k", linewidths=LEVEL_LINE_WIDTH)
    axbf.contour(csbf, colors="k", linewidths=LEVEL_LINE_WIDTH)
    axbt.contour(csbt, colors="k", linewidths=LEVEL_LINE_WIDTH)
    axsc.contour(csbt, colors="k", linewidths=LEVEL_LINE_WIDTH)

    axsc.contour(
        cssc,
        levels=[0],
        colors=SC_COLOR,
        linewidths=2 * LEVEL_LINE_WIDTH,
        linestyle="solid",
        negative_linestyles="solid",
    )
    axsc.plot([], [], color=SC_COLOR, label=r"$dB(R,Z)/d\theta$")
    axsc.legend(loc="upper right", prop={"size": 7})

    plt.colorbar(csbb, label=r"$B(R, Z)\ [Tesla]$")
    plt.colorbar(csbf, label=rf"$dB(R, Z)/d\{flux_tex}$")
    plt.colorbar(csbt, label=rf"$dB(R, Z)/d\theta$")
    plt.colorbar(cssc, label=rf"$dB(R, Z)/d\theta$")

    axis_point = (geometry.raxis, geometry.zaxis)

    for ax in [axbb, axbf, axbt, axsc]:
        ax.plot(geometry.rlab_last, geometry.zlab_last, color=LAST_COLOR)
        ax.scatter(*axis_point, marker="+", c="k", s=40, zorder=10)
        ax.set_aspect("equal")
        ax.margins(*MARGINS)

    for lower_ax in axes[1, :]:
        lower_ax.set_xlabel(r"$R\ [m]$")
    for left_ax in axes[:, 0]:
        left_ax.set_ylabel(r"$Z\ [m]$")

    if show:
        plt.show()

    return fig, axes


def _resolve_magnetic_flux_kind(
    obj: MachineObject,
    passed: MagneticFluxKind,
) -> MagneticFluxKind:
    match passed:
        case "Toroidal":
            if obj.psi_state == "Good":
                return "Toroidal"
            elif obj.psip_state == "Good":
                print('psi_state = "Bad", using poloidal flux')
                return "Poloidal"
            else:
                assert False, "unreachable"
        case "Poloidal":
            if obj.psip_state == "Good":
                return "Poloidal"
            elif obj.psi_state == "Good":
                print('psip_state = "Bad", using toroidal flux')
                return "Toroidal"
            else:
                assert False, "unreachable"

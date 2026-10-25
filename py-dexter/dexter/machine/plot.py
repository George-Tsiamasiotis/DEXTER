r"""Plotting functions for machine objects.

Some methods are also available as command line project scripts.

Functions
---------
plot_qfactor
    Plots a [`QfactorObject`][dexter.QfactorObject]'s
    $q(\psi)$, $q(\psi_p)$, $q'(\psi)$, $q'(\psi_p)$, $\psi_p(\psi)$ and $\psi(\psi_p)$.
plot_current
    Plots a [`CurrentObject`][dexter.CurrentObject]'s $g$, $I$ and their derivatives, with respect
    to $\psi$ and $\psi_p$.
plot_bfield
    Plots a [`BfieldObject`][dexter.BfieldObject]'s $B$ and its derivatives on the $R-Z$ plane.
plot_mode
    Plots a [`ModeObject`][dexter.ModeObject]'s $\alpha$, $\phi$ and its derivatives with respect
    to $\psi$ or $\psi_p$.
plot_efield
    Plots an [`EfieldObject`][dexter.EfieldObject]'s $E$, $Φ$, with respect to $\psi$ or $\psi_p$.

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from typing import cast

from dexter._utils import tex_unit

from dexter.types import Array1, Array2, MagneticFluxKind
from dexter.machine.machine import Machine
from dexter.machine.base import MachineObject
from dexter.machine.qfactors import NcQfactor
from dexter.machine.currents import NcCurrent
from dexter.machine.modes import FluteMode, ModeObject
from dexter.comspace.curve import StationaryCurve

TAU = 2 * np.pi
PI = np.pi


def plot_qfactor(
    machine: Machine,
    points: int = 1000,
    data: bool = False,
    show: bool = True,
) -> tuple[Figure, tuple[Axes, Axes, Axes]]:
    r"""Plots a [`QfactorObject`][dexter.QfactorObject]'s
    $q(\psi)$, $q(\psi_p)$, $q'(\psi)$, $q'(\psi_p)$, $\psi_p(\psi)$ and $\psi(\psi_p)$.

    The derivatives $d\psi_p/d\psi$ and $(d\psi/d\psi_p)^{-1}$ are calculated independently from
    $q(\psi)$ and $q(\psi_p)$. It is a good sanity check to ensure the two curves coincide.

    Parameters
    ----------
    machine
        The machine containing the current object.
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

    fig = plt.figure(figsize=(5.5, 4))
    axes = fig.subplots(2, 2)
    axqt: Axes = axes[0, 0]
    axqp: Axes = axes[0, 1]
    axpt: Axes = axes[1, 0]
    axtp: Axes = axes[1, 1]

    EMPTY_PLOT_COLOR = "#9E9E9E"
    QFACTOR_COLOR = "r"
    QFACTOR_DERIV_COLOR = "xkcd:bright blue"
    QFACTOR_DERIV_ALPHA = 0.6
    DERIV_COLOR = "b"
    DERIV_STYLE = (0, (3, 2))
    MARGINS = (0.001, 0.01)

    DATA_COLOR = "k"
    DATA_MARKER = "."
    DATA_SIZE = 1
    DATA_LABEL = r"$data\ points$"
    DATA_ZORDER = 2

    # ======================================================= q(ψ)

    axqt.set_xlabel(r"$\psi/\psi_{LCFS}$")
    axqt.set_ylabel(r"$q(\psi)$")
    axqt.yaxis.label.set_color(QFACTOR_COLOR)
    axqt.tick_params("y", which="both", colors=QFACTOR_COLOR)
    axdqt = axqt.twinx()
    axdqt.set_ylabel(r"$dq(\psi)/d\psi$")
    axdqt.spines["left"].set_color(QFACTOR_COLOR)
    axdqt.spines["right"].set_color(QFACTOR_DERIV_COLOR)
    axdqt.yaxis.label.set_color(QFACTOR_DERIV_COLOR)
    axdqt.tick_params("y", which="both", colors=QFACTOR_DERIV_COLOR)
    if qfactor.psi_state == "Good":
        psis = np.linspace(0, qfactor.psi_last.value, points) * 0.99999
        psis_norm = psis / qfactor.psi_last.value
        q_values = qfactor.eval_q(psi=psis)
        dq_values = qfactor.eval_q_deriv(psi=psis)
        axqt.plot(
            psis_norm,
            q_values,
            c=QFACTOR_COLOR,
            label=rf"$q(\psi)$",
        )
        axdqt.plot(
            psis_norm,
            dq_values,
            c=QFACTOR_DERIV_COLOR,
            alpha=QFACTOR_DERIV_ALPHA,
        )
        # axdqt.plot(  # debug
        #     psis_norm, np.gradient(q_values, psis), c="k", linestyle="--"
        # )
        axqt.plot(  # legend
            [],
            [],
            c=QFACTOR_DERIV_COLOR,
            alpha=QFACTOR_DERIV_ALPHA,
            label=rf"$dq(\psi)/d\psi$",
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
            qfactor = cast(NcQfactor, qfactor)
            q_data = qfactor.q_array
            psi_data = qfactor.psi_array
            psi_data_norm = psi_data / qfactor.psi_last.value
            axqt.plot(
                psi_data_norm,
                q_data,
                c=DATA_COLOR,
                linewidth=0,
                marker=DATA_MARKER,
                markersize=DATA_SIZE,
                zorder=DATA_ZORDER,
                label=DATA_LABEL,
            )
        axqt.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
        axqt.margins(*MARGINS)
        axqt.legend()
        axqt.grid()
    else:
        axqt.set_facecolor(EMPTY_PLOT_COLOR)
        axqt.text(0.5, 0.5, ha="center", va="center", s=r"$q(\psi)\ is\ not\ defined$")

    # ======================================================= q(ψp)

    axqp.set_xlabel(r"$\psi_p/\psi_{p,LCFS}$")
    axqp.set_ylabel(r"$q(\psi_p)$")
    axqp.yaxis.label.set_color(QFACTOR_COLOR)
    axqp.tick_params("y", which="both", colors=QFACTOR_COLOR)
    axdqp = axqp.twinx()
    axdqp.set_ylabel(r"$dq(\psi_p)/d\psi_p$")
    axdqp.spines["left"].set_color(QFACTOR_COLOR)
    axdqp.spines["right"].set_color(QFACTOR_DERIV_COLOR)
    axdqp.yaxis.label.set_color(QFACTOR_DERIV_COLOR)
    axdqp.tick_params("y", which="both", colors=QFACTOR_DERIV_COLOR)
    if qfactor.psip_state == "Good":
        psips = np.linspace(0, qfactor.psip_last.value, points) * 0.99999
        psips_norm = psips / qfactor.psip_last.value
        q_values = qfactor.eval_q(psip=psips)
        dq_values = qfactor.eval_q_deriv(psip=psips)
        axqp.plot(
            psips_norm,
            q_values,
            c=QFACTOR_COLOR,
            label=rf"$q(\psi_p)$",
        )
        axdqp.plot(
            psips_norm,
            dq_values,
            c=QFACTOR_DERIV_COLOR,
            alpha=QFACTOR_DERIV_ALPHA,
        )
        # axdqp.plot(  # debug
        #     psips_norm, np.gradient(q_values, psips), c="k", linestyle="--"
        # )
        axqp.plot(  # legend
            [],
            [],
            c=QFACTOR_DERIV_COLOR,
            alpha=QFACTOR_DERIV_ALPHA,
            label=rf"$dq(\psi_p)/d\psi_p$",
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
            qfactor = cast(NcQfactor, qfactor)
            q_data = qfactor.q_array
            psip_data = qfactor.psip_array
            psip_data_norm = psip_data / qfactor.psip_last.value
            axqp.plot(
                psip_data_norm,
                q_data,
                c=DATA_COLOR,
                linewidth=0,
                marker=DATA_MARKER,
                markersize=DATA_SIZE,
                zorder=DATA_ZORDER,
                label=DATA_LABEL,
            )
        axqp.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
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
            assert isinstance(qfactor, NcQfactor)
            psi_data = qfactor.psi_array
            psip_data = qfactor.psip_array
            psi_data_norm = psi_data / qfactor.psi_last.value
            psip_data_norm = psip_data / qfactor.psip_last.value
            axpt.plot(
                psi_data_norm,
                psip_data_norm,
                c=DATA_COLOR,
                linewidth=0,
                marker=DATA_MARKER,
                markersize=DATA_SIZE,
                zorder=DATA_ZORDER,
                label=DATA_LABEL,
            )
        axpt.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
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
            assert isinstance(qfactor, NcQfactor)
            psi_data = qfactor.psi_array
            psip_data = qfactor.psip_array
            psi_data_norm = psi_data / qfactor.psi_last.value
            psip_data_norm = psip_data / qfactor.psip_last.value
            axtp.plot(
                psip_data_norm,
                psi_data_norm,
                c=DATA_COLOR,
                linewidth=0,
                marker=DATA_MARKER,
                markersize=DATA_SIZE,
                zorder=DATA_ZORDER,
                label=DATA_LABEL,
            )
        axtp.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
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

    The second order derivatives of $g$ and $I$ are plotted as well, but in arbitrary units and origin.

    Parameters
    ----------
    machine
        The machine containing the current object.
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
    g_deriv2_values = current.eval_g_deriv2(**eval_arg)
    i_deriv2_values = current.eval_i_deriv2(**eval_arg)

    if data and current.machine_type == "Numerical":
        current = cast(NcCurrent, current)
        g_data = current.g_array
        i_data = current.i_array
        if flux == "Toroidal":
            flux_data = cast(Array1, current.psi_array)
            flux_data_norm = flux_data / machine.psi_last.value
        else:
            flux_data = cast(Array1, current.psip_array)
            flux_data_norm = flux_data / machine.psip_last.value
        DATA_COLOR = "k"
        DATA_MARKER = "."
        DATA_SIZE = 1
        DATA_LABEL = r"$data\ points$"
        DATA_ZORDER = 3
        axg.plot(
            flux_data_norm,
            g_data,
            c=DATA_COLOR,
            linewidth=0,
            marker=DATA_MARKER,
            markersize=DATA_SIZE,
            zorder=DATA_ZORDER,
            label=DATA_LABEL,
        )
        axi.plot(
            flux_data_norm,
            i_data,
            c=DATA_COLOR,
            linewidth=0,
            marker=DATA_MARKER,
            markersize=DATA_SIZE,
            zorder=DATA_ZORDER,
            label=DATA_LABEL,
        )

    CURRENT_COLOR = "r"
    CURRENT_DERIV_COLOR = "b"
    CURRENT_DERIV2_COLOR = "xkcd:bright blue"
    CURRENT_DERIV2_ALPHA = 0.4
    MARGINS = (0.001, 0.01)

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

    # Plot derivatives in arbitrary units
    twg_lim = twg.get_ylim()
    twi_lim = twi.get_ylim()
    twg.autoscale(False)
    twi.autoscale(False)
    twg.plot(
        fluxes_norm,
        _normalize_array(g_deriv2_values, *twg_lim),
        c=CURRENT_DERIV2_COLOR,
        alpha=CURRENT_DERIV2_ALPHA,
        label=rf"$d^2g({flux_tex})/d{flux_tex}^2$",
    )
    twi.plot(
        fluxes_norm,
        _normalize_array(i_deriv2_values, *twi_lim),
        c=CURRENT_DERIV2_COLOR,
        alpha=CURRENT_DERIV2_ALPHA,
        label=rf"$d^2I({flux_tex})/d{flux_tex}^2$",
        zorder=-2,
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
    contour_points: bool = False,
    show: bool = True,
) -> tuple[Figure, tuple[Axes, Axes, Axes]]:
    r"""Plots a [`BfieldObject`][dexter.BfieldObject]'s $B$ and its derivatives on the $R-Z$ plane.

    Parameters
    ----------
    machine
        The machine containing the current object.
    levels
        The number of contour levels.
    contour_points
        Whether or not to plot the $\partial B/\partial \theta = 0$ points as calculated from the
        [`StationaryCurve`][dexter.StationaryCurve].
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

    bfield_eval_args = lab_eval_args | {"theta": theta_grid}
    bb = bfield.eval_b(**bfield_eval_args)
    bf = bfield.eval_deriv_flux(**bfield_eval_args)
    bt = bfield.eval_deriv_theta(**bfield_eval_args)
    bb = cast(Array2, machine.quantity(bb, "NormTesla").to("Tesla").m)

    CMAP = "plasma"
    LEVEL_LINE_WIDTH = 0.5
    SC_COLOR = "xkcd:bright blue"
    LAST_COLOR = "k"
    MARGINS = (0.01, 0.01)

    csbb = axbb.contourf(rlabs, zlabs, bb, cmap=CMAP, levels=levels)
    csbf = axbf.contourf(rlabs, zlabs, bf, cmap=CMAP, levels=levels)
    csbt = axbt.contourf(rlabs, zlabs, bt, cmap=CMAP, levels=levels)
    cssc = axsc.contourf(rlabs, zlabs, bt, cmap=CMAP, levels=levels)
    axbb.contour(csbb, colors="k", linewidths=LEVEL_LINE_WIDTH)
    axbf.contour(csbf, colors="k", linewidths=LEVEL_LINE_WIDTH)
    axbt.contour(csbt, colors="k", linewidths=LEVEL_LINE_WIDTH)
    axsc.contour(csbt, colors="k", linewidths=LEVEL_LINE_WIDTH)

    axsc.contour(  # stationary_curve from matplotlib
        cssc,
        levels=[0],
        colors="r",
        linewidths=2 * LEVEL_LINE_WIDTH,
        linestyle="solid",
        negative_linestyles="solid",
    )
    sc = StationaryCurve(machine)
    for seg in sc.segments:
        lab_eval_args = {flux_arg_name: seg.flux_array.values, "theta": seg.theta_array}
        rlab_seg = geometry.eval_rlab(**lab_eval_args)
        zlab_seg = geometry.eval_zlab(**lab_eval_args)
        axsc.plot(
            rlab_seg,
            zlab_seg,
            color=SC_COLOR,
            linewidth=2 * LEVEL_LINE_WIDTH,
            zorder=2,
        )
        if contour_points:
            axsc.scatter(
                rlab_seg,
                zlab_seg,
                marker="x",
                c="k",
                s=5,
                linewidths=0.4,
                alpha=0.8,
                zorder=1,
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


def plot_mode(
    mode: ModeObject,
    flux: MagneticFluxKind = "Toroidal",
    points: int = 1000,
    data: bool = False,
    show: bool = True,
) -> tuple[Figure, tuple[Axes, Axes]]:
    r"""Plots a [`ModeObject`][dexter.ModeObject]'s $\alpha$, $\phi$ and its derivatives with
    respect to $\psi$ or $\psi_p$.

    Parameters
    ----------
    mode
        The ModeObject.
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

    ``` py title="Mode plot"
    >>> mode = dex.NcFluteMode(path, "Cubic", 3, 2)
    >>> fig, ax = dex.plot_mode(mode, data=True)

    ```

    ``` sh title="From command line"
    dexter-plot-mode ./netcdf.nc 2 1

    ```

    """

    flux = _resolve_magnetic_flux_kind(mode, flux)

    fig = plt.figure(figsize=(5, 3))
    axes = fig.subplots(2, 1)
    axa: Axes = axes[0]
    axp: Axes = axes[1]
    twa = axa.twinx()
    twp = axp.twinx()

    m = mode.m
    n = mode.n
    sub = f"{m},{n}"

    if flux == "Toroidal":
        if isinstance(mode, FluteMode):
            psi_last = mode.lcfs.value
        else:
            psi_last: float = getattr(mode, "psi_array")[-1]
        fluxes = np.linspace(1e-10, psi_last, points)
        fluxes_norm = fluxes / psi_last
        flux_arg_name = "psi"
        flux_tex = r"\psi"
        flux_last_tex = r"\psi_{LCFS}"
    else:
        if isinstance(mode, FluteMode):
            psip_last = mode.lcfs.value
        else:
            psip_last: float = getattr(mode, "psip_array")[-1]
        fluxes = np.linspace(1e-10, psip_last, points)
        fluxes_norm = fluxes / psip_last
        flux_arg_name = "psip"
        flux_tex = r"\psi_p"
        flux_last_tex = r"\psi_{p,LCFS}"

    eval_args = {flux_arg_name: fluxes, "theta": 0, "zeta": 0, "t": 0}
    alpha = mode.eval_amplitude(**eval_args)
    phase = mode.eval_phase(**eval_args)
    dalpha = mode.eval_deriv_flux(**eval_args)
    dphase = np.gradient(phase)

    ALPHA_COLOR = "r"
    PHASE_COLOR = "xkcd:pumpkin"
    ALPHA_DERIV_COLOR = "b"
    PHASE_DERIV_COLOR = "xkcd:forrest green"

    axa.plot(
        fluxes_norm,
        alpha,
        color=ALPHA_COLOR,
    )
    axa_handle = twa.plot(
        [],
        [],
        color=ALPHA_COLOR,
        label=rf"$\alpha_{{{sub}}}({flux_tex})$",
    )
    twa_handle = twa.plot(
        fluxes_norm,
        dalpha,
        color=ALPHA_DERIV_COLOR,
        label=rf"$\partial\alpha_{{{sub}}}/\partial {flux_tex}$",
    )

    axa.yaxis.label.set_color(ALPHA_COLOR)
    twa.yaxis.label.set_color(ALPHA_DERIV_COLOR)
    twa.spines["left"].set_color(ALPHA_COLOR)
    twa.spines["right"].set_color(ALPHA_DERIV_COLOR)
    axa.tick_params("y", which="both", colors=ALPHA_COLOR)
    twa.tick_params("y", which="both", colors=ALPHA_DERIV_COLOR)

    axp.plot(
        fluxes_norm,
        phase,
        color=PHASE_COLOR,
    )
    axp_handle = twp.plot(
        [],
        [],
        color=PHASE_COLOR,
        label=rf"$\phi_{{{sub}}}({flux_tex})$",
    )
    twp_handle = twp.plot(
        fluxes_norm,
        dphase,
        color=PHASE_DERIV_COLOR,
        label=rf"$\partial\phi_{{{sub}}}/\partial {flux_tex}$",
    )

    axp.yaxis.label.set_color(PHASE_COLOR)
    twp.yaxis.label.set_color(PHASE_DERIV_COLOR)
    twp.spines["left"].set_color(PHASE_COLOR)
    twp.spines["right"].set_color(PHASE_DERIV_COLOR)
    axp.tick_params("y", which="both", colors=PHASE_COLOR)
    twp.tick_params("y", which="both", colors=PHASE_DERIV_COLOR)

    if data and mode.machine_type == "Numerical":
        if mode.psi_state == "Good":
            flux_data = getattr(mode, "psi_array")
        else:
            flux_data = getattr(mode, "psip_array")

        alpha_data = getattr(mode, "alpha_array")
        phase_data = getattr(mode, "phase_array")
        flux_data_norm = flux_data / flux_data[-1]

        DATA_COLOR = "k"
        DATA_MARKER = "."
        DATA_SIZE = 1
        DATA_LABEL = r"$data\ points$"
        DATA_ZORDER = 10

        axa_data_handle = axa.plot(
            flux_data_norm,
            alpha_data,
            linewidth=0,
            marker=DATA_MARKER,
            c=DATA_COLOR,
            markersize=DATA_SIZE,
            label=DATA_LABEL,
            zorder=DATA_ZORDER,
        )
        axp_data_handle = axp.plot(
            flux_data_norm,
            phase_data,
            marker=DATA_MARKER,
            linewidth=0,
            c=DATA_COLOR,
            markersize=DATA_SIZE,
            label=DATA_LABEL,
            zorder=DATA_ZORDER,
        )
    else:
        axa_data_handle = []
        axp_data_handle = []

    MARGINS = (0, 0.01)

    axa.grid()
    axp.grid()

    axa.margins(*MARGINS)
    twa.margins(*MARGINS)
    axp.margins(*MARGINS)
    twp.margins(*MARGINS)

    axa.axhline(0, c="k", linestyle="--", linewidth=1, zorder=3)
    twp.axhline(0, c="k", linestyle="--", linewidth=1, zorder=3)
    twa.set_ybound(upper=0.1)

    axa.set_ylabel(rf"$\alpha_{{{sub}}}({flux_tex})$")
    twa.set_ylabel(rf"$\partial\alpha_{{{sub}}}/\partial {flux_tex}$")
    axp.set_ylabel(rf"$\phi_{{{sub}}}({flux_tex})$")
    twp.set_ylabel(rf"$\partial\phi_{{{sub}}}/\partial {flux_tex}$")

    axp.set_xlabel(rf"${flux_tex}/{flux_last_tex}$")

    # Bring legend to front
    axa.legend(handles=axa_handle + twa_handle + axa_data_handle)
    axp.legend(handles=axp_handle + twp_handle + axp_data_handle)
    twa.legend_ = axa.get_legend()
    twp.legend_ = axp.get_legend()

    if show:
        plt.show()

    return fig, axes


def plot_efield(
    machine: Machine,
    flux: MagneticFluxKind = "Toroidal",
    points: int = 1000,
    show: bool = True,
) -> tuple[Figure, tuple[Axes, Axes]]:
    r"""Plots an [`EfieldObject`][dexter.EfieldObject]'s $E$, $Φ$, with respect to $\psi$ or $\psi_p$.

    Parameters
    ----------
    machine
        The machine containing the electric field object.
    flux
        The kind of magnetic flux with respect to which to plot. If the toroidal flux is not a
        good coordinate, the poloidal flux is attempted.
    points
        The number of points on which to evaluate the flux in each plot.
    show
        Whether or not to call `plt.show()`.

    Example
    -------

    ``` py title="Electric field plot"
    >>> machine = dex.Machine.FromNetcdf(path, "Akima", "Bicubic")
    >>> geometry = machine.geometry
    >>> LCFS = geometry.psi_last
    >>>
    >>> strength = machine.quantity(75, "kiloVolt / meter")
    >>> rmean = 0.96 * machine.rlast
    >>> mean = dex.MagneticFlux.Toroidal(float(machine.geometry.eval_psi_of_r(rmean)))
    >>> waist = LCFS.value * 0.02
    >>> efield = dex.GaussianEfield(strength, mean, waist)
    >>> machine.efield = efield
    >>>
    >>> fig, ax = dex.plot_current(machine, data=True)

    ```

    """
    efield = machine.efield

    fig = plt.figure(figsize=(7, 2.5))
    axes = fig.subplots(1, 2)
    axe: Axes = axes[0]
    axp: Axes = axes[1]
    twe = axe.twinx()
    twp = axp.twinx()

    flux = _resolve_magnetic_flux_kind(efield, flux)

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

    _e_values = efield.eval_e(**eval_arg)
    _phi_values = efield.eval_phi(**eval_arg)
    e_values = machine.quantity(_e_values, "NormVolt/NormMeter").to("kiloVolt/meter")
    phi_values = machine.quantity(_phi_values, "NormVolt").to("Volt")

    E_COLOR = "r"
    PHI_COLOR = "b"
    MARGINS = (0, 0.01)

    e_units_label = tex_unit(e_values)
    phi_units_label = tex_unit(phi_values)
    axe.plot(
        fluxes_norm,
        e_values.m,
        c=E_COLOR,
        label=rf"$E({flux_tex})\ [{e_units_label}]$",
    )
    axp.plot(
        fluxes_norm,
        phi_values.m,
        c=PHI_COLOR,
        label=rf"$\Phi({flux_tex})\ [{phi_units_label}]$",
    )

    # Set twin ax' units to normalized
    elim = machine.quantity(np.asarray(axe.get_ylim()), "kiloVolt / meter")
    elim_nu = tuple(elim.to("NormVolt / NormMeter").m)
    twe.set_ylim(elim_nu)
    philim = machine.quantity(np.asarray(axp.get_ylim()), "Volt")
    philim_nu = tuple(philim.to("NormVolt").m)
    twp.set_ylim(philim_nu)

    axe.grid()
    axp.grid()

    axe.legend()
    axp.legend()

    axe.margins(*MARGINS)
    axp.margins(*MARGINS)

    axe.set_xlabel(rf"${flux_tex}/{flux_last_tex}$")
    axp.set_xlabel(rf"${flux_tex}/{flux_last_tex}$")

    axe.set_ylabel(rf"$E({flux_tex})\ [{e_units_label}]$")
    axp.set_ylabel(rf"$\Phi({flux_tex})\ [{phi_units_label}]$")

    twe.set_ylabel(rf"$E({flux_tex})\ [Normalized]$")
    twp.set_ylabel(rf"$\Phi({flux_tex})\ [Normalized]$")

    if show:
        plt.show()

    return fig, (axe, axp)


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


def _normalize_array(a: np.ndarray, xmin: float, xmax: float) -> np.ndarray:
    r"""Normalizes an array within the `(xmin, xmax)` range."""
    amin = a.min()
    amax = a.max()
    if amin == amax:  # avoid division by 0
        return a
    else:
        return xmin + (a - amin) * (xmax - xmin) / (a.max() - amin)

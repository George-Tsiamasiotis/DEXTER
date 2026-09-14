from typing import assert_never

from dexter.machine.base import MachineObject
from dexter.types import MagneticFluxKind
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from dexter.machine.machine import Machine
from dexter.types import Interpolation1dType, Array

plt.rcParams["figure.dpi"] = 180
plt.rcParams["savefig.dpi"] = 300
plt.rcParams["figure.autolayout"] = False
plt.rcParams["figure.constrained_layout.use"] = True


def plot_current(
    machine: Machine,
    flux: MagneticFluxKind = "Toroidal",
    points: int = 1000,
    data: bool = False,
    show: bool = True,
) -> tuple[Figure, tuple[Axes, Axes]]:
    r"""Plots a CurrentObject's `g`, `I` and its derivatives with respect to `ψ` or `ψp`.

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
        Whether or not to plot the data points, if `machine.current.machine_type == "Numerical".
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
        values_key = "psi"
        flux_tex = r"\psi"
        flux_last_tex = r"\psi_{LCFS}"
    else:
        lcfs = machine.psip_last
        values_key = "psip"
        flux_tex = r"\psi_p"
        flux_last_tex = r"\psi_{p,LCFS}"

    fluxes = np.linspace(0, lcfs.value, points)
    fluxes_norm = fluxes / lcfs.value
    values = {values_key: fluxes}

    g_values = current.eval_g(**values)
    i_values = current.eval_i(**values)
    g_deriv_values = current.eval_g_deriv(**values)
    i_deriv_values = current.eval_i_deriv(**values)

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
        DATA_SIZE = 20
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

"""Definition of the `EnergyPzetaPlane` object and its plotting methods."""

from typing import cast

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from dexter.machine.machine import Machine

from dexter._core import _PyEnergyPzetaPlane
from dexter._utils import _ReprStrImpl


class EnergyPzetaPlane(_ReprStrImpl):
    r"""Representation of the 2D $(E, P_\zeta, \mu=const)$ plane.

    Parameters
    ----------
    machine
        The machine under study.
    mu
        The magnetic moment $\mu$.

    Example
    -------

    ``` py title="$E-P_\zeta$ plane construction"
    >>> geometry = dex.LarGeometry(2, 1.75, 0.5)
    >>> LCFS = geometry.psi_last
    >>> machine = dex.Machine(
    ...     geometry=geometry,
    ...     qfactor=dex.ParabolicQfactor(qaxis=1.1, qlast=3.5, lcfs=LCFS),
    ...     current=dex.LarCurrent(),
    ...     bfield=dex.LarBfield(),
    >>> )
    >>>
    >>> plane = dex.EnergyPzetaPlane(machine, mu=2e-5)
    >>> plane.plot()
    ```

    """

    _r: _PyEnergyPzetaPlane
    machine: Machine

    def __init__(self, machine: Machine, mu: float) -> None:
        self._r = _PyEnergyPzetaPlane(
            machine.qfactor._r, machine.current._r, machine.bfield._r, mu
        )
        assert mu == self._r.mu
        self.machine = machine  # handy

    def plot(
        self,
        xlim: tuple[float, float] = (-1.6, 0.5),
        ylim: float | tuple[float, float] = 3,
        show: bool = True,
    ) -> tuple[Figure, Axes]:
        r"""Plots the plane with the orbit classification curves.

        Parameters
        ----------
        xlim
            The xaxis limits, normalized to $\psi_{p,wall}$.
        ylim
            The yaxis limits, in $E/\mu$ units. If a float, it sets the yaxis upper limit while
            the lower limit is set to `0`.
        show
            Whether or not to call `plt.show()`.
        """

        fig = plt.figure(layout="constrained", figsize=(6, 4), dpi=150)
        ax = fig.add_subplot()

        match ylim:
            case _ if isinstance(ylim, float | int):
                ymin, ymax = 0, ylim
            case [ymin, ymax]:
                pass
            case _:
                raise RuntimeError("'ylim' must be a float or a 2-tuple of floats")

        PARABOLA_XAXIS_DENSITY = 1000

        psip_last = self.machine.psip_last.value
        mu = self._r.mu

        intercept = ymax * mu
        axis_intercepts = self._r.axis_parabola.horizontal_intercepts(intercept)
        lw_intercepts = self._r.left_wall_parabola.horizontal_intercepts(intercept)
        rw_intercepts = self._r.right_wall_parabola.horizontal_intercepts(intercept)

        axis_span = np.linspace(*axis_intercepts, PARABOLA_XAXIS_DENSITY)
        lw_span = np.linspace(*lw_intercepts, PARABOLA_XAXIS_DENSITY)
        rw_span = np.linspace(*rw_intercepts, PARABOLA_XAXIS_DENSITY)

        axis = self._r.axis_parabola.eval_array(axis_span) / mu
        lw = self._r.left_wall_parabola.eval_array(lw_span) / mu
        rw = self._r.right_wall_parabola.eval_array(rw_span) / mu

        tp_pzeta = self._r.tp_pzeta_values / psip_last  # not always a linspace
        tp_lower = self._r.tp_lower_values / mu
        tp_upper = self._r.tp_upper_values / mu

        LINEWIDTH = 3
        ZORDER = 5
        MA_COLOR = "xkcd:cobalt blue"
        LW_COLOR = "xkcd:coral"
        RW_COLOR = "xkcd:forrest green"
        TP_COLOR = "xkcd:bright pink"

        ax.plot(
            axis_span / psip_last,
            axis,
            c=MA_COLOR,
            linewidth=LINEWIDTH,
            zorder=ZORDER,
            label=r"$Magnetic\ Axis$",
        )
        ax.plot(
            lw_span / psip_last,
            lw,
            c=LW_COLOR,
            linewidth=LINEWIDTH,
            zorder=ZORDER,
            label=r"$Left\ Wall$",
        )
        ax.plot(
            rw_span / psip_last,
            rw,
            c=RW_COLOR,
            linewidth=LINEWIDTH,
            zorder=ZORDER,
            label=r"$Right\ Wall$",
        )
        ax.plot(
            np.concat((tp_pzeta, tp_pzeta[::-1], [-1])),
            np.concat((tp_lower, tp_upper[::-1], [tp_lower[0]])),
            c=TP_COLOR,
            linewidth=LINEWIDTH,
            zorder=ZORDER,
            label=r"$Trapped-Passing\ Boundary$",
        )

        if self.machine._reg is not None:
            si_ax = ax.twinx()
            max_energy_nu = ymax * mu
            max_energy_si = cast(
                float, self.machine.quantity(max_energy_nu, "NormJoule").to("keV").m
            )
            si_ax.set_ybound(0, max_energy_si)
            si_ax.set_ylabel(r"$E\ [keV]$")

        ax.grid(True)
        ax.set_xlim(*xlim)
        ax.set_ylim(ymin, ymax)
        ax.set_xlabel(r"$P_\zeta/\psi_{p, LCFS}$")
        ax.set_ylabel(r"$E/\mu$")
        ax.legend(loc="lower right")

        if show:
            plt.show()

        return fig, ax

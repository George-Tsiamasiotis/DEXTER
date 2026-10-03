"""Definition of the `EnergyPzetaPlane` object and its plotting methods."""

from typing import cast

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from dexter.machine.machine import Machine
from dexter.simulate.particle import Particle
from dexter.simulate.queue import Queue
from dexter.simulate.plot import orbit_color
from dexter.types import Array1

from dexter._core import _PyEnergyPzetaPlane
from dexter._utils import _ReprStrImpl


class _ParticlePoints:
    r"""Helper container type to store the particle's E-Pζ points and their orbit color.

    Parameters
    ----------
    obj
        The particle/queue to export the data from.
    """

    pzetas: Array1
    energies: Array1
    orbit_colors: list[str]

    def __init__(self, obj: Particle | Queue | None = None) -> None:
        if obj is None:  # Empty initialization
            self.pzetas = np.asarray([])
            self.energies = np.asarray([])
            self.orbit_colors = []
        elif isinstance(obj, Particle):  # Single point
            self.pzetas = np.atleast_1d(obj.initial_conditions.pzeta0)
            self.energies = np.atleast_1d(obj.initial_energy)
            self.orbit_colors = [orbit_color(obj.orbit_type)]
        else:
            self.pzetas = obj._r._initial_pzetas()
            self.energies = obj._r._initial_energies()
            orbit_types = obj._r._orbit_types()
            self.orbit_colors = [orbit_color(orbit_type) for orbit_type in orbit_types]

    def __len__(self) -> int:
        return len(self.pzetas)


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
    _particle_points: _ParticlePoints

    def __init__(self, machine: Machine, mu: float) -> None:
        self._r = _PyEnergyPzetaPlane(
            machine.qfactor._r, machine.current._r, machine.bfield._r, mu
        )
        assert mu == self._r.mu
        self.machine = machine  # handy reference
        self._particle_points = _ParticlePoints()

    def add_particles(self, obj: Particle | Queue):
        r"""Adds particles to the plane, discarding any previously stored ones.

        At the time only certain particle attributes are stored.

        Parameters
        ----------
        obj
            The particle/queue to extract the data from.
        """
        self._particle_points = _ParticlePoints(obj)

    def clear_particles(self):
        r"""Clears all stored particles."""
        self._particle_points = _ParticlePoints()

    def show(
        self,
        xlim: tuple[float, float] | None = None,
        ylim: float | tuple[float, float] = 3,
        particles: bool = True,
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
        ZORDER = 20
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

        if particles:
            zorder = 10
            if len(self._particle_points) < 500:
                point_size = 20
                zorder = 1000
            elif 500 <= len(self._particle_points) < 5000:
                point_size = 10
            else:
                point_size = 1
            ax.scatter(
                self._particle_points.pzetas / psip_last,
                self._particle_points.energies / mu,
                c=self._particle_points.orbit_colors,
                s=point_size,
                zorder=zorder,
            )

        ax.grid(True)
        ax.set_ylim(ymin, ymax)
        if xlim is None:
            ax.margins(x=0.05)
        else:
            ax.set_xlim(*xlim)
        ax.set_xlabel(r"$P_\zeta/\psi_{p, LCFS}$")
        ax.set_ylabel(r"$E/\mu$")
        ax.legend(loc="lower right")

        if show:
            plt.show()

        return fig, ax

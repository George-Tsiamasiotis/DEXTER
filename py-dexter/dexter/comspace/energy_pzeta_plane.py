"""Definition of the `EnergyPzetaPlane` object and its plotting methods."""

from typing import cast

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from dexter.machine.machine import Machine
from dexter.comspace.gc_hamiltonian import CriticalPoint, BoundaryExtremum
from dexter.simulate.particle import Particle
from dexter.simulate.queue import Queue
from dexter.simulate.plot import orbit_color
from dexter.types import Array1

from dexter._core import _PyEnergyPzetaPlane
from dexter._utils import _ReprStrImpl, _into_pyarray1f64


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

    Attributes
    ----------
    mu
        The magnetic moment $\mu$.
    critical_points
        The calculated critical points, calculated through the
        [`calculate_critical_points`][dexter.EnergyPzetaPlane.calculate_critical_points] method.
    boundary_extrema
        The calculated boundary extrema, calculated through the
        [`calculate_boundary_extrema`][dexter.EnergyPzetaPlane.calculate_boundary_extrema] method.

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
    ... )
    >>>
    >>> plane = dex.EnergyPzetaPlane(machine, mu=2e-5)
    >>> fig, ax = plane.show(curves=True)
    >>>
    >>> # Critical points and Boundary extrema
    >>> pzeta_array = np.linspace(-1.2, 0.3, 100) * machine.psip_last.value
    >>> plane.calculate_critical_points(pzeta_array)
    >>> plane.calculate_boundary_extrema(pzeta_array)
    >>> fig, ax = plane.show()

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

    def calculate_critical_points(self, pzeta_array: Array1):
        r"""Calculates the Hamiltonian's critical points over the `pzeta_array`.

        The critical points are stored in `#!python self` and can be accessed through
        `#!python self.critical_points`.

        Parameters
        ----------
        pzeta_array
            The $P_\zeta$ values, *in Normalized Units*. It usually makes more sense to normalize
            `pzeta_array` with respect to $\psi_{p, LCFS}$.
        """
        self._r.calculate_critical_points(
            qfactor=self.machine.qfactor._r,
            current=self.machine.current._r,
            bfield=self.machine.bfield._r,
            pzeta_array=_into_pyarray1f64(pzeta_array),
        )

    def calculate_boundary_extrema(self, pzeta_array: Array1):
        r"""Calculates the Hamiltonian's boundary extrema over the `pzeta_array`.

        The critical points are stored in `#!python self` and can be accessed through
        `#!python self.boundary_extrema`.

        Parameters
        ----------
        pzeta_array
            The $P_\zeta$ values, *in Normalized Units*. It usually makes more sense to normalize
            `pzeta_array` with respect to $\psi_{p, LCFS}$.
        """
        self._r.calculate_boundary_extrema(
            qfactor=self.machine.qfactor._r,
            current=self.machine.current._r,
            bfield=self.machine.bfield._r,
            pzeta_array=_into_pyarray1f64(pzeta_array),
        )

    @property
    def mu(self) -> float:
        return self._r.mu

    @property
    def critical_points(self) -> list[CriticalPoint]:
        return [CriticalPoint._wrap(_point) for _point in self._r.critical_points]

    @property
    def boundary_extrema(self) -> list[BoundaryExtremum]:
        return [BoundaryExtremum._wrap(_point) for _point in self._r.boundary_extrema]

    def _critical_points_pzetas(self) -> Array1:
        r"""Returns an array with the calculated critical points’ `Pζ`.

        The points are visited in order of creation, as in `self.critical_points`.
        """
        return self._r.critical_points_pzetas()

    def _critical_points_energies(self) -> Array1:
        r"""Returns an array with the calculated critical points’ energies.

        The points are visited in order of creation, as in `self.critical_points`.
        """
        return self._r.critical_points_energies()

    def _boundary_extrema_pzetas(self) -> Array1:
        r"""Returns an array with the calculated boundary extrema’s `Pζ`.

        The points are visited in order of creation, as in `self.boundary_extrema`.
        """
        return self._r.boundary_extrema_pzetas()

    def _boundary_extrema_energies(self) -> Array1:
        r"""Returns an array with the calculated boundary extrema’s energies.

        The points are visited in order of creation, as in `self.boundary_extrema`.
        """
        return self._r.boundary_extrema_energies()

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
        critical_points: bool = True,
        boundary_extrema: bool = True,
        curves: bool = False,
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
        particles
            Whether or not to plot the stored particles' $(E, P_\zeta)$ points.
        critical_points
            Whether or not to plot the stored critical points.
        boundary_extrema
            Whether or not to plot the stored boundary extrema.
        curves
            Whether or not to plot the analytical curves (parabolas and approximate
            trapped-passing boundary).
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

        psip_last = self.machine.psip_last.value
        mu = self._r.mu

        if self.machine._reg is not None:
            si_ax = ax.twinx()
            min_energy_nu = ymin * mu
            max_energy_nu = ymax * mu
            min_energy_si = cast(
                float, self.machine.quantity(min_energy_nu, "NormJoule").to("keV").m
            )
            max_energy_si = cast(
                float, self.machine.quantity(max_energy_nu, "NormJoule").to("keV").m
            )
            si_ax.set_ybound(min_energy_si, max_energy_si)
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

        _had_points = False

        if critical_points:
            pzetas = self._critical_points_pzetas() / psip_last
            energies = self._critical_points_energies() / mu
            if len(pzetas) > 0:
                _had_points = True
                ax.plot(
                    pzetas,
                    energies,
                    c="xkcd:purple blue",
                    linewidth=0,
                    marker=".",
                    markeredgewidth=0,
                    markersize=4,
                    zorder=20,
                    label=r"$Critical\ Points$",
                )

        if boundary_extrema:
            pzetas = self._boundary_extrema_pzetas() / psip_last
            energies = self._boundary_extrema_energies() / mu
            if len(pzetas) > 0:
                _had_points = True
                ax.plot(
                    pzetas,
                    energies,
                    c="xkcd:coral",
                    linewidth=0,
                    marker=".",
                    markeredgewidth=0,
                    markersize=4,
                    zorder=19,
                    label=r"$Boundary\ Extrema$",
                )

        if curves:
            PARABOLA_XAXIS_DENSITY = 1000

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
            ZORDER = 15
            MA_COLOR = "xkcd:aqua"
            LW_COLOR = "xkcd:gold"
            RW_COLOR = "xkcd:jade green"
            TP_COLOR = "xkcd:pinkish purple"
            CURVES_ALPHA = 0.5 if curves and _had_points else 1

            ax.plot(
                lw_span / psip_last,
                lw,
                c=LW_COLOR,
                linewidth=LINEWIDTH,
                zorder=ZORDER,
                label=r"$Left\ Wall$",
                alpha=CURVES_ALPHA,
            )
            ax.plot(
                rw_span / psip_last,
                rw,
                c=RW_COLOR,
                linewidth=LINEWIDTH,
                zorder=ZORDER,
                label=r"$Right\ Wall$",
                alpha=CURVES_ALPHA,
            )
            ax.plot(
                axis_span / psip_last,
                axis,
                c=MA_COLOR,
                linewidth=LINEWIDTH,
                zorder=ZORDER,
                label=r"$Magnetic\ Axis$",
                alpha=CURVES_ALPHA,
            )
            ax.plot(
                np.concat((tp_pzeta, tp_pzeta[::-1], [-1])),
                np.concat((tp_lower, tp_upper[::-1], [tp_lower[0]])),
                c=TP_COLOR,
                linewidth=LINEWIDTH,
                zorder=ZORDER,
                label=r"$Trapped-Passing\ Boundary$",
                alpha=CURVES_ALPHA,
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

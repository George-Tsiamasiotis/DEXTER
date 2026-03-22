"""Defines `Equilibrium`, a container type for equilibrium objects."""

import numpy as np
import matplotlib.pyplot as plt

from math import sqrt
from typing import Optional, Any

from dexter.types import Interp1DType, Interp2DType
from .objects import (
    LarBfield,
    NcGeometry,
    NcQfactor,
    NcCurrent,
    NcBfield,
    Geometry,
    Qfactor,
    Current,
    Bfield,
    Perturbation,
)


def numerical_equilibrium(
    path: str,
    interp1d_type: Interp1DType,
    interp2d_type: Interp2DType,
) -> Equilibrium:
    """Constructs a numerical equilibrium from a netCDF file.

    Perturbations are not constructed, but can be added manually in the `Equilibrium` object.

    Parameters
    ----------
    path
        Path to netCDF file.
    interp1d_type
        The 1D interpolation type.
    interp2d_type
        The 2D interpolation type.

    Returns
    -------
    Equilibrium
        Equilibrium object containing the corresponding
        [`NcGeometry`][dexter.NcGeometry],
        [`NcQfactor`][dexter.NcQfactor],
        [`NcCurrent`][dexter.NcCurrent] and
        [`NcBfield`][dexter.NcBfield],
    """
    geometry = NcGeometry(path, interp1d_type, interp2d_type)
    qfactor = NcQfactor(path, interp1d_type)
    current = NcCurrent(path, interp1d_type)
    bfield = NcBfield(path, interp2d_type)
    return Equilibrium(
        geometry=geometry,
        qfactor=qfactor,
        current=current,
        bfield=bfield,
        perturbation=None,  # TODO:
    )


class Equilibrium:
    """An Equilibrium.

    Contains all the necessary equilibrium objects and provides plots.

    For conveniently constructing a numerical equilibrium from a netCDF file, see
    [`numerical_equilibrium`][dexter.numerical_equilibrium]

    Parameters
    ----------
    qfactor
        The equilibrium's qfactor.
    current
        The equilibrium's current.
    bfield
        The equilibrium's bfield.
    perturbation
        todo: The equilibrium's perturbation, if they exist. Defaults to None (unperturbed equilibrium).
    geometry
        The equilibrium's geometry. This field is not required for calculations. Defaults to None.

    Example
    -------

    ```python title="Equilibrium creation"
    >>> eq = Equilibrium(
    ...     geometry=dex.ParabolicQfactor(1.1, 3.9, ("Toroidal", 0.45)),
    ...     qfactor=dex.UnityQfactor(),
    ...     current=dex.LarCurrent(),
    ...     bfield=dex.LarBfield(),
    ... )

    ```
    """

    _qfactor: Qfactor
    _current: Current
    _bfield: Bfield
    _perturbation: Perturbation
    _geometry: Optional[Geometry]

    def __init__(
        self,
        qfactor: Qfactor,
        current: Current,
        bfield: Bfield,
        perturbation: Optional[Perturbation] = None,
        geometry: Optional[Geometry] = None,
    ):
        self._geometry = geometry
        self._qfactor = qfactor
        self._current = current
        self._bfield = bfield
        if perturbation is None:
            self._perturbation = Perturbation([])
        else:
            self._perturbation = perturbation

    @property
    def geometry(self) -> Optional[Geometry]:
        """The equilibrium 'Geometry'."""
        return self._geometry

    @property
    def qfactor(self) -> Qfactor:
        """The equilibrium 'Qfactor'."""
        return self._qfactor

    @property
    def current(self) -> Current:
        """The equilibrium 'Current'."""
        return self._current

    @property
    def bfield(self) -> Bfield:
        """The equilibrium 'Bfield'."""
        return self._bfield

    @property
    def perturbation(self) -> Perturbation:
        """The equilibrium 'Perturbation'."""
        return self._perturbation

    @property
    def psi_wall(self) -> float:
        r"""The toroidal flux' value at the wall, $\psi_{wall}$."""
        return self._try_getattr("psi_wall")

    @property
    def psip_wall(self) -> float:
        r"""The poloidal flux' value at the wall, $\psi_{p,wall}$."""
        return self._try_getattr("psip_wall")

    @property
    def baxis(self) -> float:
        r"""The magnetic field strength on the axis, $B_{axis}$."""
        return self._try_getattr("baxis")

    @property
    def raxis(self) -> float:
        r"""The device's major radius, $R_{axis}$."""
        return self._try_getattr("raxis")

    @property
    def rwall(self) -> float:
        r"""The radial coordinate's value at the wall, $r_{wall}$."""
        return self._try_getattr("rwall")

    def _try_getattr(self, name: str) -> Any:
        """Tries to find the attribute 'name' in all mandatory fields, raising an AttributeError if no object
        has defined it.
        """
        try:
            return getattr(self.geometry, name)
        except AttributeError:
            try:
                return getattr(self.qfactor, name)
            except AttributeError:
                try:
                    return getattr(self.current, name)
                except:
                    try:
                        return getattr(self.bfield, name)
                    except:
                        raise AttributeError(
                            f"None of the equilibrium objects define '{name}'"
                        )

    @perturbation.setter
    def perturbation(self, perturbation: Perturbation):
        """Sets the equilibrium's Perturbation."""
        self._perturbation = perturbation

    def __str__(self) -> str:
        return (
            "Equilibrium:\n"
            + str(self.geometry)
            + str(self.qfactor)
            + str(self.current)
            + str(self.bfield)
            + str(self.perturbation)
        )

    def __repr__(self) -> str:
        return self.__str__()

    def plot_b(self, levels: int = 20):
        """Plots a contour plot of the magnetic field strength for a numerical equilibrium.

        Parameters
        ----------
        levels
            The number of contour levels. Defaults to 20.
        """

        fig = plt.figure(layout="constrained", dpi=120)
        ax = fig.add_subplot(aspect="equal")

        if (not isinstance(self.geometry, NcGeometry)) or (
            not isinstance(self.bfield, NcBfield)
        ):
            raise TypeError("'geometry' and 'bfield' must be netCDF objects")

        rlab_array = self.geometry.rlab_array
        zlab_array = self.geometry.zlab_array
        rlab_wall = self.geometry.rlab_wall
        zlab_wall = self.geometry.zlab_wall
        b_array = self.bfield.b_array

        ax.set_title(r"$Magnetic$ $field$ $strength$ $B$")
        ax.set_xlabel(r"$R[m]$")
        ax.set_ylabel(r"$Z[m]$")

        contour_kw = {"levels": levels, "cmap": "plasma"}

        contour = ax.contourf(rlab_array, zlab_array, b_array, **contour_kw)
        plt.colorbar(contour, ax=ax, label=r"$B[Normalized\quad Units]$")

        ax.plot(rlab_wall, zlab_wall, color="k", linewidth=2)
        ax.use_sticky_edges = False
        ax.margins(0.04, 0.04)

        plt.show()

    def plot_db(self, levels: int = 20):
        """Plots contour plots of the magnetic field's derivatives for a numerical equilibrium.

        Parameters
        ----------
        levels
            The number of contour levels. Defaults to 20.
        """

        if (not isinstance(self.geometry, NcGeometry)) or (
            not isinstance(self.bfield, NcBfield)
        ):
            raise TypeError("'geometry' and 'bfield' must be netCDF objects")

        fig = plt.figure(layout="constrained")
        subplot_kw = {"aspect": "equal"}
        ax = fig.subplots(1, 2, subplot_kw=subplot_kw)

        rlab_array = self.geometry.rlab_array
        zlab_array = self.geometry.zlab_array
        rlab_wall = self.geometry.rlab_wall
        zlab_wall = self.geometry.zlab_wall
        shape = self.geometry.shape

        contour_kw = {"levels": levels, "cmap": "plasma"}

        if self.geometry.psi_state == "Good":
            flux = "psi"
            wall = self.geometry.psi_wall
            db_dflux = self.bfield.db_dpsi
            db_dtheta = self.bfield.db_of_psi_dtheta
        elif self.geometry.psip_state == "Good":
            flux = "psi_p"
            wall = self.geometry.psip_wall
            db_dflux = self.bfield.db_dpsip
            db_dtheta = self.bfield.db_of_psip_dtheta
        else:
            raise Exception("unreachable")

        psi_grid, theta_grid = np.meshgrid(
            np.linspace(0, wall, shape[0]),
            np.linspace(0, 2 * np.pi, shape[1]),
            indexing="ij",
        )
        db_dflux_array = db_dflux(psi_grid, theta_grid)
        db_dtheta_array = db_dtheta(psi_grid, theta_grid)

        contour1 = ax[0].contourf(rlab_array, zlab_array, db_dflux_array, **contour_kw)
        contour2 = ax[1].contourf(rlab_array, zlab_array, db_dtheta_array, **contour_kw)
        plt.colorbar(contour1, ax=ax[0], label=rf"$dB/d\{flux}[Normalized\quad Units]$")
        plt.colorbar(contour2, ax=ax[1], label=r"$dB/d\theta[Normalized\quad Units]$")
        ax[0].plot(rlab_wall, zlab_wall, color="k", linewidth=2)
        ax[1].plot(rlab_wall, zlab_wall, color="k", linewidth=2)

        ax[0].set_title(rf"$dB/d\{flux}$")
        ax[1].set_title(r"$dB/d\theta$")
        ax[0].set_xlabel(r"$R[m]$")
        ax[1].set_xlabel(r"$R[m]$")
        ax[0].set_ylabel(r"$Z[m]$")
        ax[1].set_ylabel(r"$Z[m]$")
        ax[0].use_sticky_edges = False
        ax[1].use_sticky_edges = False
        ax[0].margins(0.04, 0.04)
        ax[1].margins(0.04, 0.04)

        plt.show()

    def plot_flux_surfaces(self, number: int = 20):
        """Plots the flux surfaces of a numerical equilibrium.

        Parameters
        ----------
        number
            The number of flux surfaces to (try to) plot. Defaults to 20.",
        """

        fig = plt.figure(figsize=(8, 7), layout="constrained", dpi=120)
        ax = fig.add_subplot(aspect="equal")

        if (not isinstance(self.geometry, NcGeometry)) or (
            not isinstance(self.bfield, NcBfield)
        ):
            raise TypeError("'geometry' and 'bfield' must be netCDF objects")

        rlab_array = self.geometry.rlab_array
        zlab_array = self.geometry.zlab_array
        rlab_wall = self.geometry.rlab_wall
        zlab_wall = self.geometry.zlab_wall

        shape = self.geometry.shape
        step = max([1, int(shape[0] / number)])
        print(f"Displaying {int(shape[0]/step)} surfaces.")
        for i in range(0, shape[0], step):
            ax.plot(rlab_array[i], zlab_array[i], color="blue", zorder=-1)

        # Cursor
        geom_center = (self.geometry.rgeo, self.geometry.zaxis)
        axis_point = (self.geometry.raxis, self.geometry.zaxis)

        def format_coord(x, y):
            r = sqrt((axis_point[0] - x) ** 2 + (axis_point[1] - y) ** 2)
            return f"(R, Z) = ({x:.5g}, {y:.5g}), r={r:.5g} "

        setattr(ax, "format_coord", format_coord)

        ax.set_title(r"$Poloidal\ flux\ surfaces$")
        ax.set_xlabel(r"$R[m]$")
        ax.set_ylabel(r"$Z[m]$")

        ax.plot(*axis_point, "ko", markersize=4, label="$R_{axis}$")
        ax.plot(*geom_center, "ro", markersize=4, label="$R_{geometric}$")

        ax.plot(rlab_wall, zlab_wall, color="k", linewidth=2)
        ax.legend()

        plt.show()

    def plot_boozer_theta(self, number: int = 80):
        r"""Plots the $\theta_B = const$ lines on a numerical equilibrium.

        Parameters
        ----------
        number
            The number of lines to (try to) plot. Defaults to 80.",
        """

        fig = plt.figure(figsize=(8, 7), layout="constrained", dpi=120)
        ax = fig.add_subplot(aspect="equal")

        if (not isinstance(self.geometry, NcGeometry)) or (
            not isinstance(self.bfield, NcBfield)
        ):
            raise TypeError("'geometry' and 'bfield' must be netCDF objects")

        rlab_array = self.geometry.rlab_array.T
        zlab_array = self.geometry.zlab_array.T
        rlab_wall = self.geometry.rlab_wall
        zlab_wall = self.geometry.zlab_wall

        shape = self.geometry.shape
        step = max([1, int(shape[1] / number)])
        print(f"Displaying {int(shape[1]/step)} surfaces.")
        for i in range(
            0, shape[1] - step, step
        ):  # last one lands too close to the first one
            ax.plot(rlab_array[i], zlab_array[i], color="blue", zorder=-1)

        # Cursor
        geom_center = (self.geometry.rgeo, self.geometry.zaxis)
        axis_point = (self.geometry.raxis, self.geometry.zaxis)

        def format_coord(x, y):
            r = sqrt((axis_point[0] - x) ** 2 + (axis_point[1] - y) ** 2)
            return f"(R, Z) = ({x:.5g}, {y:.5g}), r={r:.5g} "

        setattr(ax, "format_coord", format_coord)

        ax.set_title(r"$Boozer\ theta\ '\theta_B=const'\ lines$")
        ax.set_xlabel(r"$R[m]$")
        ax.set_ylabel(r"$Z[m]$")

        ax.plot(*axis_point, "ko", markersize=4, label="$R_{axis}$")
        ax.plot(*geom_center, "ro", markersize=4, label="$R_{geometric}$")

        ax.plot(rlab_wall, zlab_wall, color="k", linewidth=2)
        ax.legend()

        plt.show()

    def plot_midplane(self):
        r"""Plots $B$, $dB/d(\psi/\psi_p)$ and $dB/d\theta$ on the midplane.

        The midplane is defined as two $\theta=const$ lines:

        + $\theta = \pi$ and $\psi=[0, \psi_{wall}]$ (or $\psi_p=[0, \psi_{p,wall}]$)
        + $\theta = 0$ and $\psi=[0, \psi_{wall}]$ (or $\psi_p=[0, \psi_{p,wall}]$)
        """
        length = 1000

        if isinstance(self.bfield, LarBfield) or self.bfield.psi_state == "Good":
            b_of_flux = self.bfield.b_of_psi
            db_dflux = self.bfield.db_dpsi
            db_dtheta = self.bfield.db_of_psi_dtheta
            flux_wall = self.psi_wall
            flux_str = r"\psi"
            flux_wall_str = r"\psi_{wall}"
        else:
            b_of_flux = self.bfield.b_of_psip
            db_dflux = self.bfield.db_dpsip
            db_dtheta = self.bfield.db_of_psip_dtheta
            flux_wall = self.psip_wall
            flux_str = r"\psi_p"
            flux_wall_str = r"\psi_{p,wall}"

        fluxes = np.concat(
            (
                np.linspace(flux_wall, 0, length // 2),
                np.linspace(0, flux_wall, length // 2),
            )
        )
        thetas = np.concat(
            (
                np.full(length // 2, np.pi),
                np.full(length // 2, 0),
            ),
        )

        fig = plt.figure(figsize=(8, 7), layout="constrained", dpi=120)
        fig.suptitle(r"$Magnetic\ field\ quantities\ on\ midplane$")
        ax = fig.add_subplot()

        bs_of_flux = b_of_flux(fluxes, thetas)
        dbs_dflux = db_dflux(fluxes, thetas)
        dbs_dtheta = db_dtheta(fluxes, thetas)
        xs = np.linspace(-1, 1, length)
        ax.plot(
            xs,
            bs_of_flux,
            linewidth=2,
            color="b",
            label=rf"$B$",
        )
        ax.plot(
            xs,
            dbs_dflux,
            linewidth=2,
            color="r",
            label=rf"$dB/d{flux_str}$",
        )
        ax.plot(
            xs,
            dbs_dtheta,
            linewidth=2,
            color="k",
            label=rf"$dB/d\theta$",
        )

        ax.set_xlabel(rf"${flux_str} / {flux_wall_str}$")
        ax.set_ylabel(rf"$Amplitude$")
        ax.grid(True)
        ax.margins(0)
        ax.legend()

        plt.show()

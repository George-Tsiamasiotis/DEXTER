"""Defines `Equilibrium`, a container type for equilibrium objects."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

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
    *,
    padding: int = 15,
) -> Equilibrium:
    r"""Constructs a numerical equilibrium from a netCDF file.

    Perturbations are not constructed, but can be added manually in the `Equilibrium` object.

    Parameters
    ----------
    path
        Path to netCDF file.
    interp1d_type
        The 1D interpolation type.
    interp2d_type
        The 2D interpolation type.
    padding
        The left-right $\theta$ (per-side) padding width of the $B$ array. Defaults to 15.

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
    bfield = NcBfield(path, interp2d_type, padding=padding)
    return Equilibrium(
        geometry=geometry,
        qfactor=qfactor,
        current=current,
        bfield=bfield,
        perturbation=Perturbation([]),  # TODO:
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
    ...     geometry=dex.LarGeometry(baxis=2, raxis=1.75, rlast=0.5),
    ...     qfactor=dex.ParabolicQfactor(1.1, 3.9, dex.LastClosedFluxSurface("Toroidal", 0.45)),
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

    fig: Figure
    ax: Axes
    axes: Any  # list[Axes]

    def __init__(
        self,
        qfactor: Qfactor,
        current: Current,
        bfield: Bfield,
        perturbation: Optional[Perturbation] = Perturbation([]),
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
    def geometry(self) -> Geometry:
        """The equilibrium 'Geometry'."""
        if self._geometry is None:
            raise AttributeError("Geometry has not been defined")
        else:
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
    def psi_last(self) -> float:
        r"""The value of the last closed toroidal flux surface, $\psi_{LCFS}$."""
        return self._try_getattr("psi_last")

    @property
    def psip_last(self) -> float:
        r"""The value of the last closed poloidal flux surface, $\psi_{p,LCFS}$."""
        return self._try_getattr("psip_last")

    @property
    def baxis(self) -> float:
        r"""The magnetic field strength on the axis, $B_{axis}$."""
        return self._try_getattr("baxis")

    @property
    def raxis(self) -> float:
        r"""The device's major radius, $R_{axis}$."""
        return self._try_getattr("raxis")

    @property
    def rlast(self) -> float:
        r"""The radial coordinate's value at the last closed flux surface, $r_{LCFS}$."""
        return self._try_getattr("rlast")

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
        string = "Equilibrium:\n"
        string += str(getattr(self, "geometry", "")) + "\n"
        string += str(getattr(self, "qfactor", "")) + "\n"
        string += str(getattr(self, "current", "")) + "\n"
        string += str(getattr(self, "bfield", "")) + "\n"
        string += str(getattr(self, "perturbation", "")) + "\n"
        return string

    def __repr__(self) -> str:
        return self.__str__()

    def plot_b(self, levels: int = 20, show: bool = True):
        """Plots a contour plot of the magnetic field strength for a numerical equilibrium.

        Parameters
        ----------
        levels
            The number of contour levels. Defaults to 20.
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """

        self.fig = plt.figure(layout="constrained", dpi=120)
        self.ax = self.fig.add_subplot(aspect="equal")

        if (not isinstance(self.geometry, NcGeometry)) or (
            not isinstance(self.bfield, NcBfield)
        ):
            raise TypeError("'geometry' and 'bfield' must be netCDF objects")

        rlab_array = self.geometry.rlab_array
        zlab_array = self.geometry.zlab_array
        rlab_last = self.geometry.rlab_last
        zlab_last = self.geometry.zlab_last
        b_array = self.bfield.b_array

        self.ax.set_title(r"$Magnetic$ $field$ $strength$ $B$")
        self.ax.set_xlabel(r"$R[m]$")
        self.ax.set_ylabel(r"$Z[m]$")

        contour_kw = {"levels": levels, "cmap": "plasma"}

        contour = self.ax.contourf(rlab_array, zlab_array, b_array, **contour_kw)
        plt.colorbar(contour, ax=self.ax, label=r"$B[Normalized\quad Units]$")

        self.ax.plot(rlab_last, zlab_last, color="k", linewidth=2)
        self.ax.use_sticky_edges = False
        self.ax.margins(0.04, 0.04)

        if show:
            plt.show()

    def plot_db(self, levels: int = 20, show: bool = True):
        """Plots contour plots of the magnetic field's derivatives for a numerical equilibrium.

        Parameters
        ----------
        levels
            The number of contour levels. Defaults to 20.
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """

        if (not isinstance(self.geometry, NcGeometry)) or (
            not isinstance(self.bfield, NcBfield)
        ):
            raise TypeError("'geometry' and 'bfield' must be netCDF objects")

        self.fig = plt.figure(layout="constrained")
        subplot_kw = {"aspect": "equal"}
        self.axes = self.fig.subplots(1, 2, subplot_kw=subplot_kw)

        rlab_array = self.geometry.rlab_array
        zlab_array = self.geometry.zlab_array
        rlab_last = self.geometry.rlab_last
        zlab_last = self.geometry.zlab_last
        shape = self.geometry.shape

        contour_kw = {"levels": levels, "cmap": "plasma"}

        if self.geometry.psi_state == "Good":
            flux = "psi"
            lcfs = self.geometry.psi_last
            db_dflux = self.bfield.db_dpsi
            db_dtheta = self.bfield.db_of_psi_dtheta
        elif self.geometry.psip_state == "Good":
            flux = "psi_p"
            lcfs = self.geometry.psip_last
            db_dflux = self.bfield.db_dpsip
            db_dtheta = self.bfield.db_of_psip_dtheta
        else:
            raise Exception("unreachable")

        psi_grid, theta_grid = np.meshgrid(
            np.linspace(0, lcfs, shape[0]),
            np.linspace(0, 2 * np.pi, shape[1]),
            indexing="ij",
        )
        db_dflux_array = db_dflux(psi_grid, theta_grid)
        db_dtheta_array = db_dtheta(psi_grid, theta_grid)

        contour1 = self.axes[0].contourf(
            rlab_array, zlab_array, db_dflux_array, **contour_kw
        )
        contour2 = self.axes[1].contourf(
            rlab_array, zlab_array, db_dtheta_array, **contour_kw
        )
        plt.colorbar(
            contour1, ax=self.axes[0], label=rf"$dB/d\{flux}[Normalized\quad Units]$"
        )
        plt.colorbar(
            contour2, ax=self.axes[1], label=r"$dB/d\theta[Normalized\quad Units]$"
        )
        self.axes[0].plot(rlab_last, zlab_last, color="k", linewidth=2)
        self.axes[1].plot(rlab_last, zlab_last, color="k", linewidth=2)

        self.axes[0].set_title(rf"$dB/d\{flux}$")
        self.axes[1].set_title(r"$dB/d\theta$")
        self.axes[0].set_xlabel(r"$R[m]$")
        self.axes[1].set_xlabel(r"$R[m]$")
        self.axes[0].set_ylabel(r"$Z[m]$")
        self.axes[1].set_ylabel(r"$Z[m]$")
        self.axes[0].use_sticky_edges = False
        self.axes[1].use_sticky_edges = False
        self.axes[0].margins(0.04, 0.04)
        self.axes[1].margins(0.04, 0.04)

        if show:
            plt.show()

    def plot_flux_surfaces(self, number: int = 20, show: bool = True):
        """Plots the flux surfaces of a numerical equilibrium.

        Parameters
        ----------
        number
            The number of flux surfaces to (try to) plot. Defaults to 20.",
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """

        self.fig = plt.figure(figsize=(8, 7), layout="constrained", dpi=120)
        self.ax = self.fig.add_subplot(aspect="equal")

        if (not isinstance(self.geometry, NcGeometry)) or (
            not isinstance(self.bfield, NcBfield)
        ):
            raise TypeError("'geometry' and 'bfield' must be netCDF objects")

        rlab_array = self.geometry.rlab_array
        zlab_array = self.geometry.zlab_array
        rlab_last = self.geometry.rlab_last
        zlab_last = self.geometry.zlab_last

        shape = self.geometry.shape
        step = max([1, int(shape[0] / number)])
        print(f"Displaying {int(shape[0]/step)} surfaces.")
        for i in range(0, shape[0], step):
            self.ax.plot(rlab_array[i], zlab_array[i], color="blue", zorder=-1)

        # Cursor
        geom_center = (self.geometry.rgeo, self.geometry.zaxis)
        axis_point = (self.geometry.raxis, self.geometry.zaxis)

        def format_coord(x, y):
            r = sqrt((axis_point[0] - x) ** 2 + (axis_point[1] - y) ** 2)
            return f"(R, Z) = ({x:.5g}, {y:.5g}), r={r:.5g} "

        setattr(self.ax, "format_coord", format_coord)

        self.ax.set_title(r"$Poloidal\ flux\ surfaces$")
        self.ax.set_xlabel(r"$R[m]$")
        self.ax.set_ylabel(r"$Z[m]$")

        self.ax.plot(*axis_point, "ko", markersize=4, label="$R_{axis}$")
        self.ax.plot(*geom_center, "ro", markersize=4, label="$R_{geometric}$")

        self.ax.plot(rlab_last, zlab_last, color="k", linewidth=2)
        self.ax.legend()

        if show:
            plt.show()

    def plot_boozer_theta(self, number: int = 80, show: bool = True):
        r"""Plots the $\theta_B = const$ lines on a numerical equilibrium.

        Parameters
        ----------
        number
            The number of lines to (try to) plot. Defaults to 80.",
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """

        self.fig = plt.figure(figsize=(8, 7), layout="constrained", dpi=120)
        self.ax = self.fig.add_subplot(aspect="equal")

        if (not isinstance(self.geometry, NcGeometry)) or (
            not isinstance(self.bfield, NcBfield)
        ):
            raise TypeError("'geometry' and 'bfield' must be netCDF objects")

        rlab_array = self.geometry.rlab_array.T
        zlab_array = self.geometry.zlab_array.T
        rlab_last = self.geometry.rlab_last
        zlab_last = self.geometry.zlab_last

        shape = self.geometry.shape
        step = max([1, int(shape[1] / number)])
        print(f"Displaying {int(shape[1]/step)} surfaces.")
        for i in range(
            0, shape[1] - step, step
        ):  # last one lands too close to the first one
            self.ax.plot(rlab_array[i], zlab_array[i], color="blue", zorder=-1)

        # Cursor
        geom_center = (self.geometry.rgeo, self.geometry.zaxis)
        axis_point = (self.geometry.raxis, self.geometry.zaxis)

        def format_coord(x, y):
            r = sqrt((axis_point[0] - x) ** 2 + (axis_point[1] - y) ** 2)
            return f"(R, Z) = ({x:.5g}, {y:.5g}), r={r:.5g} "

        setattr(self.ax, "format_coord", format_coord)

        self.ax.set_title(r"$Boozer\ theta\ '\theta_B=const'\ lines$")
        self.ax.set_xlabel(r"$R[m]$")
        self.ax.set_ylabel(r"$Z[m]$")

        self.ax.plot(*axis_point, "ko", markersize=4, label="$R_{axis}$")
        self.ax.plot(*geom_center, "ro", markersize=4, label="$R_{geometric}$")

        self.ax.plot(rlab_last, zlab_last, color="k", linewidth=2)
        self.ax.legend()

        if show:
            plt.show()

    def plot_midplane(self, show: bool = True):
        r"""Plots $B$, $dB/d(\psi/\psi_p)$ and $dB/d\theta$ on the midplane.

        The midplane is defined as two $\theta=const$ lines:

        + $\theta = \pi$ and $\psi=[0, \psi_{LCFS}]$ (or $\psi_p=[0, \psi_{p,LCFS}]$)
        + $\theta = 0$ and $\psi=[0, \psi_{LCFS}]$ (or $\psi_p=[0, \psi_{p,LCFS}]$)

        Parameters
        ----------
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """
        length = 1000

        if isinstance(self.bfield, LarBfield) or self.bfield.psi_state == "Good":
            b_of_flux = self.bfield.b_of_psi
            db_dflux = self.bfield.db_dpsi
            db_dtheta = self.bfield.db_of_psi_dtheta
            lcfs = self.psi_last
            flux_str = r"\psi"
            lcfs_str = r"\psi_{LCFS}"
        else:
            b_of_flux = self.bfield.b_of_psip
            db_dflux = self.bfield.db_dpsip
            db_dtheta = self.bfield.db_of_psip_dtheta
            lcfs = self.psip_last
            flux_str = r"\psi_p"
            lcfs_str = r"\psi_{p,LCFS}"

        fluxes = np.concat(
            (
                np.linspace(lcfs, 0, length // 2),
                np.linspace(0, lcfs, length // 2),
            )
        )
        thetas = np.concat(
            (
                np.full(length // 2, np.pi),
                np.full(length // 2, 0),
            ),
        )

        self.fig = plt.figure(figsize=(8, 7), layout="constrained", dpi=120)
        self.fig.suptitle(r"$Magnetic\ field\ quantities\ on\ midplane$")
        self.ax = self.fig.add_subplot()

        bs_of_flux = b_of_flux(fluxes, thetas)
        dbs_dflux = db_dflux(fluxes, thetas)
        dbs_dtheta = db_dtheta(fluxes, thetas)
        xs = np.linspace(-1, 1, length)
        self.ax.plot(
            xs,
            bs_of_flux,
            linewidth=2,
            color="b",
            label=rf"$B$",
        )
        self.ax.plot(
            xs,
            dbs_dflux,
            linewidth=2,
            color="r",
            label=rf"$dB/d{flux_str}$",
        )
        self.ax.plot(
            xs,
            dbs_dtheta,
            linewidth=2,
            color="k",
            label=rf"$dB/d\theta$",
        )

        self.ax.set_xlabel(rf"${flux_str} / {lcfs_str}$")
        self.ax.set_ylabel(rf"$Amplitude$")
        self.ax.grid(True)
        self.ax.margins(0)
        self.ax.legend()

        if show:
            plt.show()

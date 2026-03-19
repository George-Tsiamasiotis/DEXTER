"""Defines `Equilibrium`, a container type for equilibrium objects."""

import numpy as np
import matplotlib.pyplot as plt

from typing import Optional, Any

from dexter.types import Interp1DType, Interp2DType
from .objects import (
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

    def _try_getattr(self, name: str) -> Any:
        """Tries to find the attribute 'name' in all mandatory fields, raising an AttributeError if no object
        has defined it.
        """
        try:
            return getattr(self.qfactor, name)
        except AttributeError:
            try:
                return getattr(self.current, name)
            except AttributeError:
                try:
                    return getattr(self.current, name)
                except:
                    raise AttributeError(
                        f"None of the equilibrium objects define '{name}'"
                    )

    @perturbation.setter
    def perturbation(self, perturbation: Perturbation):
        """Sets the equilibrium's Perturbation."""
        self._perturbation = perturbation

    def plot_b(self, levels: int = 20):
        """Plots a contour plot of the magnetic field strength.

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
        """Plots contour plots of the magnetic field's derivatives.

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

        db_dflux_array = np.full(shape, np.nan)
        db_dtheta_array = np.full(shape, np.nan)
        psi_grid, theta_grid = np.meshgrid(
            np.linspace(0, wall, shape[0]),
            np.linspace(0, 2 * np.pi, shape[1]),
            indexing="ij",
        )
        for i in range(shape[0]):
            for j in range(shape[1]):
                db_dflux_array[i, j] = db_dflux(psi_grid[i, j], theta_grid[i, j])
                db_dtheta_array[i, j] = db_dtheta(psi_grid[i, j], theta_grid[i, j])

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

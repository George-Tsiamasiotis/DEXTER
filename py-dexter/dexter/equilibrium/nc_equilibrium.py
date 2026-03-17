"""Defines various helper functions associated with the 'dexter-equilibrium' crate."""

import numpy as np
import matplotlib.pyplot as plt

from dexter.types import Interp1DType, Interp2DType
from .objects import (
    NcGeometry,
    NcQfactor,
    NcCurrent,
    NcBfield,
)


class NcEquilibrium:
    """Numerical equilibrium from an netCDF file.

    This object is simply a container for the corresponding equilibrium objects, but also offers plots that
    need more than one object.

    Parameters
    ----------
    path
        The path to the NetCDF file.
    interp1d_type
        The type of Interpolation of the 1D quantities.
    interp2d_type
        The type of Interpolation of the 2D quantities.

    Attributes
    ----------
    geometry
        The equilibrium's geometry.
    qfactor
        The equilibrium's qfactor.
    current
        The equilibrium's current.
    bfield
        The equilibrium's bfield.
    perturbation
        todo: The equilibrium's perturbation, if they exist.

    Example
    -------

    ```python title="NcEquilibrium creation"
    >>> eq = dex.NcEquilibrium(path, "Cubic", "Bicubic")
    >>> (geometry, qfactor, current, bfield) = eq.objects()

    ```
    """

    geometry: NcGeometry
    qfactor: NcQfactor
    current: NcCurrent
    bfield: NcBfield
    perturbation: None

    def __init__(
        self,
        path: str,
        interp1d_type: Interp1DType,
        interp2d_type: Interp2DType,
    ):
        self.geometry = NcGeometry(path, interp1d_type, interp2d_type)
        self.qfactor = NcQfactor(path, interp1d_type)
        self.current = NcCurrent(path, interp1d_type)
        self.bfield = NcBfield(path, interp2d_type)
        self.perturbation = None  # TODO:

    def objects(self) -> tuple[NcGeometry, NcQfactor, NcCurrent, NcBfield]:
        """Returns a tuple with the constructed objects."""
        return (self.geometry, self.qfactor, self.current, self.bfield)

    def plot_b(self, levels: int = 20):
        """Plots a contour plot of the magnetic field strength.

        Parameters
        ----------
        levels
            The number of contour levels. Defaults to 20.
        """

        fig = plt.figure(layout="constrained", dpi=120)
        ax = fig.add_subplot(aspect="equal")

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

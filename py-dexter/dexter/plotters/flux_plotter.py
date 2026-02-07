"""Plotter Parent class that provides plotting methods between the two fluxes."""

from dexter.types import Array1, EquilibriumType
import numpy as np
import matplotlib.pyplot as plt

from collections.abc import Callable

FIG_KW = {"figsize": (7, 6), "dpi": 120, "layout": "constrained"}
SUBPLOT_KW = {"xmargin": 0, "ymargin": 0}
PLOT_KW = {"c": "r", "linewidth": 1}
SCATTER_KW = {"c": "k", "s": 4, "zorder": 2}
PSI_WALL_BOUND = 0.5
PSIP_WALL_BOUND = 0.5


class _FluxPlotter:
    """Provides plotting functions between the two flux coordinates."""

    equilibrium_type: EquilibriumType
    psip_of_psi: Callable
    psi_of_psip: Callable

    # Not guaranteed to exist
    psi_wall: float | None
    psip_wall: float | None
    psi_array: Array1
    psip_array: Array1

    def plot_psip_of_psi(self, points: int = 1000, data: bool = False):
        r"""Plots $\psi(\psi_p)$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $\psi_p(\psi)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        psi_wall = getattr(self, "psi_wall", PSI_WALL_BOUND)
        psis = np.linspace(0, psi_wall, points)
        qs = np.asarray([self.psip_of_psi(psi) for psi in psis])

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$\psi_p(\psi)$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$\psi_p(\psi)$ $[Normalized]$")

        ax.plot(psis, qs, label=ax.get_ylabel(), **PLOT_KW)
        if self.equilibrium_type == "Numerical" and data:
            ax.scatter(
                self.psi_array,
                self.psip_array,
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.spines["bottom"].set_position("zero")
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_psi_of_psip(self, points: int = 1000, data: bool = False):
        r"""Plots $\psi(\psi_p)$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $\psi(\psi_p)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        psip_wall = getattr(self, "psip_wall", PSIP_WALL_BOUND)
        psips = np.linspace(0, psip_wall, points)
        qs = np.asarray([self.psi_of_psip(psip) for psip in psips])

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$\psi(\psi_p)$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$\psi(\psi_p)$ $[Normalized]$")

        ax.plot(psips, qs, label=ax.get_ylabel(), **PLOT_KW)
        if self.equilibrium_type == "Numerical" and data:
            ax.scatter(
                self.psip_array,
                self.psi_array,
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.spines["bottom"].set_position("zero")
        ax.grid(True)
        ax.legend()
        plt.show()

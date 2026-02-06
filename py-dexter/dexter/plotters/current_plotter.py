"""Current Parent class that provides plotting methods."""

from dexter.types import Array1, EquilibriumType
import numpy as np
import matplotlib.pyplot as plt

from collections.abc import Callable

FIG_KW = {"figsize": (7, 6), "dpi": 120, "layout": "constrained"}
SUBPLOT_KW = {"xmargin": 0, "ymargin": 0}
PLOT_KW = {"c": "r", "linewidth": 1}
OVERLAY_PLOT_KW = PLOT_KW | {"c": "b", "linestyle": "--"}
SCATTER_KW = {"c": "k", "s": 4, "zorder": 2}
SCATTER_KW = {"c": "k", "s": 4, "zorder": 2}
PSI_WALL_BOUND = 0.5
PSIP_WALL_BOUND = 0.5


class _CurrentPlotter:
    """Provides plotting functions for a Current's evaluation methods."""

    equilibrium_type: EquilibriumType
    g_of_psi: Callable
    g_of_psip: Callable
    i_of_psi: Callable
    i_of_psip: Callable
    dg_dpsi: Callable
    dg_dpsip: Callable
    di_dpsi: Callable
    di_dpsip: Callable

    # Not guaranteed to exist
    psi_wall: float | None
    psip_wall: float | None
    psi_array: Array1
    psip_array: Array1
    g_array: Array1
    i_array: Array1

    def plot_g_of_psi(self, points: int = 1000, data: bool = False):
        r"""Plots $g(\psi)$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $g(\psi)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        psi_wall = getattr(self, "psi_wall", PSI_WALL_BOUND)
        psis = np.linspace(0, psi_wall, points)
        gs = np.asarray([self.g_of_psi(psi) for psi in psis])

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$g(\psi)$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$g(\psi)$ $[Normalized]$")

        ax.plot(psis, gs, label=ax.get_ylabel(), **PLOT_KW)
        if self.equilibrium_type == "Numerical" and data:
            ax.scatter(
                self.psi_array, self.g_array, label=r"$data$ $points$", **SCATTER_KW
            )

        ax.spines["bottom"].set_position("zero")
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_g_of_psip(self, points: int = 1000, data: bool = False):
        r"""Plots $g(\psi_p)$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $g(\psi_p)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        psip_wall = getattr(self, "psip_wall", PSIP_WALL_BOUND)
        psips = np.linspace(0, psip_wall, points)
        gs = np.asarray([self.g_of_psip(psip) for psip in psips])

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$g(\psi_p)$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$g(\psi_p)$ $[Normalized]$")

        ax.plot(psips, gs, label=ax.get_ylabel(), **PLOT_KW)
        if self.equilibrium_type == "Numerical" and data:
            ax.scatter(
                self.psip_array, self.g_array, label=r"$data$ $points$", **SCATTER_KW
            )

        ax.spines["bottom"].set_position("zero")
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_i_of_psi(self, points: int = 1000, data: bool = False):
        r"""Plots $I(\psi)$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $I(\psi)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        psi_wall = getattr(self, "psi_wall", PSI_WALL_BOUND)
        psis = np.linspace(0, psi_wall, points)
        i = np.asarray([self.i_of_psi(psi) for psi in psis])

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$I(\psi)$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$I(\psi)$ $[Normalized]$")

        ax.plot(psis, i, label=ax.get_ylabel(), **PLOT_KW)
        if self.equilibrium_type == "Numerical" and data:
            ax.scatter(
                self.psi_array, self.i_array, label=r"$data$ $points$", **SCATTER_KW
            )

        ax.spines["bottom"].set_position("zero")
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_i_of_psip(self, points: int = 1000, data: bool = False):
        r"""Plots $I(\psi_p)$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $I(\psi_p)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        psip_wall = getattr(self, "psip_wall", PSIP_WALL_BOUND)
        psips = np.linspace(0, psip_wall, points)
        i = np.asarray([self.i_of_psip(psip) for psip in psips])

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$I(\psi_p)$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$I(\psi_p)$ $[Normalized]$")

        ax.plot(psips, i, label=ax.get_ylabel(), **PLOT_KW)
        if self.equilibrium_type == "Numerical" and data:
            ax.scatter(
                self.psip_array, self.i_array, label=r"$data$ $points$", **SCATTER_KW
            )

        ax.spines["bottom"].set_position("zero")
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_dg_dpsi(self, points: int = 1000):
        r"""Plots $dg(\psi)/d\psi$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $dg(\psi)/d\psi$. Defaults to 1000.
        """

        psi_wall = getattr(self, "psi_wall", PSI_WALL_BOUND)
        psis = np.linspace(0, psi_wall, points)
        ds = np.asarray([self.dg_dpsi(psi) for psi in psis])

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$dg(\psi)/\psi$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$dg(\psi)/d\psi$ $[Normalized]$")

        ax.plot(psis, ds, label=ax.get_ylabel(), **PLOT_KW)

        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_dg_dpsip(self, points: int = 1000):
        r"""Plots $dg(\psi_p)/d\psi_p$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $dg(\psi_p)/d\psi_p$. Defaults to 1000.
        """

        psip_wall = getattr(self, "psip_wall", PSIP_WALL_BOUND)
        psips = np.linspace(0, psip_wall, points)
        ds = np.asarray([self.dg_dpsip(psip) for psip in psips])

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$dg(\psi_p)/\psi_p$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$dg(\psi_p)/d\psi_p$ $[Normalized]$")

        ax.plot(psips, ds, label=ax.get_ylabel(), **PLOT_KW)

        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_di_dpsi(self, points: int = 1000):
        r"""Plots $dI(\psi)/d\psi$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $dI(\psi)/d\psi$. Defaults to 1000.
        """

        psi_wall = getattr(self, "psi_wall", PSI_WALL_BOUND)
        psis = np.linspace(0, psi_wall, points)
        ds = np.asarray([self.di_dpsi(psi) for psi in psis])

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$dI(\psi)/\psi$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$dI(\psi)/d\psi$ $[Normalized]$")

        ax.plot(psis, ds, label=ax.get_ylabel(), **PLOT_KW)

        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_di_dpsip(self, points: int = 1000):
        r"""Plots $dI(\psi_p)/d\psi_p$.


        Parameters
        ----------
        points
            The number of points in which to evaluate $dI(\psi_p)/d\psi_p$. Defaults to 1000.
        """

        psip_wall = getattr(self, "psip_wall", PSIP_WALL_BOUND)
        psips = np.linspace(0, psip_wall, points)
        ds = np.asarray([self.di_dpsip(psip) for psip in psips])

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$dI(\psi_p)/\psi_p$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$dI(\psi_p)/d\psi_p$ $[Normalized]$")

        ax.plot(psips, ds, label=ax.get_ylabel(), **PLOT_KW)

        ax.grid(True)
        ax.legend()
        plt.show()

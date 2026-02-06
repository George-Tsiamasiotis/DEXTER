"""Qfactor Parent class that provides plotting methods."""

from dexter.types import Array1, EquilibriumType
import numpy as np
import matplotlib.pyplot as plt

from collections.abc import Callable

FIG_KW = {"figsize": (7, 6), "dpi": 120, "layout": "constrained"}
SUBPLOT_KW = {"xmargin": 0, "ymargin": 0}
PLOT_KW = {"c": "r", "linewidth": 1}
OVERLAY_PLOT_KW = PLOT_KW | {"c": "b", "linestyle": "--"}
SCATTER_KW = {"c": "k", "s": 4, "zorder": 2}
PSI_WALL_BOUND = 0.5
PSIP_WALL_BOUND = 0.5


class _QfactorPlotter:
    """Provides plotting functions for a Qfactor's evaluation methods."""

    equilibrium_type: EquilibriumType
    q_of_psi: Callable
    q_of_psip: Callable
    psip_of_psi: Callable
    psi_of_psip: Callable
    dpsip_dpsi: Callable
    dpsi_dpsip: Callable
    iota_of_psi: Callable
    iota_of_psip: Callable

    # Not guaranteed to exist
    psi_wall: float | None
    psip_wall: float | None
    psi_array: Array1
    psip_array: Array1
    q_array: Array1

    def plot_q_of_psi(self, points: int = 1000, data: bool = False):
        r"""Plots $q(\psi)$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $q(\psi)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        psi_wall = getattr(self, "psi_wall", PSI_WALL_BOUND)
        psis = np.linspace(0, psi_wall, points)
        qs = np.asarray([self.q_of_psi(psi) for psi in psis])

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$q(\psi)$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$q(\psi)$")

        ax.plot(psis, qs, label=ax.get_ylabel(), **PLOT_KW)
        if self.equilibrium_type == "Numerical" and data:
            ax.scatter(
                self.psi_array,
                self.q_array,
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.spines["bottom"].set_position("zero")
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_q_of_psip(self, points: int = 1000, data: bool = False):
        r"""Plots $q(\psi_p)$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $q(\psi_p)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        psip_wall = getattr(self, "psip_wall", PSIP_WALL_BOUND)
        psips = np.linspace(0, psip_wall, points)
        qs = np.asarray([self.q_of_psip(psip) for psip in psips])

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$q(\psi_p)$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$q(\psi_p)$")

        ax.plot(psips, qs, label=ax.get_ylabel(), **PLOT_KW)
        if self.equilibrium_type == "Numerical" and data:
            ax.scatter(
                self.psip_array,
                self.q_array,
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.spines["bottom"].set_position("zero")
        ax.grid(True)
        ax.legend()
        plt.show()

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

    def plot_dpsip_dpsi(self, points: int = 1000):
        r"""Plots $d\psi_p(\psi)/d\psi$ and $\iota(\psi)$.

        This is a check to make sure the two quantities do indeed overlap.

        Parameters
        ----------
        points
            The number of points in which to evaluate the two splines. Defaults to 1000.
        """

        psi_wall = getattr(self, "psi_wall", PSI_WALL_BOUND)
        psis = np.linspace(0, psi_wall, points)
        ds = np.asarray([self.dpsip_dpsi(psi) for psi in psis])
        iotas = np.asarray([self.iota_of_psi(psi) for psi in psis])

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$d\psi_p(\psi)/d\psi, \quad \iota(\psi)$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$d\psi_p(\psi)/d\psi, \quad \iota(\psi)$")

        ax.plot(psis, ds, label=r"$d\psi_p(\psi)/d\psi$", **PLOT_KW)
        ax.plot(psis, iotas, label=r"$\iota(\psi)$", **OVERLAY_PLOT_KW)

        ax.spines["bottom"].set_position("zero")
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_dpsi_dpsip(self, points: int = 1000):
        r"""Plots $d\psi(\psi_p)/d\psi_p$ and $q(\psi_p)$.

        This is a check to make sure the two quantities do indeed overlap.

        Parameters
        ----------
        points
            The number of points in which to evaluate the two splines. Defaults to 1000.
        """

        psip_wall = getattr(self, "psip_wall", PSIP_WALL_BOUND)
        psips = np.linspace(0, psip_wall, points)
        ds = np.asarray([self.dpsi_dpsip(psip) for psip in psips])
        qs = np.asarray([self.q_of_psip(psip) for psip in psips])

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$d\psi(\psi_p)/d\psi_p, \quad q(\psi_p)$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$d\psi(\psi_p)/d\psi_p, \quad q(\psi_p)$")

        ax.plot(psips, ds, label=r"$d\psi(\psi_p)/d\psi_p$", **PLOT_KW)
        ax.plot(psips, qs, label=r"$q(\psi_p)$", **OVERLAY_PLOT_KW)

        ax.spines["bottom"].set_position("zero")
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_iota_of_psi(self, points: int = 1000):
        r"""Plots $\iota(\psi)$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $\iota(\psi)$. Defaults to 1000.
        """

        psi_wall = getattr(self, "psi_wall", PSI_WALL_BOUND)
        psis = np.linspace(0, psi_wall, points)
        iotas = np.asarray([self.iota_of_psi(psi) for psi in psis])

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$\iota(\psi)$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$\iota(\psi)$")

        ax.plot(psis, iotas, label=ax.get_ylabel(), **PLOT_KW)

        ax.spines["bottom"].set_position("zero")
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_iota_of_psip(self, points: int = 1000):
        r"""Plots $\iota(\psi_p)$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $\iota(\psi_p)$. Defaults to 1000.
        """

        psip_wall = getattr(self, "psip_wall", PSIP_WALL_BOUND)
        psips = np.linspace(0, psip_wall, points)
        iotas = np.asarray([self.iota_of_psip(psip) for psip in psips])

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$\iota(\psi_p)$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$\iota(\psi_p)$")

        ax.plot(psips, iotas, label=ax.get_ylabel(), **PLOT_KW)

        ax.spines["bottom"].set_position("zero")
        ax.grid(True)
        ax.legend()
        plt.show()

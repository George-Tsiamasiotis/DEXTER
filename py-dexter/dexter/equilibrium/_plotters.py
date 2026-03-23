"""Plotter Parent classes that provide simple plotting methods

Method are inherited by the final rust wrapped objects.

These plotting methods are restricted to whatever resources each wrapped type has and are meant
for quick and simple diagnostic plots. They should not count to external objects to work.
"""

from dexter._core import _PyLarGeometry, _PyNcGeometry
from dexter._core import _PyUnityQfactor, _PyParabolicQfactor, _PyNcQfactor
from dexter._core import _PyLarCurrent, _PyNcCurrent
from dexter._core import _PyCosHarmonic, _PyNcHarmonic
import numpy as np
import matplotlib.pyplot as plt

from math import sqrt
from collections.abc import Callable

FIG_KW = {"figsize": (7, 6), "dpi": 120, "layout": "constrained"}
SUBPLOT_KW = {"xmargin": 0, "ymargin": 0}
PLOT_KW = {"c": "r", "linewidth": 1}
OVERLAY_PLOT_KW = PLOT_KW | {"c": "b", "linestyle": "--"}
SCATTER_KW = {"c": "k", "s": 4, "zorder": 2}
XAXIS_KW = {"c": "k", "linewidth": 1.5}
FLUX_SURFACE_KW = {"c": "b", "zorder": 2}
JACOBIAN_KW = {"levels": None, "cmap": "plasma", "zorder": 2}
LCFS_KW = {"c": "k", "linewidth": 2, "linestyle": "-", "zorder": 2}
RESONANCE_KW = {"c": "g", "linewidth": 1.5, "linestyle": ":", "zorder": 2}
PSI_LAST_BOUND = 0.5
PSIP_LAST_BOUND = 0.5


class _FluxPlotter:
    """Provides plotting functions between the two flux coordinates."""

    psip_of_psi: Callable
    psi_of_psip: Callable

    _rust: _PyNcGeometry | _PyUnityQfactor | _PyParabolicQfactor | _PyNcQfactor

    def plot_psip_of_psi(self, points: int = 1000, data: bool = False):
        r"""Plots $\psi(\psi_p)$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $\psi_p(\psi)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        psi_last = getattr(self, "psi_last", PSI_LAST_BOUND)
        psis = np.linspace(0, psi_last, points)
        qs = self.psip_of_psi(psis)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$\psi_p(\psi)$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$\psi_p(\psi)$ $[Normalized]$")

        ax.plot(psis, qs, label=ax.get_ylabel(), **PLOT_KW)
        if self._rust.equilibrium_type == "Numerical" and data:
            ax.scatter(
                getattr(self._rust, "psi_array"),
                getattr(self._rust, "psip_array"),
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.axhline(y=0, **XAXIS_KW)
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

        psip_last = getattr(self, "psip_last", PSIP_LAST_BOUND)
        psips = np.linspace(0, psip_last, points)
        qs = self.psi_of_psip(psips)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$\psi(\psi_p)$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$\psi(\psi_p)$ $[Normalized]$")

        ax.plot(psips, qs, label=ax.get_ylabel(), **PLOT_KW)
        if self._rust.equilibrium_type == "Numerical" and data:
            ax.scatter(
                getattr(self._rust, "psip_array"),
                getattr(self._rust, "psi_array"),
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.axhline(y=0, **XAXIS_KW)
        ax.grid(True)
        ax.legend()
        plt.show()


class _GeometryPlotter:
    """Provides plotting functions for a Geometry's evaluation methods."""

    r_of_psi: Callable
    r_of_psip: Callable
    psi_of_r: Callable
    psip_of_r: Callable

    _rust: _PyLarGeometry | _PyNcGeometry

    def plot_r_of_psi(self, points: int = 1000, data: bool = False):
        r"""Plots $r(\psi)$, where $r$ is in $[m]$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $r(\psi)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        psi_last = getattr(self, "psi_last")
        psis = np.linspace(0, psi_last, points)
        rs = self.r_of_psi(psis)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$r(\psi)$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$r(\psi)$ $[Normalized]$")

        ax.plot(psis, rs, label=ax.get_ylabel(), **PLOT_KW)
        if self._rust.equilibrium_type == "Numerical" and data:
            ax.scatter(
                getattr(self._rust, "psi_array"),
                getattr(self._rust, "r_array"),
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.axhline(y=0, **XAXIS_KW)
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_r_of_psip(self, points: int = 1000, data: bool = False):
        r"""Plots $r(\psi_p)$, where $r$ is in $[m]$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $r(\psi_p)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        psip_last = getattr(self, "psip_last")
        psips = np.linspace(0, psip_last, points)
        rs = self.r_of_psip(psips)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$r(\psi_p)$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$r(\psi_p)$ $[Normalized]$")

        ax.plot(psips, rs, label=ax.get_ylabel(), **PLOT_KW)
        if self._rust.equilibrium_type == "Numerical" and data:
            ax.scatter(
                getattr(self._rust, "psip_array"),
                getattr(self._rust, "r_array"),
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.axhline(y=0, **XAXIS_KW)
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_psi_of_r(self, points: int = 1000, data: bool = False):
        r"""Plots $\psi(r)$, where $r$ is in $[m]$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $\psi(r)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        rlast = getattr(self, "rlast")
        rs = np.linspace(0, rlast, points)
        psis = self.psi_of_r(rs)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$\psi(r)$ $profile$")
        ax.set_xlabel(r"$r$ $[m]$")
        ax.set_ylabel(r"$\psi(r)$ $[Normalized]$")

        ax.plot(rs, psis, label=ax.get_ylabel(), **PLOT_KW)
        if self._rust.equilibrium_type == "Numerical" and data:
            ax.scatter(
                getattr(self._rust, "r_array"),
                getattr(self._rust, "psi_array"),
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.axhline(y=0, **XAXIS_KW)
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_psip_of_r(self, points: int = 1000, data: bool = False):
        r"""Plots $\psi_p(r)$, where $r$ is in $[m]$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $\psi_p(r)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        rlast = self._rust.rlast
        rs = np.linspace(0, rlast, points)
        psips = self.psip_of_r(rs)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$\psi_p(r)$ $profile$")
        ax.set_xlabel(r"$r$ $[m]$")
        ax.set_ylabel(r"$\psi_p(r)$ $[Normalized]$")

        ax.plot(rs, psips, label=ax.get_ylabel(), **PLOT_KW)
        if self._rust.equilibrium_type == "Numerical" and data:
            ax.scatter(
                getattr(self._rust, "r_array"),
                getattr(self._rust, "psip_array"),
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.axhline(y=0, **XAXIS_KW)
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_last(self):
        r"""Plots the device's Last Closed Flux Surface (LCFS) in the $R, Z$ frame."""

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW | {"aspect": "equal"})
        ax.set_title(r"$LCFS$")
        ax.set_xlabel(r"$R[m]$")
        ax.set_ylabel(r"$Z[m]$")

        zlab_last = self._rust.zlab_last
        rlab_last = self._rust.rlab_last

        ax.plot(rlab_last, zlab_last, **LCFS_KW)

        ax.grid(True)
        plt.show()


class _NumericalGeometryPlotter:
    """Provides plotting functions for a Numerical Geometry's evaluation methods."""

    _rust: _PyNcGeometry

    def plot_flux_surfaces(self, number: int = 20):
        r"""Plots the flux surfaces in the $R-Z$ frame.

        Parameters
        ----------
        number
            The number of flux surfaces to (try to) plot. Defaults to 20.
        """

        rlab_array = self._rust.rlab_array
        zlab_array = self._rust.zlab_array

        step = max([1, int(self._rust.shape[0] / number)])
        print(f"Displaying {int(self._rust.shape[0] / step)} surfaces.")

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW | {"aspect": "equal"})
        ax.set_title(r"$Flux$ $surfaces$")
        ax.set_xlabel(r"$R[m]$")
        ax.set_ylabel(r"$Z[m]$")

        for i in range(0, self._rust.shape[0], step):
            ax.plot(rlab_array[i], zlab_array[i], **FLUX_SURFACE_KW)

        geom_center = (self._rust.rgeo, self._rust.zaxis)
        axis = (self._rust.raxis, self._rust.zaxis)

        ax.plot(*geom_center, "ko", markersize=4, label="$R_{axis}$")
        ax.plot(*axis, "ro", markersize=4, label="$R_{geo}$")

        def format_coord(x, y):
            r = sqrt((axis[0] - x) ** 2 + (axis[1] - y) ** 2)
            return f"(R, Z) = ({x:.5g}, {y:.5g}), r={r:.5g}[m]"

        setattr(ax, "format_coord", format_coord)

        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_jacobian(self, levels: int = 20):
        r"""Plots the Jacobian $J(R, Z)$

        Parameters
        ----------
        levels
            The number of contour levels. Defaults to 20.
        """

        rlab_array = self._rust.rlab_array
        zlab_array = self._rust.zlab_array
        jacobian_array = self._rust.jacobian_array

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW | {"aspect": "equal"})
        ax.set_title(r"$Jacobian$")
        ax.set_xlabel(r"$R[m]$")
        ax.set_ylabel(r"$Z[m]$")

        c = ax.contourf(
            rlab_array,
            zlab_array,
            jacobian_array,
            **JACOBIAN_KW | {"levels": levels},
        )
        plt.colorbar(c, ax=ax)

        ax.plot(rlab_array[-1], zlab_array[-1], **LCFS_KW)

        ax.grid(True)
        plt.show()


class _QfactorPlotter:
    """Provides plotting functions for a Qfactor's evaluation methods."""

    q_of_psi: Callable
    q_of_psip: Callable
    dpsip_dpsi: Callable
    dpsi_dpsip: Callable
    iota_of_psi: Callable
    iota_of_psip: Callable

    _rust: _PyUnityQfactor | _PyParabolicQfactor | _PyNcQfactor

    def plot_q_of_psi(self, points: int = 1000, data: bool = False):
        r"""Plots $q(\psi)$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $q(\psi)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        psi_last = getattr(self._rust, "psi_last", PSI_LAST_BOUND)
        psis = np.linspace(0, psi_last, points)
        qs = self.q_of_psi(psis)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$q(\psi)$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$q(\psi)$")

        ax.plot(psis, qs, label=ax.get_ylabel(), **PLOT_KW)
        if self._rust.equilibrium_type == "Numerical" and data:
            ax.scatter(
                getattr(self._rust, "psi_array"),
                getattr(self._rust, "q_array"),
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.axhline(y=0, **XAXIS_KW)
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

        psip_last = getattr(self, "psip_last", PSIP_LAST_BOUND)
        psips = np.linspace(0, psip_last, points)
        qs = self.q_of_psip(psips)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$q(\psi_p)$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$q(\psi_p)$")

        ax.plot(psips, qs, label=ax.get_ylabel(), **PLOT_KW)
        if self._rust.equilibrium_type == "Numerical" and data:
            ax.scatter(
                getattr(self._rust, "psip_array"),
                getattr(self._rust, "q_array"),
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.axhline(y=0, **XAXIS_KW)
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

        psi_last = getattr(self, "psi_last", PSI_LAST_BOUND)
        psis = np.linspace(0, psi_last, points)
        ds = self.dpsip_dpsi(psis)
        iotas = self.iota_of_psi(psis)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$d\psi_p(\psi)/d\psi, \quad \iota(\psi)$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$d\psi_p(\psi)/d\psi, \quad \iota(\psi)$")

        ax.plot(psis, ds, label=r"$d\psi_p(\psi)/d\psi$", **PLOT_KW)
        ax.plot(psis, iotas, label=r"$\iota(\psi)$", **OVERLAY_PLOT_KW)

        ax.axhline(y=0, **XAXIS_KW)
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

        psip_last = getattr(self, "psip_last", PSIP_LAST_BOUND)
        psips = np.linspace(0, psip_last, points)
        ds = self.dpsi_dpsip(psips)
        qs = self.q_of_psip(psips)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$d\psi(\psi_p)/d\psi_p, \quad q(\psi_p)$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$d\psi(\psi_p)/d\psi_p, \quad q(\psi_p)$")

        ax.plot(psips, ds, label=r"$d\psi(\psi_p)/d\psi_p$", **PLOT_KW)
        ax.plot(psips, qs, label=r"$q(\psi_p)$", **OVERLAY_PLOT_KW)

        ax.axhline(y=0, **XAXIS_KW)
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

        psi_last = getattr(self, "psi_last", PSI_LAST_BOUND)
        psis = np.linspace(0, psi_last, points)
        iotas = self.iota_of_psi(psis)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$\iota(\psi)$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$\iota(\psi)$")

        ax.plot(psis, iotas, label=ax.get_ylabel(), **PLOT_KW)

        ax.axhline(y=0, **XAXIS_KW)
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

        psip_last = getattr(self, "psip_last", PSIP_LAST_BOUND)
        psips = np.linspace(0, psip_last, points)
        iotas = self.iota_of_psip(psips)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$\iota(\psi_p)$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$\iota(\psi_p)$")

        ax.plot(psips, iotas, label=ax.get_ylabel(), **PLOT_KW)

        ax.axhline(y=0, **XAXIS_KW)
        ax.grid(True)
        ax.legend()
        plt.show()


class _CurrentPlotter:
    """Provides plotting functions for a Current's evaluation methods."""

    g_of_psi: Callable
    g_of_psip: Callable
    i_of_psi: Callable
    i_of_psip: Callable
    dg_dpsi: Callable
    dg_dpsip: Callable
    di_dpsi: Callable
    di_dpsip: Callable

    _rust: _PyLarCurrent | _PyNcCurrent

    def plot_g_of_psi(self, points: int = 1000, data: bool = False):
        r"""Plots $g(\psi)$.

        Parameters
        ----------
        points
            The number of points in which to evaluate $g(\psi)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        psi_last = getattr(self, "psi_last", PSI_LAST_BOUND)
        psis = np.linspace(0, psi_last, points)
        gs = self.g_of_psi(psis)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$g(\psi)$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$g(\psi)$ $[Normalized]$")

        ax.plot(psis, gs, label=ax.get_ylabel(), **PLOT_KW)
        if self._rust.equilibrium_type == "Numerical" and data:
            ax.scatter(
                getattr(self._rust, "psi_array"),
                getattr(self._rust, "g_array"),
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.axhline(y=0, **XAXIS_KW)
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

        psip_last = getattr(self, "psip_last", PSIP_LAST_BOUND)
        psips = np.linspace(0, psip_last, points)
        gs = self.g_of_psip(psips)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$g(\psi_p)$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$g(\psi_p)$ $[Normalized]$")

        ax.plot(psips, gs, label=ax.get_ylabel(), **PLOT_KW)
        if self._rust.equilibrium_type == "Numerical" and data:
            ax.scatter(
                getattr(self._rust, "psip_array"),
                getattr(self._rust, "g_array"),
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.axhline(y=0, **XAXIS_KW)
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

        psi_last = getattr(self, "psi_last", PSI_LAST_BOUND)
        psis = np.linspace(0, psi_last, points)
        i = self.i_of_psi(psis)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$I(\psi)$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$I(\psi)$ $[Normalized]$")

        ax.plot(psis, i, label=ax.get_ylabel(), **PLOT_KW)
        if self._rust.equilibrium_type == "Numerical" and data:
            ax.scatter(
                getattr(self._rust, "psi_array"),
                getattr(self._rust, "i_array"),
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.axhline(y=0, **XAXIS_KW)
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

        psip_last = getattr(self, "psip_last", PSIP_LAST_BOUND)
        psips = np.linspace(0, psip_last, points)
        i = self.i_of_psip(psips)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$I(\psi_p)$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$I(\psi_p)$ $[Normalized]$")

        ax.plot(psips, i, label=ax.get_ylabel(), **PLOT_KW)
        if self._rust.equilibrium_type == "Numerical" and data:
            ax.scatter(
                getattr(self._rust, "psip_array"),
                getattr(self._rust, "i_array"),
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.axhline(y=0, **XAXIS_KW)
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

        psi_last = getattr(self, "psi_last", PSI_LAST_BOUND)
        psis = np.linspace(0, psi_last, points)
        ds = self.dg_dpsi(psis)

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

        psip_last = getattr(self, "psip_last", PSIP_LAST_BOUND)
        psips = np.linspace(0, psip_last, points)
        ds = self.dg_dpsip(psips)

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

        psi_last = getattr(self, "psi_last", PSI_LAST_BOUND)
        psis = np.linspace(0, psi_last, points)
        ds = self.di_dpsi(psis)

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

        psip_last = getattr(self, "psip_last", PSIP_LAST_BOUND)
        psips = np.linspace(0, psip_last, points)
        ds = self.di_dpsip(psips)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW)
        ax.set_title(r"$dI(\psi_p)/\psi_p$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$dI(\psi_p)/d\psi_p$ $[Normalized]$")

        ax.plot(psips, ds, label=ax.get_ylabel(), **PLOT_KW)

        ax.grid(True)
        ax.legend()
        plt.show()


class _HarmonicPlotter:
    """Provides plotting functions for a Harmonic's evaluation methods."""

    alpha_of_psi: Callable
    alpha_of_psip: Callable
    phase_of_psi: Callable
    phase_of_psip: Callable

    _rust: _PyCosHarmonic | _PyNcHarmonic

    def plot_alpha_of_psi(self, points: int = 1000, data: bool = False):
        r"""Plots the harmonic's amplitude $\alpha(\psi)$.

        Note
        ----

        It is assumed that the amplitude $\alpha$ is only a function of the flux.

        Parameters
        ----------
        points
            The number of points in which to evaluate $\alpha(\psi)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        psi_last = getattr(self, "psi_last", PSI_LAST_BOUND)
        psis = np.linspace(0, psi_last, points)
        alphas = self.alpha_of_psi(psis, 0, 0, 0)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW | {"ymargin": 0.1})
        ax.set_title(r"$\alpha(\psi)$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$\alpha(\psi)$ $[Normalized]$")

        ax.plot(psis, alphas, label=ax.get_ylabel(), **PLOT_KW)
        if self._rust.equilibrium_type == "Numerical" and data:
            ax.scatter(
                getattr(self._rust, "psi_array"),
                getattr(self._rust, "alpha_array"),
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.axhline(y=0, **XAXIS_KW)
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_alpha_of_psip(self, points: int = 1000, data: bool = False):
        r"""Plots the harmonic's amplitude $\alpha(\psi_p)$.

        Note
        ----

        It is assumed that the amplitude $\alpha$ is only a function of the flux.

        Parameters
        ----------
        points
            The number of points in which to evaluate $\alpha(\psi_p)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        """

        psip_last = getattr(self, "psip_last", PSIP_LAST_BOUND)
        psips = np.linspace(0, psip_last, points)
        alphas = self.alpha_of_psip(psips, 0, 0, 0)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW | {"ymargin": 0.1})
        ax.set_title(r"$\alpha(\psi_p)$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$\alpha(\psi_p)$ $[Normalized]$")

        ax.plot(psips, alphas, label=ax.get_ylabel(), **PLOT_KW)
        if self._rust.equilibrium_type == "Numerical" and data:
            ax.scatter(
                getattr(self._rust, "psip_array"),
                getattr(self._rust, "alpha_array"),
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        ax.axhline(y=0, **XAXIS_KW)
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_phase_of_psi(
        self, points: int = 1000, data: bool = False, resonance: bool = True
    ):
        r"""Plots the harmonic's phase $\phi(\psi)$.

        Note
        ----

        It is assumed that the amplitude $\phi$ is only a function of the flux.

        Parameters
        ----------
        points
            The number of points in which to evaluate $\phi(\psi)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        resonance
            Whether or not to plot the resonance location, if `phase_method` is `Resonance`.
            Defaults to True.
        """

        psi_last = getattr(self, "psi_last", PSI_LAST_BOUND)
        psis = np.linspace(0, psi_last, points)
        phis = self.phase_of_psi(psis, 0, 0, 0)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW | {"ymargin": 0.1})
        ax.set_title(r"$\phi(\psi)$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$\phi(\psi)$ $[Normalized]$")

        ax.plot(psis, phis, label=ax.get_ylabel(), **PLOT_KW)
        if self._rust.equilibrium_type == "Numerical" and data:
            ax.scatter(
                getattr(self._rust, "psi_array"),
                getattr(self._rust, "phase_array"),
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        if resonance and (getattr(self, "phase_method", False) == "Resonance"):
            psi_res = getattr(self, "psi_phase_resonance")
            ax.axvline(
                x=psi_res,
                label=f"$Resonance$ $(n/m = {self._rust.n}/{self._rust.m})$",
                **RESONANCE_KW,
            )

        ax.axhline(y=0, **XAXIS_KW)
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_phase_of_psip(
        self, points: int = 1000, data: bool = False, resonance: bool = True
    ):
        r"""Plots the harmonic's phase $\phi(\psi_p)$.

        Note
        ----

        It is assumed that the amplitude $\phi$ is only a function of the flux.

        Parameters
        ----------
        points
            The number of points in which to evaluate $\phi(\psi_p)$. Defaults to 1000.
        data
            Whether or not to plot the data array points (numerical equilibria only). Defaults to False.
        resonance
            Whether or not to plot the resonance location, if `phase_method` is `Resonance`.
            Defaults to True.
        """

        psip_last = getattr(self, "psip_last", PSIP_LAST_BOUND)
        psips = np.linspace(0, psip_last, points)
        phis = self.phase_of_psip(psips, 0, 0, 0)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW | {"ymargin": 0.1})
        ax.set_title(r"$\phi(\psi_p)$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$\phi(\psi_p)$ $[Normalized]$")

        ax.plot(psips, phis, label=ax.get_ylabel(), **PLOT_KW)
        if self._rust.equilibrium_type == "Numerical" and data:
            ax.scatter(
                getattr(self._rust, "psip_array"),
                getattr(self._rust, "phase_array"),
                label=r"$data$ $points$",
                **SCATTER_KW,
            )

        if resonance and (getattr(self, "phase_method", False) == "Resonance"):
            psip_res = getattr(self, "psip_phase_resonance")
            ax.axvline(
                x=psip_res,
                label=f"$Resonance$ $(n/m = {self._rust.n}/{self._rust.m})$",
                **RESONANCE_KW,
            )

        ax.axhline(y=0, **XAXIS_KW)
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_dalpha_of_psi(self, points: int = 1000):
        r"""Plots the harmonic's amplitude's deritave $d\alpha(\psi)/d\psi$.

        Note
        ----

        It is assumed that the amplitude $\alpha$ is only a function of the flux.

        Parameters
        ----------
        points
            The number of points in which to evaluate $d\alpha(\psi)/d\psi$. Defaults to 1000.
        """

        psi_last = getattr(self, "psi_last", PSI_LAST_BOUND)
        psis = np.linspace(0, psi_last, points)
        alphas = self.alpha_of_psi(psis, 0, 0, 0)
        dalphas = np.gradient(alphas)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW | {"ymargin": 0.1})
        ax.set_title(r"$d\alpha(\psi)/d\psi$ $profile$")
        ax.set_xlabel(r"$\psi$ $[Normalized]$")
        ax.set_ylabel(r"$d\alpha(\psi)/d\psi$ $[Normalized]$")

        ax.plot(psis, dalphas, label=ax.get_ylabel(), **PLOT_KW)

        ax.axhline(y=0, **XAXIS_KW)
        ax.grid(True)
        ax.legend()
        plt.show()

    def plot_dalpha_of_psip(self, points: int = 1000):
        r"""Plots the harmonic's amplitude's deritave $d\alpha(\psi_p)/d\psi_p$.

        Note
        ----

        It is assumed that the amplitude $\alpha$ is only a function of the flux.

        Parameters
        ----------
        points
            The number of points in which to evaluate $d\alpha(\psi_p)/d\psi_p$. Defaults to 1000.
        """

        psip_last = getattr(self, "psip_last", PSIP_LAST_BOUND)
        psips = np.linspace(0, psip_last, points)
        alphas = self.alpha_of_psip(psips, 0, 0, 0)
        dalphas = np.gradient(alphas)

        fig = plt.figure(**FIG_KW)
        ax = fig.add_subplot(**SUBPLOT_KW | {"ymargin": 0.1})
        ax.set_title(r"$d\alpha(\psi_p)/d\psi_p$ $profile$")
        ax.set_xlabel(r"$\psi_p$ $[Normalized]$")
        ax.set_ylabel(r"$d\alpha(\psi_p)/d\psi_p$ $[Normalized]$")

        ax.plot(psips, dalphas, label=ax.get_ylabel(), **PLOT_KW)

        ax.axhline(y=0, **XAXIS_KW)
        ax.grid(True)
        ax.legend()
        plt.show()

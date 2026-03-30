"""Plotter Parent classes that provide simple plotting methods"""

import numpy as np
import matplotlib.pyplot as plt
from math import floor, log10, sqrt
from math import pi as PI
from warnings import warn
from cycler import cycler
from matplotlib.patches import Patch
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from typing import Any

from dexter._core import _PyParticle, _PyQueue
from dexter import Equilibrium, LarGeometry, OrbitType
from dexter.types import Array1

dpi = 120
figsize = (10, 7)
target_points = 50_000
SCATTER_KW = {"s": 0.8, "c": "blue"}
LABEL_KW = {"labelpad": 10, "rotation": 0, "fontsize": 15}
CARTESIAN_POINCARE_FIG_KW = {"figsize": (9, 6), "layout": "constrained", "dpi": 120}
CARTESIAN_POINCARE_SCATTER_KW = {"s": 0.3, "marker": "o"}
CARTESIAN_POINCARE_INITIAL_KW = {
    "c": "k",
    "marker": "x",
    "markersize": 3,
    "alpha": 0.4,
    "zorder": -2,
}
RZ_POINCARE_PLOT_KW = {
    "marker": ".",
    "markersize": 0.8,
    "color": "blue",
    "alpha": 0.6,
    "linestyle": "",
}
RZ_POINCARE_INITIAL_KW = {
    "c": "k",
    "marker": "x",
    "markersize": 3,
    "alpha": 0.4,
    "zorder": -2,
}
QKINETIC_POINCARE_FIG_KW = {"figsize": (9, 6), "layout": "constrained", "dpi": 120}
COPASSING_COLOR = "xkcd:light purple"
CUPASSING_COLOR = "xkcd:navy blue"
TRAPPED_COLOR = "xkcd:blue"
UNDEFINED_COLOR = "xkcd:coral"
UNCLASSIFIED_COLOR = "xkcd:crimson"


class _ParticlePlotter:
    """Provides plotting functions for the Particle Class."""

    fig: Figure
    ax: Axes
    axes: Any  # list[Axes]

    _rust: _PyParticle

    def plot_evolution(
        self,
        percentage: float = 100,
        downsample: bool = True,
        show: bool = True,
    ):
        """Plots the time evolution of an integrated Particle.

        Parameters
        ----------
        percentage: float
            The percentage of the evolution to plot. Defaults to 100.
        downsample: bool
            Whether or not to downsample the evolution arrays. This can be
            really helpful when plotting arrays with a lot of points, since it
            drastically improves both figure creation and interaction
            performance. Downsampling is done by increasing the [::step] just
            enough so that the final number of points is more than 50.000.
            Defaults to True
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """

        if self._rust.steps_stored == 0:
            raise Exception("Particle hasn't been integrated")

        if percentage < 0 or percentage > 100:
            raise ValueError("Percentage must be between 0 and 100.")

        step = 1
        if downsample:
            length = len(self._rust.t_array)
            oom = floor(log10(length))  # order of magnitude of number of points
            target_oom = floor(log10(target_points))
            if oom > target_oom:
                step = 10 ** (oom - target_oom)

        t_array = self._rust.t_array
        points = min(
            int(np.floor(len(t_array) * percentage / 100)),
            len(t_array),
        )

        t = t_array[:points][::step]
        psi = self._rust.psi_array[:points][::step]
        psip = self._rust.psip_array[:points][::step]
        theta = self._rust.theta_array[:points][::step]
        zeta = self._rust.zeta_array[:points][::step]
        rho = self._rust.rho_array[:points][::step]
        # mu = self.mu_array[:points][::step]
        ptheta = self._rust.ptheta_array[:points][::step]
        pzeta = self._rust.pzeta_array[:points][::step]
        energy = self._rust.energy_array[:points][::step]

        if downsample and len(t) > target_points * 10:
            warn("Downsampling did not work..")

        # ===========================

        self.fig = plt.figure(figsize=figsize, layout="constrained", dpi=dpi)

        self.axes = self.fig.subplots(8, 1, sharex=True)
        for ax in self.axes:
            ax.yaxis.set_ticks_position("right")
            ax.yaxis.set_label_position("left")

        self.axes[0].scatter(t, psi, **SCATTER_KW)
        self.axes[1].scatter(t, psip, **SCATTER_KW)
        self.axes[2].scatter(t, theta, **SCATTER_KW)
        self.axes[3].scatter(t, zeta, **SCATTER_KW)
        self.axes[4].scatter(t, rho, **SCATTER_KW)
        self.axes[5].scatter(t, ptheta, **SCATTER_KW)
        self.axes[6].scatter(t, pzeta, **SCATTER_KW)
        self.axes[7].scatter(t, energy, **SCATTER_KW)

        self.axes[0].set_ylabel(r"$\psi$", **LABEL_KW)
        self.axes[1].set_ylabel(r"$\psi_p$", **LABEL_KW)
        self.axes[2].set_ylabel(r"$\theta$", **LABEL_KW)
        self.axes[3].set_ylabel(r"$\zeta$", **LABEL_KW)
        self.axes[4].set_ylabel(r"$\rho_{||}$", **LABEL_KW)
        self.axes[5].set_ylabel(r"$P_\theta$", **LABEL_KW)
        self.axes[6].set_ylabel(r"$P_\zeta$", **LABEL_KW)
        self.axes[7].set_ylabel(r"$E$", **LABEL_KW)
        self.axes[7].set_xlabel(r"$t[Normalized]$", **LABEL_KW)

        # Zoom out Pzeta plot
        if abs(np.std(pzeta)) < 1e-6:
            ylim = np.array(self.axes[6].get_ylim())
            self.axes[6].set_ylim(np.sort([ylim[0] / 3, ylim[1] * 3]).tolist())
        # Zoom out Energy plot
        if abs(np.std(energy)) < 1e-6:
            ylim = np.array(self.axes[7].get_ylim())
            self.axes[7].set_ylim(np.sort([ylim[0] / 2, ylim[1] * 2]).tolist())

        if show:
            plt.show()
            plt.close()

    def plot_db_drift(
        self,
        equilibrium: Equilibrium,
        show: bool = True,
    ):
        r"""Plots the magnetic field strength $B(\psi/\psi_p, \theta)$ along the particle's orbit.

        Parameters
        ----------
        equilibrium
            The equilibrium in which the particle was integrated.
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """

        if self._rust.steps_stored == 0:
            raise Exception("Particle hasn't been integrated")

        t = self._rust.t_array
        psi = self._rust.psi_array
        psip = self._rust.psip_array
        theta = self._rust.theta_array

        bfield = equilibrium.bfield
        coordinate = self._rust.initial_conditions.flux0.kind
        if coordinate == "toroidal":
            b = bfield.b_of_psi(psi, theta % (2 * PI))
        else:
            b = bfield.b_of_psip(psip, theta % (2 * PI))

        # ===========================

        self.fig = plt.figure(figsize=figsize, layout="constrained", dpi=dpi)
        self.fig.suptitle(r"$\Delta B(\psi/\psi_p, \theta)\ drift\ during\ transit$")

        self.axes = self.fig.subplots(4, 1, sharex=True)
        for ax in self.axes:
            ax.yaxis.set_ticks_position("right")
            ax.yaxis.set_label_position("left")

        self.axes[0].scatter(t, b, **SCATTER_KW)
        self.axes[1].scatter(t, psi, **SCATTER_KW)
        self.axes[2].scatter(t, psip, **SCATTER_KW)
        self.axes[3].scatter(t, theta, **SCATTER_KW)

        self.axes[0].set_ylabel(r"$B$", **LABEL_KW)
        self.axes[1].set_ylabel(r"$\psi$", **LABEL_KW)
        self.axes[2].set_ylabel(r"$\psi_p$", **LABEL_KW)
        self.axes[3].set_ylabel(r"$\theta$", **LABEL_KW)
        self.axes[3].set_xlabel(r"$t[Normalized]$", **LABEL_KW)

        if show:
            plt.show()
            plt.close()

    def plot_poloidal_drift(
        self,
        percentage: float = 100,
        show: bool = True,
    ):
        r"""Plots the $\theta-\psi$ drift of the particle.

        Note
        ----

        This method simply plots the $\theta$ and $\psi$ time arrays in a polar projection. It is not aware
        of the geometry of the device so the plot limit is not the actual $\psi_{LCFS}$.

        Parameters
        ----------
        percentage: float
            The percentage of the evolution to plot. Defaults to 100.
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """

        if self._rust.steps_stored == 0:
            raise Exception("Particle hasn't been integrated")

        if percentage < 0 or percentage > 100:
            raise ValueError("Percentage must be between 0 and 100.")

        # ===========================

        self.fig = plt.figure(figsize=figsize, layout="constrained", dpi=dpi)
        self.ax = self.fig.add_subplot(polar=True)

        self.ax.scatter(self._rust.theta_array, self._rust.psi_array, **SCATTER_KW)
        self.ax.set_title(r"$\theta-\psi$ drift")

        self.ax.set_rorigin(0)  # type: ignore
        self.ax.set_rlabel_position(22.5)  # type: ignore

        if show:
            plt.show()
            plt.close()


class _QueuePlotter:
    """Provides plotting functions for the Queue Class."""

    fig: Figure
    ax: Axes
    axes: Any  # list[Axes]

    _rust: _PyQueue

    def plot_energies(self, show=True):
        """Plots a stem plot of all the particles' energies.

        Particles are visited in order of instantiation.

        Useful to visualize the energy spans of the particles under study.

        Parameters
        ----------
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """

        self.fig = plt.figure(figsize=(9, 7), layout="constrained", dpi=120)
        self.ax = self.fig.add_subplot()

        self.ax.stem(self._rust.energy_array, linefmt="k--", markerfmt="blue")
        self.ax.set_title(r"$Particle\ Energies$")
        self.ax.set_xlabel(r"$Particle\ \#$")
        self.ax.set_ylabel(r"$Energy\ [Normalized]$")
        self.ax.set_ybound(lower=0)

        if show:
            plt.show()

    def plot_steps_taken(self, show=True):
        """Plots a stem plot of all the number of steps each particle has taken.

        Particles are visited in order of instantiation.

        Parameters
        ----------
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """

        self.fig = plt.figure(figsize=(9, 7), layout="constrained", dpi=120)
        self.ax = self.fig.add_subplot()

        self.ax.stem(self._rust.steps_taken_array, linefmt="k--", markerfmt="blue")
        self.ax.set_title(r"$Particle\ steps\ taken$")
        self.ax.set_xlabel(r"$Particle\ \#$")
        self.ax.set_ylabel(r"$Steps\ taken$")
        self.ax.set_ybound(lower=0)

        if show:
            plt.show()

    def plot_steps_stored(self, show=True):
        """Plots a stem plot of all the number of steps each particle has stored.

        Particles are visited in order of instantiation.

        Parameters
        ----------
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """

        self.fig = plt.figure(figsize=(9, 7), layout="constrained", dpi=120)
        self.ax = self.fig.add_subplot()

        self.ax.stem(self._rust.steps_stored_array, linefmt="k--", markerfmt="blue")
        self.ax.set_title(r"$Particle\ steps\ stored$")
        self.ax.set_xlabel(r"$Particle\ \#$")
        self.ax.set_ylabel(r"$Steps\ stored$")
        self.ax.set_ybound(lower=0)

        if show:
            plt.show()

    def plot_const_zeta_cartesian_poincare(
        self,
        initial: bool = False,
        color: bool = False,
        show: bool = True,
    ):
        r"""Plots the $\zeta=const$ poincare map on the $\theta-\psi$ cartesian plane.

        Parameters
        ----------
        color
            Whether or not to use a different color for each orbit. Defaults to False.
        initial
            Whether or not to plot the initial points. Note that the initial $\theta_0$ should be equal to the
            intersection angle, and the initial flux coordinate should be $\psi$. Defaults to False.
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """

        try:
            intersect_params = getattr(self, "_intersect_params")._rust
        except:
            raise RuntimeError("Queue must have run the 'intersect' routine.")

        if intersect_params.intersection != "ConstZeta":
            raise RuntimeError("Intersection surface must be 'ConstZeta'")

        if initial and self._rust[0].initial_conditions.flux0.kind == "Poloidal":
            raise ValueError("Cannot plot 'ψp0' in an 'theta-psi' plot")

        self.fig = plt.figure(**CARTESIAN_POINCARE_FIG_KW)
        self.ax = self.fig.add_subplot()

        if color:
            self.ax.set_prop_cycle(cycler(color="brcmk"))
        else:
            self.ax.set_prop_cycle(cycler(color=["blue"]))

        self.ax.set_xlabel(r"$\theta\ [rads]$")
        self.ax.set_ylabel(r"$\psi\ [Normalized]$")
        self.ax.set_xlim(-PI, PI)
        self.ax.set_xticks(
            np.linspace(-np.pi, np.pi, 5),
            [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"],
        )
        self.ax.set_title(
            rf"$\psi$ - $\theta$ cross section at $\zeta$ = {intersect_params.angle}"
        )

        for particle in self._rust.particles:
            theta = pi_mod(particle.theta_array)
            psi = particle.psi_array
            self.ax.scatter(theta, psi, **CARTESIAN_POINCARE_SCATTER_KW)
            if initial:
                init = particle.initial_conditions
                initial_point = (init.theta0, init.flux0.value)
                self.ax.plot(*initial_point, **CARTESIAN_POINCARE_INITIAL_KW)

        if show:
            plt.show()
            plt.close()

    def plot_const_theta_cartesian_poincare(
        self,
        initial: bool = False,
        color: bool = False,
        show: bool = True,
    ):
        r"""Plots the $\theta=const$ poincare map on the $\zeta-\psi_p$ cartesian plane.

        Parameters
        ----------
        color
            Whether or not to use a different color for each orbit. Defaults to False.
        initial
            Whether or not to plot the initial points. Note that the initial $\zeta_0$ should be equal to the
            intersection angle, and the initial flux coordinate should be $\psi_p$. Defaults to False.
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """

        try:
            intersect_params = getattr(self, "_intersect_params")._rust
        except:
            raise RuntimeError("Queue must have run the 'intersect' routine.")

        if intersect_params.intersection != "ConstTheta":
            raise RuntimeError("Intersection surface must be 'ConstTheta'")

        if initial and self._rust[0].initial_conditions.flux0.kind == "Toroidal":
            raise ValueError("Cannot plot 'ψ0' in an 'zeta-psip' plot")

        self.fig = plt.figure(**CARTESIAN_POINCARE_FIG_KW)
        self.ax = self.fig.add_subplot()

        if color:
            self.ax.set_prop_cycle(cycler(color="brcmk"))
        else:
            self.ax.set_prop_cycle(cycler(color=["blue"]))

        self.ax.set_xlabel(r"$\zeta\ [rads]$")
        self.ax.set_ylabel(r"$\psi_p\ [Normalized]$")
        self.ax.set_xlim(-PI, PI)
        self.ax.set_xticks(
            np.linspace(-np.pi, np.pi, 5),
            [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"],
        )
        self.ax.set_title(
            rf"$\psi_p$ - $\zeta$ cross section at $\theta$ = {intersect_params.angle}"
        )

        for particle in self._rust.particles:
            zeta = pi_mod(particle.zeta_array)
            psip = particle.psip_array
            self.ax.scatter(zeta, psip, **CARTESIAN_POINCARE_SCATTER_KW)
            if initial:
                init = particle.initial_conditions
                initial_point = (init.theta0, init.flux0.value)
                self.ax.plot(*initial_point, **CARTESIAN_POINCARE_INITIAL_KW)

        if show:
            plt.show()
            plt.close()

    def plot_const_zeta_rz_poincare(
        self,
        equilibrium: Equilibrium,
        initial: bool = False,
        color: bool = False,
        show: bool = True,
    ):
        r"""Plots the $\zeta=const$ poincare map on the $R-Z$ coordinate system.

        Parameters
        ----------
        color
            Whether or not to use a different color for each orbit. Defaults to False.
        initial
            Whether or not to plot the initial points. Note that the initial $\theta_0$ should be equal to the
            intersection angle, and the initial flux coordinate should be $\psi$. Defaults to False.
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """

        try:
            intersect_params = getattr(self, "_intersect_params")._rust
        except:
            raise RuntimeError("Queue must have run the 'intersect' routine.")

        if equilibrium.geometry is None:
            raise RuntimeError("Equilibrium's 'Geometry' must be initialized.")

        if intersect_params.intersection != "ConstZeta":
            raise RuntimeError("Intersection surface must be 'ConstZeta'")

        if initial and self._rust[0].initial_conditions.flux0.kind == "Poloidal":
            raise ValueError("Cannot plot 'ψp0' in an 'theta-psi' plot")

        self.fig = plt.figure(**CARTESIAN_POINCARE_FIG_KW)
        self.ax = self.fig.add_subplot(aspect="equal")

        if color:
            self.ax.set_prop_cycle(cycler(color="brcmk"))
        else:
            self.ax.set_prop_cycle(cycler(color=["blue"]))

        self.ax.set_xlabel(r"$R\ [m]$")
        self.ax.set_xlabel(r"$Z\ [m]$")
        self.ax.set_title(
            rf"Poincare map at $\zeta$ = {intersect_params.angle} cross section"
        )

        geom = equilibrium.geometry
        if isinstance(geom, LarGeometry) or geom.psi_state == "Good":
            rlab_of_flux = geom.rlab_of_psi
            zlab_of_flux = geom.zlab_of_psi
        else:
            rlab_of_flux = geom.rlab_of_psip
            zlab_of_flux = geom.zlab_of_psip

        for particle in self._rust.particles:
            if particle.steps_stored == 0:  # Empty arrays mess with vectorized methods
                continue
            theta = particle.theta_array % (2 * PI)
            if particle.initial_conditions.flux0.kind == "Toroidal":
                flux = particle.psi_array
            else:
                flux = particle.psip_array
            rlab = rlab_of_flux(flux, theta)
            zlab = zlab_of_flux(flux, theta)
            self.ax.plot(rlab, zlab, **RZ_POINCARE_PLOT_KW)
            if initial:
                init = particle.initial_conditions
                psi0 = init.flux0.value
                theta0 = init.theta0 % (2 * PI)
                rlab0 = rlab_of_flux(psi0, theta0)
                zlab0 = zlab_of_flux(psi0, theta0)
                initial_point = (rlab0, zlab0)
                self.ax.plot(*initial_point, **RZ_POINCARE_INITIAL_KW)

        # Cursor
        geom_center = (geom.rgeo, geom.zaxis)
        axis_point = (geom.raxis, geom.zaxis)

        def format_coord(x, y):
            r = sqrt((axis_point[0] - x) ** 2 + (axis_point[1] - y) ** 2)
            return f"(R, Z) = ({x:.5g}, {y:.5g}), r={r:.5g} "

        setattr(self.ax, "format_coord", format_coord)
        self.ax.plot(*axis_point, "ko", markersize=4, label="$R_{axis}$")
        self.ax.plot(*geom_center, "ro", markersize=4, label="$R_{geometric}$")

        rlab_last = equilibrium.geometry.rlab_last
        zlab_last = equilibrium.geometry.zlab_last
        self.ax.plot(rlab_last, zlab_last, color="k", linewidth=2)
        self.ax.legend()

        if show:
            plt.show()
            plt.close()

    def plot_qkinetic_radial_sweep(self, show: bool = True):
        r"""Plots the contained particles' calculated $q_{kinetic}$ with respect to their $\psi_0/\psi_{p,0}$.

        This is useful for low energy particles and magnetic field lines, where
        $\Delta\psi \approx $\Delta\psi_p \approx 0$.

        Parameters
        ----------
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """

        self.fig = plt.figure(**QKINETIC_POINCARE_FIG_KW)
        self.ax = self.fig.add_subplot()

        if self._rust.particles[0].initial_conditions.flux0.kind == "Toroidal":
            flux_coord = r"$\psi$"
        else:
            flux_coord = r"$\psi_p$"

        fluxes = np.asarray(
            [p.initial_conditions.flux0.value for p in self._rust.particles]
        )
        colors = np.asarray([orbit_color(p.orbit_type) for p in self._rust.particles])
        self.ax.scatter(fluxes, self._rust.qkinetic_array, c=colors, s=1)

        copassing = Patch(color=COPASSING_COLOR, label="Copassing")
        cupassing = Patch(color=CUPASSING_COLOR, label="CuPassing")
        trapped = Patch(color=TRAPPED_COLOR, label="Trapped")
        undefined = Patch(color=UNDEFINED_COLOR, label="Undefined")
        unclassified = Patch(color=UNCLASSIFIED_COLOR, label="Unclassified")

        self.ax.set_xlabel(rf"{flux_coord} $[Normalized]$")
        self.ax.set_ylabel("$q_{kinetic}$")
        self.ax.set_title("$q_{kinetic}-$" + flux_coord)
        self.ax.legend(handles=[trapped, copassing, cupassing, undefined, unclassified])
        self.ax.axhline(y=0, ls="--", lw=1.5, c="k")
        self.ax.grid(True)

        if show:
            plt.show()
            plt.close()

    def plot_qkinetic_pzeta_sweep(self, show: bool = True):
        r"""Plots the contained particles' calculated $q_{kinetic}$ with respect to their $P_\zeta$.

        Parameters
        ----------
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """

        self.fig = plt.figure(**QKINETIC_POINCARE_FIG_KW)
        self.ax = self.fig.add_subplot()

        pzetas = np.asarray([p.initial_conditions.pzeta0 for p in self._rust.particles])
        colors = np.asarray([orbit_color(p.orbit_type) for p in self._rust.particles])
        self.ax.scatter(pzetas, self._rust.qkinetic_array, c=colors, s=1)

        copassing = Patch(color=COPASSING_COLOR, label="Copassing")
        cupassing = Patch(color=CUPASSING_COLOR, label="CuPassing")
        trapped = Patch(color=TRAPPED_COLOR, label="Trapped")
        undefined = Patch(color=UNDEFINED_COLOR, label="Undefined")
        unclassified = Patch(color=UNCLASSIFIED_COLOR, label="Unclassified")

        self.ax.set_xlabel(r"$P_\zeta\ [Normalized]$")
        self.ax.set_ylabel("$q_{kinetic}$")
        self.ax.set_title(r"$q_{kinetic}-P_\zeta$")
        self.ax.legend(handles=[trapped, copassing, cupassing, undefined, unclassified])
        self.ax.axhline(y=0, ls="--", lw=1.5, c="k")
        self.ax.grid(True)

        if show:
            plt.show()
            plt.close()

    def plot_qkinetic_energy_sweep(self, show: bool = True):
        """TODO:

        Parameters
        ----------
        show
            Whether or not to call `plt.show()`. Defaults to True.
        """

        self.fig = plt.figure(**QKINETIC_POINCARE_FIG_KW)
        self.ax = self.fig.add_subplot()

        energies = np.asarray([p.initial_energy for p in self._rust.particles])
        colors = np.asarray([orbit_color(p.orbit_type) for p in self._rust.particles])
        self.ax.scatter(energies, self._rust.qkinetic_array, c=colors, s=1)

        copassing = Patch(color=COPASSING_COLOR, label="Copassing")
        cupassing = Patch(color=CUPASSING_COLOR, label="CuPassing")
        trapped = Patch(color=TRAPPED_COLOR, label="Trapped")
        undefined = Patch(color=UNDEFINED_COLOR, label="Undefined")
        unclassified = Patch(color=UNCLASSIFIED_COLOR, label="Unclassified")

        self.ax.set_xlabel(rf"$Energy\ [Normalized]$")
        self.ax.set_ylabel("$q_{kinetic}$")
        self.ax.set_title("$q_{kinetic}-Energy$")
        self.ax.legend(handles=[trapped, copassing, cupassing, undefined, unclassified])
        self.ax.axhline(y=0, ls="--", lw=1.5, c="k")
        self.ax.grid(True)

        if show:
            plt.show()
            plt.close()


# ================================================================================================


def pi_mod(arr: Array1) -> Array1:
    """Mods an angle time series in the interval [-π, π]."""
    a: Array1 = np.mod(arr, 2 * np.pi)
    a = a - 2 * np.pi * (a > np.pi)
    return a


def orbit_color(orbit_type: OrbitType) -> str:
    match orbit_type:
        case "CoPassing":
            return COPASSING_COLOR
        case "CuPassing":
            return CUPASSING_COLOR
        case "Trapped":
            return TRAPPED_COLOR
        case "Unclassified":
            return UNCLASSIFIED_COLOR
        case "Undefined":
            return UNDEFINED_COLOR

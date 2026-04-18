"""Plotter Parent classes that provide simple plotting methods"""

from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
from math import floor, log10, sqrt
from math import pi as PI
from warnings import warn
from cycler import cycler
from matplotlib.patches import Patch

from dexter._core import _PyParticle, _PyQueue
from dexter import Equilibrium, LarGeometry, OrbitType
from dexter.types import Array1, Canvas, MultiCanvas

dpi = 120
figsize = (10, 7)
EVOLUTION_TARGET_POINTS = 50_000
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


class _ParticlePlotter:
    """Provides plotting functions for the Particle Class."""

    _rust: _PyParticle

    def plot_evolution(
        self,
        percentage: float = 100,
        downsample: bool = True,
        show: bool = True,
    ) -> MultiCanvas:
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

        Returns
        -------
        MultiCanvas
            The produced `Figure` and `Axes`.
        """

        if self._rust.steps_stored == 0:
            raise Exception("Particle hasn't been integrated")

        if percentage < 0 or percentage > 100:
            raise ValueError("Percentage must be between 0 and 100.")

        step = 1
        if downsample:
            length = len(self._rust.t_array)
            oom = floor(log10(length))  # order of magnitude of number of points
            target_oom = floor(log10(EVOLUTION_TARGET_POINTS))
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

        if downsample and len(t) > EVOLUTION_TARGET_POINTS * 10:
            warn("Downsampling did not work..")

        # ===========================

        fig = plt.figure(figsize=figsize, layout="constrained", dpi=dpi)

        axes = fig.subplots(8, 1, sharex=True)
        for ax in axes:
            ax.yaxis.set_ticks_position("right")
            ax.yaxis.set_label_position("left")

        axes[0].scatter(t, psi, **SCATTER_KW)
        axes[1].scatter(t, psip, **SCATTER_KW)
        axes[2].scatter(t, theta, **SCATTER_KW)
        axes[3].scatter(t, zeta, **SCATTER_KW)
        axes[4].scatter(t, rho, **SCATTER_KW)
        axes[5].scatter(t, ptheta, **SCATTER_KW)
        axes[6].scatter(t, pzeta, **SCATTER_KW)
        axes[7].scatter(t, energy, **SCATTER_KW)

        axes[0].set_ylabel(r"$\psi$", **LABEL_KW)
        axes[1].set_ylabel(r"$\psi_p$", **LABEL_KW)
        axes[2].set_ylabel(r"$\theta$", **LABEL_KW)
        axes[3].set_ylabel(r"$\zeta$", **LABEL_KW)
        axes[4].set_ylabel(r"$\rho_{||}$", **LABEL_KW)
        axes[5].set_ylabel(r"$P_\theta$", **LABEL_KW)
        axes[6].set_ylabel(r"$P_\zeta$", **LABEL_KW)
        axes[7].set_ylabel(r"$E$", **LABEL_KW)
        axes[7].set_xlabel(r"$t[Normalized]$", **LABEL_KW)

        # Zoom out Pzeta plot
        if abs(np.std(pzeta)) < 1e-6:
            ylim = np.array(axes[6].get_ylim())
            axes[6].set_ylim(np.sort([ylim[0] / 3, ylim[1] * 3]).tolist())
        # Zoom out Energy plot
        if abs(np.std(energy)) < 1e-6:
            ylim = np.array(axes[7].get_ylim())
            axes[7].set_ylim(np.sort([ylim[0] / 2, ylim[1] * 2]).tolist())

        if show:
            plt.show()
            plt.close()

        return (fig, axes)

    def plot_db_drift(
        self,
        equilibrium: Equilibrium,
        show: bool = True,
    ) -> MultiCanvas:
        r"""Plots the magnetic field strength $B(\psi/\psi_p, \theta)$ along the particle's orbit.

        Parameters
        ----------
        equilibrium
            The equilibrium in which the particle was integrated.
        show
            Whether or not to call `plt.show()`. Defaults to True.

        Returns
        -------
        MultiCanvas
            The produced `Figure` and `Axes`.
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

        fig = plt.figure(figsize=figsize, layout="constrained", dpi=dpi)
        fig.suptitle(r"$\Delta B(\psi/\psi_p, \theta)\ drift\ during\ transit$")

        axes = fig.subplots(4, 1, sharex=True)
        for ax in axes:
            ax.yaxis.set_ticks_position("right")
            ax.yaxis.set_label_position("left")

        axes[0].scatter(t, b, **SCATTER_KW)
        axes[1].scatter(t, psi, **SCATTER_KW)
        axes[2].scatter(t, psip, **SCATTER_KW)
        axes[3].scatter(t, theta, **SCATTER_KW)

        axes[0].set_ylabel(r"$B$", **LABEL_KW)
        axes[1].set_ylabel(r"$\psi$", **LABEL_KW)
        axes[2].set_ylabel(r"$\psi_p$", **LABEL_KW)
        axes[3].set_ylabel(r"$\theta$", **LABEL_KW)
        axes[3].set_xlabel(r"$t[Normalized]$", **LABEL_KW)

        if show:
            plt.show()
            plt.close()

        return (fig, axes)


class _QueuePlotter:
    """Provides plotting functions for the Queue Class."""

    _rust: _PyQueue

    def plot_energies(self, show=True) -> Canvas:
        """Plots a stem plot of all the particles' energies.

        Particles are visited in order of instantiation.

        Useful to visualize the energy spans of the particles under study.

        Parameters
        ----------
        show
            Whether or not to call `plt.show()`. Defaults to True.

        Returns
        -------
        Canvas
            The produced `Figure` and `Ax`.
        """

        fig = plt.figure(figsize=(9, 7), layout="constrained", dpi=120)
        ax = fig.add_subplot()

        ax.stem(self._rust.energy_array, linefmt="k--", markerfmt="blue")
        ax.set_title(r"$Particle\ Energies$")
        ax.set_xlabel(r"$Particle\ \#$")
        ax.set_ylabel(r"$Energy\ [Normalized]$")
        ax.set_ybound(lower=0)

        if show:
            plt.show()
            plt.close()

        return (fig, ax)

    def plot_steps_taken(self, show=True) -> Canvas:
        """Plots a stem plot of all the number of steps each particle has taken.

        Particles are visited in order of instantiation.

        Parameters
        ----------
        show
            Whether or not to call `plt.show()`. Defaults to True.

        Returns
        -------
        Canvas
            The produced `Figure` and `Ax`.
        """

        fig = plt.figure(figsize=(9, 7), layout="constrained", dpi=120)
        ax = fig.add_subplot()

        ax.stem(self._rust.steps_taken_array, linefmt="k--", markerfmt="blue")
        ax.set_title(r"$Particle\ steps\ taken$")
        ax.set_xlabel(r"$Particle\ \#$")
        ax.set_ylabel(r"$Steps\ taken$")
        ax.set_ybound(lower=0)

        if show:
            plt.show()
            plt.close()

        return (fig, ax)

    def plot_steps_stored(self, show=True) -> Canvas:
        """Plots a stem plot of all the number of steps each particle has stored.

        Particles are visited in order of instantiation.

        Parameters
        ----------
        show
            Whether or not to call `plt.show()`. Defaults to True.

        Returns
        -------
        Canvas
            The produced `Figure` and `Ax`.
        """

        fig = plt.figure(figsize=(9, 7), layout="constrained", dpi=120)
        ax = fig.add_subplot()

        ax.stem(self._rust.steps_stored_array, linefmt="k--", markerfmt="blue")
        ax.set_title(r"$Particle\ steps\ stored$")
        ax.set_xlabel(r"$Particle\ \#$")
        ax.set_ylabel(r"$Steps\ stored$")
        ax.set_ybound(lower=0)

        if show:
            plt.show()
            plt.close()

        return (fig, ax)

    def plot_const_zeta_cartesian_poincare(
        self,
        initial: bool = False,
        color: bool = False,
        show: bool = True,
    ) -> Canvas:
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

        Returns
        -------
        Canvas
            The produced `Figure` and `Ax`.
        """

        try:
            intersect_params = getattr(self, "_intersect_params")._rust
        except:
            raise RuntimeError("Queue must have run the 'intersect' routine.")

        if intersect_params.intersection != "ConstZeta":
            raise RuntimeError("Intersection surface must be 'ConstZeta'")

        if initial and self._rust[0].initial_conditions.flux0.kind == "Poloidal":
            raise ValueError("Cannot plot 'ψp0' in an 'theta-psi' plot")

        fig = plt.figure(**CARTESIAN_POINCARE_FIG_KW)
        ax = fig.add_subplot()

        if color:
            ax.set_prop_cycle(cycler(color="brcmk"))
        else:
            ax.set_prop_cycle(cycler(color=["blue"]))

        ax.set_xlabel(r"$\theta\ [rads]$")
        ax.set_ylabel(r"$\psi\ [Normalized]$")
        ax.set_xlim(-PI, PI)
        ax.set_xticks(
            np.linspace(-np.pi, np.pi, 5),
            [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"],
        )
        ax.set_title(
            rf"$\psi$ - $\theta$ cross section at $\zeta$ = {intersect_params.angle}"
        )

        for particle in self._rust.particles:
            theta = pi_mod(particle.theta_array)
            psi = particle.psi_array
            ax.scatter(theta, psi, **CARTESIAN_POINCARE_SCATTER_KW)
            if initial:
                init = particle.initial_conditions
                initial_point = (init.theta0, init.flux0.value)
                ax.plot(*initial_point, **CARTESIAN_POINCARE_INITIAL_KW)

        if show:
            plt.show()
            plt.close()

        return (fig, ax)

    def plot_const_theta_cartesian_poincare(
        self,
        initial: bool = False,
        color: bool = False,
        show: bool = True,
    ) -> Canvas:
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

        Returns
        -------
        Canvas
            The produced `Figure` and `Ax`.
        """

        try:
            intersect_params = getattr(self, "_intersect_params")._rust
        except:
            raise RuntimeError("Queue must have run the 'intersect' routine.")

        if intersect_params.intersection != "ConstTheta":
            raise RuntimeError("Intersection surface must be 'ConstTheta'")

        if initial and self._rust[0].initial_conditions.flux0.kind == "Toroidal":
            raise ValueError("Cannot plot 'ψ0' in an 'zeta-psip' plot")

        fig = plt.figure(**CARTESIAN_POINCARE_FIG_KW)
        ax = fig.add_subplot()

        if color:
            ax.set_prop_cycle(cycler(color="brcmk"))
        else:
            ax.set_prop_cycle(cycler(color=["blue"]))

        ax.set_xlabel(r"$\zeta\ [rads]$")
        ax.set_ylabel(r"$\psi_p\ [Normalized]$")
        ax.set_xlim(-PI, PI)
        ax.set_xticks(
            np.linspace(-np.pi, np.pi, 5),
            [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"],
        )
        ax.set_title(
            rf"$\psi_p$ - $\zeta$ cross section at $\theta$ = {intersect_params.angle}"
        )

        for particle in self._rust.particles:
            zeta = pi_mod(particle.zeta_array)
            psip = particle.psip_array
            ax.scatter(zeta, psip, **CARTESIAN_POINCARE_SCATTER_KW)
            if initial:
                init = particle.initial_conditions
                initial_point = (init.theta0, init.flux0.value)
                ax.plot(*initial_point, **CARTESIAN_POINCARE_INITIAL_KW)

        if show:
            plt.show()
            plt.close()

        return (fig, ax)

    def plot_const_zeta_rz_poincare(
        self,
        equilibrium: Equilibrium,
        initial: bool = False,
        color: bool = False,
        show: bool = True,
    ) -> Canvas:
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

        Returns
        -------
        Canvas
            The produced `Figure` and `Ax`.
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

        fig = plt.figure(**CARTESIAN_POINCARE_FIG_KW)
        ax = fig.add_subplot(aspect="equal")

        if color:
            ax.set_prop_cycle(cycler(color="brcmk"))
        else:
            ax.set_prop_cycle(cycler(color=["blue"]))

        ax.set_xlabel(r"$R\ [m]$")
        ax.set_xlabel(r"$Z\ [m]$")
        ax.set_title(
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
            ax.plot(rlab, zlab, **RZ_POINCARE_PLOT_KW)
            if initial:
                init = particle.initial_conditions
                psi0 = init.flux0.value
                theta0 = init.theta0 % (2 * PI)
                rlab0 = rlab_of_flux(psi0, theta0)
                zlab0 = zlab_of_flux(psi0, theta0)
                initial_point = (rlab0, zlab0)
                ax.plot(*initial_point, **RZ_POINCARE_INITIAL_KW)

        # Cursor
        geom_center = (geom.rgeo, geom.zaxis)
        axis_point = (geom.raxis, geom.zaxis)

        def format_coord(x, y):
            r = sqrt((axis_point[0] - x) ** 2 + (axis_point[1] - y) ** 2)
            return f"(R, Z) = ({x:.5g}, {y:.5g}), r={r:.5g} "

        setattr(ax, "format_coord", format_coord)
        ax.plot(*axis_point, "ko", markersize=4, label="$R_{axis}$")
        ax.plot(*geom_center, "ro", markersize=4, label="$R_{geometric}$")

        rlab_last = equilibrium.geometry.rlab_last
        zlab_last = equilibrium.geometry.zlab_last
        ax.plot(rlab_last, zlab_last, color="k", linewidth=2)
        ax.legend()

        if show:
            plt.show()
            plt.close()

        return (fig, ax)

    def plot_qkinetic_radial_sweep(
        self,
        flux_last: float | None = None,
        show: bool = True,
    ) -> Canvas:
        r"""Plots the contained particles' calculated $q_{kinetic}$ with respect to their $\psi_0/\psi_{p,0}$.

        This is useful for low energy particles and magnetic field lines, where
        $\Delta\psi \approx $\Delta\psi_p \approx 0$.

        Parameters
        ----------
        flux_last
            The value of the magnetic flux $/psi/\psi_p$ at the Last Closed Flux Surface. If passed, the
            magnetic flux values are displayed with respect to it.
        show
            Whether or not to call `plt.show()`. Defaults to True.

        Returns
        -------
        Canvas
            The produced `Figure` and `Ax`.
        """

        fig = plt.figure(**QKINETIC_POINCARE_FIG_KW)
        ax = fig.add_subplot()

        particles = self._rust.particles
        if particles[0].initial_conditions.flux0.kind == "Toroidal":
            flux_coord = r"$\psi$"
            flux_coord_last = r"$\psi_{last}$"
        else:
            flux_coord = r"$\psi_p$"
            flux_coord_last = r"$\psi_{p,last}$"

        fluxes = np.asarray([p.initial_conditions.flux0.value for p in particles])
        fluxes[fluxes == None] = np.nan
        colors = np.asarray([orbit_color(p.orbit_type) for p in particles])

        if flux_last is None:
            ax.scatter(fluxes, self._rust.qkinetic_array, c=colors, s=1)
            ax.set_xlabel(rf"${flux_coord}\ [Normalized]$")
        else:
            ax.scatter(fluxes / flux_last, self._rust.qkinetic_array, c=colors, s=1)
            ax.set_xlabel(rf"{flux_coord}$/${flux_coord_last}")

        ax.set_ylabel("$q_{kinetic}$")
        ax.set_title("$q_{kinetic}-$" + flux_coord)
        found_orbit_types = Counter([particle.orbit_type for particle in particles])
        ax.legend(handles=orbit_color_legend_handles(found_orbit_types))
        ax.axhline(y=0, ls="--", lw=1.5, c="k")
        ax.grid(True)

        if show:
            plt.show()
            plt.close()

        return (fig, ax)

    def plot_qkinetic_pzeta_sweep(
        self,
        psip_last: float | None = None,
        show: bool = True,
    ) -> Canvas:
        r"""Plots the contained particles' calculated $q_{kinetic}$ with respect to their $P_\zeta$.

        Parameters
        ----------
        psip_last
            The value of the poloidal flux $\psi_p$ at the Last Closed Flux Surface. If passed, the
            $P_\zeta$ values are displayed with respect to $\psi_{p,last}$
        show
            Whether or not to call `plt.show()`. Defaults to True.

        Returns
        -------
        Canvas
            The produced `Figure` and `Ax`.
        """

        fig = plt.figure(**QKINETIC_POINCARE_FIG_KW)
        ax = fig.add_subplot()

        particles = self._rust.particles
        pzetas = np.asarray([p.initial_conditions.pzeta0 for p in particles])
        pzetas[pzetas == None] = np.nan
        colors = np.asarray([orbit_color(p.orbit_type) for p in particles])
        if psip_last is None:
            ax.scatter(pzetas, self._rust.qkinetic_array, c=colors, s=1)
            ax.set_xlabel(r"$P_\zeta\ [Normalized]$")
        else:
            ax.scatter(pzetas / psip_last, self._rust.qkinetic_array, c=colors, s=1)
            ax.set_xlabel(r"$P_\zeta/\psi_{p,last}$")

        ax.set_ylabel("$q_{kinetic}$")
        ax.set_title(r"$q_{kinetic}-P_\zeta$")
        found_orbit_types = Counter([particle.orbit_type for particle in particles])
        ax.legend(handles=orbit_color_legend_handles(found_orbit_types))
        ax.axhline(y=0, ls="--", lw=1.5, c="k")
        ax.grid(True)

        if show:
            plt.show()
            plt.close()

        return (fig, ax)

    def plot_qkinetic_energy_sweep(self, show: bool = True) -> Canvas:
        """TODO:

        Parameters
        ----------
        show
            Whether or not to call `plt.show()`. Defaults to True.

        Returns
        -------
        Canvas
            The produced `Figure` and `Ax`.
        """

        fig = plt.figure(**QKINETIC_POINCARE_FIG_KW)
        ax = fig.add_subplot()

        particles = self._rust.particles
        energies = np.asarray([p.initial_energy for p in particles])
        energies[energies == None] = np.nan
        colors = np.asarray([orbit_color(p.orbit_type) for p in particles])
        ax.scatter(energies, self._rust.qkinetic_array, c=colors, s=1)

        ax.set_xlabel(rf"$Energy\ [Normalized]$")
        ax.set_ylabel("$q_{kinetic}$")
        ax.set_title("$q_{kinetic}-Energy$")
        found_orbit_types = Counter([particle.orbit_type for particle in particles])
        ax.legend(handles=orbit_color_legend_handles(found_orbit_types))
        ax.axhline(y=0, ls="--", lw=1.5, c="k")
        ax.grid(True)

        if show:
            plt.show()
            plt.close()

        return (fig, ax)


# ================================================================================================


def pi_mod(arr: Array1) -> Array1:
    """Mods an angle time series in the interval [-π, π]."""
    a: Array1 = np.mod(arr, 2 * np.pi)
    a = a - 2 * np.pi * (a > np.pi)
    return a


def orbit_color(orbit_type: OrbitType) -> str:
    r"""Confined = 'bright', Lost = 'deep'"""
    match orbit_type:
        case "Undefined":
            return "xkcd:coral"
        case "TrappedConfined":
            return "xkcd:bright red"
        case "TrappedLost":
            return "xkcd:deep red"
        case "CoPassingConfined":
            return "xkcd:bright blue"
        case "CoPassingLost":
            return "xkcd:deep blue"
        case "CuPassingConfined":
            return "xkcd:bright green"
        case "CuPassingLost":
            return "xkcd:deep green"
        case "Potato":
            return "xkcd:tan"
        case "Stagnated":
            return "xkcd:sky"
        case "Unclassified":
            return "xkcd:bright purple"
        case _:  # Failed(..), also update legend handlers
            return "xkcd:indigo"


def orbit_color_legend_handles(counter: Counter) -> list[Patch]:
    res = []
    for orbit_type, counts in counter.items():
        if orbit_type.startswith("Failed"):
            continue
        res.append(
            Patch(color=orbit_color(orbit_type), label=f"{orbit_type} ({counts})")
        )
    res.append(Patch(color="xkcd:indigo", label="Failed"))
    return res

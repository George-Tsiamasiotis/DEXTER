"""Plotter Parent classes that provide simple plotting methods"""

import numpy as np
import matplotlib.pyplot as plt
from math import floor, log10
from warnings import warn

from dexter.types import Array1, IntegrationStatus

dpi = 120
figsize = (10, 7)
target_points = 50_000
SCATTER_KW = {"s": 0.8, "c": "blue"}
LABEL_KW = {"labelpad": 10, "rotation": 0, "fontsize": 15}


class _ParticlePlotter:
    """Provides plotting functions for the Particle Class."""

    integration_status: IntegrationStatus
    t_array: Array1
    psi_array: Array1
    psip_array: Array1
    theta_array: Array1
    zeta_array: Array1
    rho_array: Array1
    mu_array: Array1
    ptheta_array: Array1
    pzeta_array: Array1
    energy_array: Array1

    def plot_evolution(self, percentage: float = 100, downsample: bool = True):
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
        """

        if self.integration_status in ["Initialized", "OutOfBoundsInitialization"]:
            raise Exception("Particle hasn't been integrated")

        if percentage < 0 or percentage > 100:
            raise ValueError("Percentage must be between 0 and 100.")

        step = 1
        if downsample:
            length = len(self.t_array)
            oom = floor(log10(length))  # order of magnitude of number of points
            target_oom = floor(log10(target_points))
            if oom > target_oom:
                step = 10 ** (oom - target_oom)

        points = min(
            int(np.floor(len(self.t_array) * percentage / 100)),
            len(self.t_array),
        )

        t = self.t_array[:points][::step]
        psi = self.psi_array[:points][::step]
        psip = self.psip_array[:points][::step]
        theta = self.theta_array[:points][::step]
        zeta = self.zeta_array[:points][::step]
        rho = self.rho_array[:points][::step]
        # mu = self.mu_array[:points][::step]
        ptheta = self.ptheta_array[:points][::step]
        pzeta = self.pzeta_array[:points][::step]
        energy = self.energy_array[:points][::step]

        if downsample and len(t) > target_points * 10:
            warn("Downsampling did not work..")

        # ===========================

        fig = plt.figure(figsize=figsize, layout="constrained", dpi=dpi)

        axes: list = fig.subplots(8, 1, sharex=True)
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

        plt.show()
        plt.close()

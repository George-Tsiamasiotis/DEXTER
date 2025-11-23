"""Plotting functions for a Particle's evolution."""

from warnings import warn
from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from math import floor, log10
from dexter import Particle

s = 0.3
c = "blue"
dpi = 120
figsize = (9, 6)
labelpad = 10
target_points = 50_000


def evolution_plot(
    fig: Figure,
    particle: Particle,
    percentage: float = 100,
    downsample: bool = True,
):
    """Plots the time evolution of an integrated Particle.

    Parameters
    ----------
    fig: Figure
        The figure to plot on
    particle: Particle
        The Particle object.

    Other Parameters
    ----------------
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
    if percentage < 0 or percentage > 100:
        raise ValueError("Percentage must be between 0 and 100.")

    step = 1
    if downsample:
        length = len(particle.evolution.time)
        oom = floor(log10(length))  # order of magnitude of number of points
        target_oom = floor(log10(target_points))
        if oom > target_oom:
            step = 10 ** (oom - target_oom)

    points = min(
        int(np.floor(particle.evolution.time.shape[0] * percentage / 100)),
        len(particle.evolution.time),
    )

    time = particle.evolution.time[:points][::step]
    theta = particle.evolution.theta[:points][::step]
    psip = particle.evolution.psip[:points][::step]
    rho = particle.evolution.rho[:points][::step]
    zeta = particle.evolution.zeta[:points][::step]
    pzeta = particle.evolution.pzeta[:points][::step]
    ptheta = particle.evolution.ptheta[:points][::step]
    energy = particle.evolution.energy[:points][::step]

    if downsample and len(time) > target_points * 10:
        warn("Downsampling did not work..")

    axes: list[Axes] = fig.subplots(7, 1, sharex=True)
    for ax in axes:
        ax.yaxis.set_ticks_position("left")
        ax.yaxis.set_label_position("right")

    axes[0].scatter(time, theta, s, c)
    axes[1].scatter(time, psip, s, c)
    axes[2].scatter(time, rho, s, c)
    axes[3].scatter(time, zeta, s, c)
    axes[4].scatter(time, ptheta, s, c)
    axes[5].scatter(time, pzeta, s, c)
    axes[6].scatter(time, energy, s, c)
    # Zoom out Pzeta plot
    if abs(np.std(pzeta)) < 1e-6:
        current_ylim = np.array(axes[5].get_ylim())
        axes[5].set_ylim(np.sort([current_ylim[0] / 3, current_ylim[1] * 3]).tolist())
    # Zoom out Energy plot
    if abs(np.std(energy)) < 1e-6:
        current_ylim = np.array(axes[6].get_ylim())
        axes[6].set_ylim(np.sort([current_ylim[0] / 2, current_ylim[1] * 2]).tolist())

    axes[0].set_ylabel(r"$\theta$", rotation=0, labelpad=labelpad)
    axes[1].set_ylabel(r"$\psi_p$", rotation=0, labelpad=labelpad)
    axes[2].set_ylabel(r"$\rho_{||}$", rotation=0, labelpad=labelpad)
    axes[3].set_ylabel(r"$\zeta$", rotation=0, labelpad=labelpad)
    axes[4].set_ylabel(r"$P_\theta$", rotation=0, labelpad=labelpad)
    axes[5].set_ylabel(r"$P_\zeta$", rotation=0, labelpad=labelpad)
    axes[6].set_ylabel(r"$E$", rotation=0, labelpad=labelpad)

    plt.show()
    plt.close()

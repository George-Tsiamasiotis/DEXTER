r"""Plotting functions for simulated objects.

Functions
---------
plot_evolution
    Plots the time evolution of a particle's dynamical variables.
"""

import warnings

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from math import floor, log10

from dexter.machine.machine import Machine
from dexter.simulate.particle import Particle
from dexter._utils import _tex_unit


def plot_evolution(
    machine: Machine,
    particle: Particle,
    downsample: bool = True,
    show: bool = True,
) -> tuple[Figure, tuple[Axes]]:
    r"""Plots the time evolution of a particle's dynamical variables.

    Parameters
    ----------
    machine
        The machine in which the particle was integrated. It is used to convert specific
        quantities to SI units.
    particle
        The integrated particle.
    downsample
        Whether or not to downsample the evolution arrays. This can be
        really helpful when plotting arrays with a lot of points, since it
        drastically improves both figure creation and interaction
        performance. Downsampling is done by increasing the [::step] just
        enough so that the final number of points is more than 50.000.
    show
        Whether or not to call `plt.show()`.

    Raises
    ------
    RuntimeError
        If the particle has not been integrated (`#!python particle.steps_taken == 0`).

    Example
    -------

    ``` py title="Particle integration and evolution plotting"
    >>> machine = dex.Machine.FromNetcdf(path, "Akima", "Bicubic")
    >>>
    >>> # Initial conditions setup
    >>> initial_conditions = dex.InitialConditions.Mixed(
    ...     t0=0,
    ...     flux0=dex.MagneticFlux.Toroidal(0.1),
    ...     theta0=3.14,
    ...     zeta0=0,
    ...     pzeta0=-0.5*machine.psip_last.value,
    ...     mu0=7e-6,
    ... )
    >>>
    >>> # Particle setup and integration
    >>> particle = dex.Particle(initial_conditions)
    >>> particle.integrate(machine, (0, 500))
    >>>
    >>> fig, axes = dex.plot_evolution(machine, particle)

    ```
    """
    if particle.steps_taken == 0:
        raise RuntimeError("Particle has not been integrated")

    fig = plt.figure(figsize=(9, 5), dpi=140)
    axes = fig.subplots(4, 2, sharex=True)
    axpsi = axes[0, 0]
    axtheta = axes[1, 0]
    axptheta = axes[2, 0]
    axrho = axes[3, 0]
    axpsip = axes[0, 1]
    axzeta = axes[1, 1]
    axpzeta = axes[2, 1]
    axenergy = axes[3, 1]

    DOWNSAMPLE_TARGET = 50_000

    t_array = particle.t_array

    step = 1
    if downsample:
        length = len(t_array)
        oom = floor(log10(length))  # order of magnitude of number of points
        target_oom = floor(log10(DOWNSAMPLE_TARGET))
        if oom > target_oom:
            step = 10 ** (oom - target_oom)

    if downsample and not len(t_array) < DOWNSAMPLE_TARGET * 10:
        warnings.warn("Downsampling did not work..")

    steps_stored = particle.steps_stored
    points = min(int(np.floor(steps_stored)), steps_stored)
    initial_energy = particle.initial_energy

    t = t_array[::step]
    psi = particle.psi_array[::step] / machine.psi_last.value
    psip = particle.psip_array[::step] / machine.psip_last.value
    theta = particle.theta_array[::step]
    zeta = particle.zeta_array[::step]
    rho = particle.rho_array[::step]
    ptheta = particle.ptheta_array[::step] / machine.psi_last.value
    pzeta = particle.pzeta_array[::step] / machine.psip_last.value
    energy = particle.energy_array[::step]

    FLUX_COLOR = "xkcd:forrest green"
    ANGLE_COLOR = "xkcd:royal blue"
    RHO_COLOR = "xkcd:crimson"
    PTHETA_COLOR = "xkcd:cerulean"
    PZETA_COLOR = "xkcd:cerulean"
    ENERGY_COLOR = "xkcd:dark orange"

    PLOT_KW = dict(linewidth=0, marker=".", markersize=2)

    axpsi.plot(t, psi, c=FLUX_COLOR, **PLOT_KW)
    axpsip.plot(t, psip, c=FLUX_COLOR, **PLOT_KW)
    axtheta.plot(t, theta, c=ANGLE_COLOR, **PLOT_KW)
    axzeta.plot(t, zeta, c=ANGLE_COLOR, **PLOT_KW)
    axrho.plot(t, rho, c=RHO_COLOR, **PLOT_KW)
    axptheta.plot(t, ptheta, c=PTHETA_COLOR, **PLOT_KW)
    axpzeta.plot(t, pzeta, c=PZETA_COLOR, **PLOT_KW)
    axenergy.plot(t, energy, c=ENERGY_COLOR, **PLOT_KW)

    MARGINS = (0, 0.03)
    for ax in axes.flatten():
        ax.margins(*MARGINS)

    for ax in axes[:, 0]:
        ax.yaxis.set_ticks_position("right")
        ax.yaxis.set_label_position("left")
    for ax in axes[:, 1]:
        ax.yaxis.set_ticks_position("right")
        ax.yaxis.set_label_position("left")

    LEFT_LABEL_KW = dict(fontsize=12)
    RIGHT_LABEL_KW = dict(fontsize=12, rotation=270, labelpad=18)
    RIGHT_LABEL_KW = dict(fontsize=12)

    if machine.geometry is not None:
        # Possible divide by zero at t=0 during conversion
        with warnings.catch_warnings(action="ignore"):
            tu = machine.quantity(t, "NormSecond").to("second").to_compact()
            t = tu.m
            eu = machine.quantity(energy, "NormJoule").to("kiloelectronvolt")
            e = eu.m

        tunits = _tex_unit(tu)
        axrho.set_xlabel(rf"$t\ [{tunits}]$")
        axenergy.set_xlabel(rf"$t\ [{tunits}]$")
        axenergy.set_ylabel(rf"$E(t)\ [keV]$", **RIGHT_LABEL_KW)

    else:
        axrho.set_xlabel(rf"$t\ [Normalized]$")
        axenergy.set_xlabel(rf"$t\ [Normalized]$")
        axenergy.set_ylabel(r"$E(t)' \[Normalized]$", **RIGHT_LABEL_KW)

    axpsi.set_ylabel(r"$\psi(t)/\psi_{LCFS}$", **LEFT_LABEL_KW)
    axtheta.set_ylabel(r"$\theta(t)\ [rad]$", **LEFT_LABEL_KW)
    axptheta.set_ylabel(r"$P_\theta(t)/\psi_{LCFS}$", **LEFT_LABEL_KW)
    axrho.set_ylabel(r"$\rho_{||}(t)\ [Normalized]$", **LEFT_LABEL_KW)
    axpsip.set_ylabel(r"$\psi_p(t)/\psi_{p,LCFS}$", **RIGHT_LABEL_KW)
    axzeta.set_ylabel(r"$\zeta(t)\ [rad]$", **RIGHT_LABEL_KW)
    axpzeta.set_ylabel(r"$P_\zeta(t)/\psi_{p,LCFS}$", **RIGHT_LABEL_KW)

    if show:
        plt.show()

    return fig, axes

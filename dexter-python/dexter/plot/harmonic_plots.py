"""Plotting functions for a Harmonic's amplitude α and phase φ."""

import numpy as np
from matplotlib.axes import Axes
from dexter import Harmonic


def alpha_plot(ax: Axes, harmonic: Harmonic):
    """Plots the amplitude α extraced data, as a function of ψp.

    Parameters
    ----------
    ax: Axes
        The Axes object to plot on.
    harmonic: Harmonic
        The Harmonic object.
    """
    psip_data = harmonic.psip_data
    a_data = harmonic.a_data
    # Smooth derivative curve
    psips = np.linspace(harmonic.psip_data[0], harmonic.psip_data[-1], 1000)
    a = np.zeros((len(psips)))
    da_dpsip = np.zeros((len(psips)))
    for n in range(len(da_dpsip)):
        a[n] = harmonic.a(psips[n])
        da_dpsip[n] = harmonic.da_dpsip(psips[n])

    ax.set_xlabel(r"$\psi_p$")
    ax.set_ylabel(r"$\alpha(\psi_p)$")
    ax.set_title(r"$Amplitude$  $\alpha(\psi_p)$")
    ax.grid(True)
    ax.margins(0)

    ax.scatter(
        psip_data, a_data, c="k", s=4, zorder=2, alpha=0.8, label=r"$data$ $points$"
    )
    ax.plot(psips, a, c="b", label=r"$\alpha(\psi_p)$")
    ax.plot([], [], c="r", label=r"$\partial \alpha(\psi_p)/\partial \psi_p$")

    dax = ax.twinx()
    dax.set_ylabel(r"$\partial \alpha(\psi_p)/\partial \psi_p$")
    dax.plot(psips, da_dpsip, c="r")

    ax.legend()


def phase_plot(ax: Axes, harmonic: Harmonic):
    """Plots the phase φ extraced data, as a function of ψp.

    Parameters
    ----------
    ax: Axes
        The Axes object to plot on.
    harmonic: Harmonic
        The Harmonic object.
    """
    psip_data = harmonic.psip_data
    phase_data = harmonic.phase_data
    # Smooth derivative curve
    psips = np.linspace(harmonic.psip_data[0], harmonic.psip_data[-1], 1000)
    phase = np.zeros((len(psips)))
    for n in range(len(phase)):
        phase[n] = harmonic.phase(psips[n])

    ax.set_xlabel(r"$\psi_p$")
    ax.set_ylabel(r"$\phi(\psi_p)$")
    ax.set_title(r"$Phase$  $\phi(\psi_p)$")
    ax.grid(True)
    ax.margins(0)

    ax.scatter(
        psip_data, phase_data, c="k", s=4, zorder=2, alpha=0.8, label=r"$data$ $points$"
    )
    ax.plot(psips, phase, c="b", label=r"$\phi(\psi_p)$")

    ax.legend()

"""Plotting functions for a Currents' g(ψp), dg/dψp, I(ψp) and dI/dψp."""

import numpy as np
from matplotlib.axes import Axes
from dexter import Currents


def g_plot(ax: Axes, current: Currents):
    """Plots the toroidal current g extraced data, as a function of ψp or r.

    Parameters
    ----------
    ax: Axes
        The Axes object to plot on.
    currents: Currents
        The Currents object.
    """
    psip_data = current.psip_data
    g_data = current.g_data

    # Smooth derivative curve
    psips = np.linspace(current.psip_data[0], current.psip_data[-1], 1000)
    g = np.zeros((len(psips)))
    dg_dpsip = np.zeros((len(psips)))
    for n in range(len(dg_dpsip)):
        g[n] = current.g(psips[n])
        dg_dpsip[n] = current.dg_dpsip(psips[n])

    ax.set_xlabel(r"$\psi_p$")
    ax.set_ylabel(r"$g(\psi_p)$")
    ax.grid(True)
    ax.margins(0)

    ax.scatter(
        psip_data, g_data, c="k", s=2, zorder=2, alpha=0.8, label=r"$data$ $points$"
    )
    ax.plot(psips, g, c="b", label=r"$g(\psi_p)$")
    ax.plot([], [], c="r", label=r"$\partial g(\psi_p)/\partial \psi_p$")

    dax = ax.twinx()
    dax.set_ylabel(r"$\partial g(\psi_p)/\partial \psi_p$")
    dax.plot(psips, dg_dpsip, c="r")

    ax.legend()


def i_plot(ax: Axes, current: Currents):
    """Plots the poloidal current I extraced data, as a function of ψp or r.

    Parameters
    ----------
    ax: Axes
        The Axes object to plot on.
    currents: Currents
        The Currents object.
    """
    psip_data = current.psip_data
    i_data = current.i_data
    # Smooth derivative curve
    psips = np.linspace(current.psip_data[0], current.psip_data[-1], 1000)
    di_dpsip = np.zeros((len(psips)))
    i = np.zeros((len(psips)))
    for n in range(len(di_dpsip)):
        i[n] = current.i(psips[n])
        di_dpsip[n] = current.di_dpsip(psips[n])

    ax.set_ylabel(r"$I(\psi_p)$")
    ax.set_xlabel(r"$\psi_p$")
    ax.grid(True)
    ax.margins(0)

    ax.scatter(
        psip_data, i_data, c="k", s=2, zorder=2, alpha=0.8, label=r"$data$ $points$"
    )
    ax.plot(psips, i, c="b", label=r"$I(\psi_p)$")
    ax.plot([], [], c="r", label=r"$\partial I(\psi_p)/\partial \psi_p$")

    dax = ax.twinx()
    dax.set_ylabel(r"$\partial I(\psi_p)/\partial \psi_p$")
    dax.plot(psips, di_dpsip, c="r")
    ax.legend()

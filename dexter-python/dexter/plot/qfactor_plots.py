"""Plotting functions for a Qfactors's q(ψp), dψ/dψp, and ψ(ψp)."""

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from dexter import Qfactor

plt.rcParams["text.usetex"] = True


def q_plot(ax: Axes, qfactor: Qfactor, radial: bool = True):
    """Plots the q factor extraced and derived data, as a function of ψp or r.

    Parameters
    ----------
    ax: Axes
        The Axes object to plot on.
    qfactor: Qfactor
        The Qfactor object.

    Other Parameters
    ----------------
    radial: bool
        Whether or not to use r(ψp) as the x-axis coordinate. Defaults to True.
    """
    x = qfactor.psip_data
    q = qfactor.q_data
    q_derived = qfactor.q_data_derived

    ax.grid(True)
    ax.margins(0)

    if radial:
        q_label = r"$q(r)$"
        ax.set_xlabel(r"$r(\psi_p)$")
        ax.set_ylabel(r"$q(r)$")
        ax.set_title(r"$q(r)$")
        for i in range(len(x)):
            x[i] = qfactor.r(x[i])
    else:
        q_label = r"$q(\psi_p)$"
        ax.set_xlabel(r"$\psi_p$")
        ax.set_ylabel(r"$q(\psi_p)$")
        ax.set_title(r"$q(\psi_p)$")

    ax.scatter(x, q, c="k", s=2, zorder=2, alpha=0.8, label=r"$data$ $points$")
    ax.plot(x, q, c="b", label=q_label)
    ax.plot(x, q_derived, c="r", label=r"$d\psi / d\psi_p$")
    ax.legend()


def psi_plot(ax: Axes, qfactor: Qfactor):
    """Plots the ψ extraced data, as a function of ψp.

    Parameters
    ----------
    ax: Axes
        The Axes object to plot on.
    qfactor: Qfactor
        The Qfactor object.
    """
    psip = qfactor.psip_data
    psi = qfactor.psi_data

    ax.set_xlabel(r"$\psi_p$")
    ax.set_ylabel(r"$\psi(\psi_p)$")
    ax.set_title(r"$\psi(\psi_p)$")
    ax.grid(True)
    ax.margins(0)

    ax.scatter(psip, psi, c="k", s=2, zorder=2, alpha=0.8, label=r"$data$ $points$")
    ax.plot(psip, psi, c="b", label=r"$\psi(\psi_p)$")
    ax.legend()

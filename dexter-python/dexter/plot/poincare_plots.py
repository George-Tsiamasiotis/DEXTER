"""Plotting functions for Poincare maps."""

from cycler import cycler
from matplotlib.axes import Axes
import numpy as np
import matplotlib.pyplot as plt
from dexter import Heap, MappingParameters, Qfactor

s = 0.3
marker = "o"


def poincare_plot(
    ax: Axes,
    heap: Heap,
    params: MappingParameters,
    qfactor: Qfactor,
    radial: bool = False,
    color: bool = False,
    wall: bool = True,
):

    if not color:
        c = "blue"
    else:
        c = None  # default
        ax.set_prop_cycle(cycler(color="bgrcmyk"))

    if params.section == "ConstTheta":
        xs = heap.zetas
        ys = heap.psips
        _wall = qfactor.psip_wall
        ax.set_xlabel(r"$\zeta$")
        ax.set_ylabel(r"$\psi_p$", rotation=0)
        ax.set_title(
            rf"$\zeta-\psi_p,$ $cross$ $section$ $at$ $\theta={params.alpha:.4g}$"
        )
    else:
        xs = heap.thetas
        ys = heap.psis
        _wall = qfactor.psi_wall
        ax.set_xlabel(r"$\theta$")
        ax.set_ylabel(r"$\psi$", rotation=0)
        ax.set_title(
            rf"$\theta-\psi,$ $cross$ $section$ $at$ $\zeta={params.alpha:.4g}$"
        )

    if radial:  # ψp -> r(ψp)
        _wall = qfactor.r(_wall)
        ax.set_ylabel(r"$r(\psi_p)[m]$", rotation=90)
        for i in range(ys.shape[0]):
            for j in range(ys.shape[1]):
                ys[i, j] = (
                    qfactor.r(float(heap.psips[i, j]))
                    if not np.isnan(ys[i, j])
                    else np.nan
                )

    for i in range(len(xs)):
        ax.scatter(pi_mod(xs[i]), ys[i], s, c, marker=marker)

    ax.set_xlim(-np.pi, np.pi)
    ax.set_xticks(
        np.linspace(-np.pi, np.pi, 5),
        [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"],
    )

    if wall:
        ax.axhline(y=_wall, c="r")
    ax.set_ylim(np.clip(ax.get_ylim(), a_min=0, a_max=None).tolist())

    plt.show()
    plt.close()


def pi_mod(arr: np.ndarray):
    a = np.mod(arr, 2 * np.pi)
    a = a - 2 * np.pi * (a > np.pi)
    return a

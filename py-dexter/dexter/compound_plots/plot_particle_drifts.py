from dexter.simulate.colors import orbit_color
import numpy as np
import matplotlib.pyplot as plt

from dexter.equilibrium.equilibrium import Equilibrium
from dexter.simulate.objects import Particle
from dexter.types import Canvas

SCATTER_KW = {"s": 0.8, "c": "red"}
LOG_LOCATOR_BASE = 1 + 1e-10


def plot_particle_drifts(
    particles: list[Particle],
    equilibrium: Equilibrium,
    *,
    show: bool = True,
) -> Canvas:
    r"""Creates a contour plot of the energy, calculated on a $(\psi/\psi_p, \theta)$ grid.

    Parameters
    ----------
    particles
        The Particles to be plotted.
    equilibrium
        The equilibrium in which the particle was integrated.

    Other Parameters
    ----------------
    show
        Whether or not to call `plt.show()`. Defaults to True.

    Returns
    -------
    Canvas
        The produced `Figure` and `Ax`.
    """
    # for p in particles:
    #     if p.integration_status != "ClosedPeriods(1)":
    #         raise ValueError("Particle has not been integrated")

    fig, ax = equilibrium.geometry.plot_last(show=False)

    for p in particles:
        theta = np.mod(p.theta_array, 2 * np.pi)
        psi = p.psi_array
        r = equilibrium.geometry.rlab_of_psi(psi, theta)
        z = equilibrium.geometry.zlab_of_psi(psi, theta)
        c = orbit_color(p.orbit_type)
        ax.plot(r, z, color=c, linewidth=1)

    ax.grid(False)

    if show:
        plt.show()
        plt.close()

    return (fig, ax)

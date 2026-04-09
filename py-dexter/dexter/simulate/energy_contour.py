import numpy as np
import matplotlib.pyplot as plt
from math import pi as PI
from matplotlib.ticker import LogLocator, MaxNLocator

from dexter.equilibrium import Equilibrium, LarGeometry
from dexter.simulate import Particle, COMs
from dexter.types import Array1, Array2, Canvas, Locator

SCATTER_KW = {"s": 0.8, "c": "red"}
LOG_LOCATOR_BASE = 1 + 1e-10


def plot_energy_contour(
    flux_array: Array1,
    theta_array: Array1,
    energy_array: Array2,
    *,
    levels: int = 30,
    show: bool = True,
) -> Canvas:
    r"""Creates a contour plot of the energy, calculated on a $(\psi/\psi_p, \theta)$ grid.

    Parameters
    ----------
    theta_array
        1D array containing the $\theta$ values.
    flux_array
        1D array containing the $\psi/psi_p$ values.
    energy_array
        2D array containing the Energy values.

    Other Parameters
    ----------------
    levels
        The number of energy levels on the contour. Defaults to 30.
    show
        Whether or not to call `plt.show()`. Defaults to True.

    Returns
    -------
    Canvas
        The produced `Figure` and `Ax`.
    """
    fig = plt.figure(figsize=(10, 8), layout="constrained", dpi=120)
    ax = fig.subplots(1)

    kw = {
        "levels": levels,
        "cmap": "plasma",
    }
    contour = ax.contourf(theta_array, flux_array, energy_array, **kw)
    fig.colorbar(contour)

    ax.set_title(r"$Energy\ contour$")
    ax.set_xlabel(r"$\theta\ [rads]$")
    ax.set_ylabel(r"$flux\ [Normalized]$")

    ax.set_xticks(
        np.linspace(-np.pi, np.pi, 5),
        [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"],
    )

    if show:
        plt.show()
        plt.close()

    return (fig, ax)


def plot_particle_poloidal_drift(
    particle: Particle,
    equilibrium: Equilibrium,
    *,
    flux_span: tuple[float, float] = (0, 1),
    levels: int = 30,
    density: int = 200,
    locator: Locator = "MaxN",
    show: bool = True,
) -> Canvas:
    r"""Creates a contour plot of the energy, calculated on a $(\psi/\psi_p, \theta)$ grid.

    Parameters
    ----------
    particle
        The Particle.
    equilibrium
        The equilibrium in which the particle was integrated.

    Other Parameters
    ----------------
    flux_span
        The $\psi/\psi_p$ span in which to calculate the energy contours, with respect to
        $\psi_{wall}/\psi_{p,wall}$. Defaults to (0, 1).
    levels
        The number of energy levels on the contour. Defaults to 30.
    density
        The energy contour grid density. Defaults to 200
    show
        Whether or not to call `plt.show()`. Defaults to True.

    Returns
    -------
    Canvas
        The produced `Figure` and `Ax`.
    """
    fig, ax = equilibrium.geometry.plot_last(show=False)

    thetas = np.mod(particle.theta_array, 2 * PI)

    coms = COMs(
        pzeta=particle.initial_conditions.pzeta0,
        mu=particle.initial_conditions.mu0,
    )

    geom = equilibrium.geometry
    if isinstance(geom, LarGeometry) or geom.psi_state == "Good":
        rlab_of_flux = geom.rlab_of_psi
        zlab_of_flux = geom.zlab_of_psi
        energy_of_flux = coms.energy_of_psi_grid
        lcfs = equilibrium.psi_last
        flux = particle.psi_array
        title = r"$\theta-\psi\ drift$"
    else:
        rlab_of_flux = geom.rlab_of_psip
        zlab_of_flux = geom.zlab_of_psip
        energy_of_flux = coms.energy_of_psip_grid
        lcfs = equilibrium.psip_last
        flux = particle.psip_array
        title = r"$\theta-\psi_p\ drift$"

    rlab = rlab_of_flux(flux, thetas)
    zlab = zlab_of_flux(flux, thetas)

    flux_array = np.linspace(flux_span[0], flux_span[1], density) * lcfs
    theta_array = np.linspace(0, 2 * np.pi, density)
    psi_grid, theta_grid = np.meshgrid(flux_array, theta_array)
    energy_grid = energy_of_flux(equilibrium, flux_array, theta_array)
    rlab_grid = rlab_of_flux(psi_grid, theta_grid)
    zlab_grid = zlab_of_flux(psi_grid, theta_grid)

    _locator = locator.lower()
    _locator = (
        LogLocator(base=LOG_LOCATOR_BASE, numticks=levels)
        if _locator == "log"
        else MaxNLocator(nbins=levels)
    )

    kw = {
        "levels": levels,
        "locator": _locator,
        "cmap": "plasma",
    }

    contour = ax.contourf(rlab_grid, zlab_grid, energy_grid.T, **kw)

    fig.colorbar(contour, label=r"$Energy\ [Normalized]$")
    ax.scatter(rlab, zlab, **SCATTER_KW)
    ax.set_title(title)
    ax.grid(False)

    if show:
        plt.show()
        plt.close()

    return (fig, ax)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from dexter.types import Array1, Array2


def energy_contour(
    flux_array: Array1,
    theta_array: Array1,
    energy_array: Array2,
    levels: int = 30,
    show: bool = True,
) -> tuple[Figure, Axes]:
    r"""Creates a contour plot of the energy, calculated on a $(\psi/\psi_p, \theta$ grid.

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

    return (fig, ax)

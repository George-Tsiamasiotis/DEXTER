"""Plotting functions for Poincare maps."""

from dexter.types import NDArray2D
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib.axes import Axes
from matplotlib.patches import Circle
from dexter import Geometry, Heap, MappingParameters

s = 0.3
marker = "o"


def poincare_plot(
    ax: Axes,
    heap: Heap,
    params: MappingParameters,
    geometry: Geometry,
    radial: bool = False,
    lab: bool = False,
    color: bool = True,
    wall: bool = True,
):
    # ============== Color

    if lab and params.section == "ConstTheta":
        raise Exception("Cannot plot on a polar axis with 'ConstTheta'")

    if not color:
        c = "blue"
    else:
        c = None  # default
        ax.set_prop_cycle(cycler(color="brcmk"))

    # ============== Radial Conversion

    match (radial, params.section):
        case (False, "ConstTheta"):  # plot (ζ, ψp)
            xs = heap.zetas
            ys = heap.psips
            _wall = geometry.psip_wall
            x_coord = r"\zeta"
            y_coord = r"\psi_p"
        case (False, "ConstZeta"):  # plot (θ, ψ)
            xs = heap.thetas
            ys = heap.psis
            _wall = geometry.psi_wall
            x_coord = r"\theta"
            y_coord = r"\psi"
        case (True, "ConstTheta"):  # plot (ζ, r(ψp))
            xs = heap.zetas
            ys = np.asarray(
                [
                    [geometry.r(psip) if not np.isnan(psip) else np.nan for psip in row]
                    for row in heap.psips
                ]
            )
            _wall = geometry.r_wall
            x_coord = r"\zeta"
            y_coord = r"r(\psi_p)"
        case (True, "ConstZeta"):  # plot (θ, r(ψp))
            xs = heap.thetas
            # Use the ψs calculated from the map, since we can't convert ψ to r
            ys = np.asarray(
                [
                    [geometry.r(psip) if not np.isnan(psip) else np.nan for psip in row]
                    for row in heap.psips
                ]
            )
            _wall = geometry.r_wall
            x_coord = r"\theta"
            y_coord = r"r(\psi_p)"

    # ============== Polar projection

    if lab:  # Plot R(ψp, θ) and Z(ψp, θ), ignore everything else
        xs, ys = to_lab(geometry, heap.thetas, heap.psips)
        ax.set_aspect("equal")
        ax.set_xlabel("$R[m]$")
        ax.set_ylabel("$Z[m]$")
        ax.grid(True)
        if wall:
            wall_circle = Circle((geometry.raxis, 0), _wall, fill=False, color="r")
            ax.add_patch(wall_circle)
        ax.set_title(rf"$Cross$ $section$ $at$ $\zeta={params.alpha:.4g}$")

        # Cursor
        def format_coord(x, y):
            r = np.sqrt((geometry.raxis - x) ** 2 + (y) ** 2)
            return f"(R, Z) = ({x:.5g}, {y:.5g}), r={r:.5g} "

        ax.format_coord = format_coord  # type: ignore

    # ============== Cartesian projection

    if not lab:
        ax.set_xlim(-np.pi, np.pi)
        ax.set_xticks(
            np.linspace(-np.pi, np.pi, 5),
            [r"$-\pi$", r"$-\pi/2$", r"$0$", r"$\pi/2$", r"$\pi$"],
        )
        if wall:
            ax.axhline(y=_wall, c="r")
        ax.set_xlabel(rf"${x_coord} [rads]$")
        ax.set_ylabel(f"${y_coord} [m]$" if radial else f"${y_coord}$", rotation=0)
        ax.set_title(
            rf"${x_coord}-{y_coord},$ $cross$ $section$ $at$ ${x_coord}={params.alpha:.4g}$"
        )

    for i in range(len(xs)):
        ax.scatter(pi_mod(xs[i]), ys[i], s, c, marker=marker)

    if not lab:
        ax.set_ylim(bottom=max(0, ax.get_ylim()[0]))

    plt.show()
    plt.close()


def pi_mod(arr: np.ndarray):
    a = np.mod(arr, 2 * np.pi)
    a = a - 2 * np.pi * (a > np.pi)
    return a


def to_lab(geometry: Geometry, thetas: NDArray2D, psips: NDArray2D):

    new_xs = np.zeros(thetas.shape)
    new_ys = np.zeros(psips.shape)
    for i in range(thetas.shape[0]):
        for j in range(psips.shape[1]):
            theta = thetas[i, j]
            psip = psips[i, j]
            if np.isnan(theta) or np.isnan(psip):
                pass
                new_xs[i, j] = np.nan
                new_ys[i, j] = np.nan
            else:
                new_xs[i, j] = geometry.rlab(psip, theta)
                new_ys[i, j] = geometry.zlab(psip, theta)
    return new_xs, new_ys

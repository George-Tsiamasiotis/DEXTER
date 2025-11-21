"""Plotting functions for a Bfields's B(ψp, θ), dB/dψp, and dB/dθ."""

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from dexter import Bfield

plt.rcParams["text.usetex"] = True

levels = 30
cmap = "managua"


def b_plot(ax: Axes, bfield: Bfield):
    """Plots the b extraced data, as a function of ψp and θ.

    Parameters
    ----------
    ax: Axes
        The Axes object to plot on.
    bfield: Bfield
        The Bfield object.
    """
    r_data = bfield.rlab_data
    z_data = bfield.zlab_data
    b_data = bfield.b_data

    ax.axis("equal")
    ax.set_title(r"$Magnetic$ $field$ $strength$ $B(\psi_p,\theta)$")
    ax.set_xlabel(r"$R[m]$")
    ax.set_ylabel(r"$Z[m]$")

    contour = ax.contourf(r_data, z_data, b_data, **{"levels": levels, "cmap": cmap})
    plt.colorbar(contour, ax=ax)


def db_plot(axx: Axes, axy: Axes, bfield: Bfield):
    """Plots the db_dpsip and dp_dtheta extraced data, as a function of ψp and θ.

    Parameters
    ----------
    axx: Axes
        The Axes object to plot db_dpsip on.
    axy: Axes
        The Axes object to plot db_dtheta on.
    bfield: Bfield
        The Bfield object.
    """

    r_data = bfield.rlab_data
    z_data = bfield.zlab_data
    db_dpsip_grid = bfield.db_dpsip_data
    db_dtheta_grid = bfield.db_dtheta_data

    contour_kw = {"levels": levels, "cmap": cmap}

    axx.axis("equal")
    axy.axis("equal")
    axx.set_xlabel(r"$R[m]$")
    axx.set_ylabel(r"$Z[m]$")
    axy.set_xlabel(r"$R[m]$")
    axy.set_ylabel(r"$Z[m]$")
    axx.set_title(r"$\partial B/\partial\theta$")
    axy.set_title(r"$\partial B/\partial\psi_p$")

    contourx = axx.contourf(r_data, z_data, db_dpsip_grid, **contour_kw)
    contoury = axy.contourf(r_data, z_data, db_dtheta_grid, **contour_kw)
    plt.colorbar(contourx, ax=axx)
    plt.colorbar(contoury, ax=axy)

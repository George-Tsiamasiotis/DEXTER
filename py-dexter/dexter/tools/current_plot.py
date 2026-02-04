"""Plots the equilbrium's plasma currents as functions of ψ, ψp, or r."""

import argparse
from typing import assert_never, get_args
from dexter.types import Interp1DType

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "nc_file",
    help="the netCDF file",
)
parser.add_argument(
    "-t",
    "--typ",
    help="the interpolation type (default=steffen)",
    choices=get_args(Interp1DType),
    default="steffen",
)
parser.add_argument(
    "-c",
    "--coord",
    help="plot with respect to r[m]",
    choices=["toroidal", "poloidal", "radial"],
    default="radial",
)
parser.add_argument(
    "-d",
    dest="plot_data_points",
    help="plot data points",
    action="store_true",
)
parser.add_argument(
    "-i",
    dest="interp_points",
    help="the number of points to interpolate",
    type=int,
    default=1000,
)
args = parser.parse_args()

# ==========================================================================

import numpy as np
import matplotlib.pyplot as plt
from dexter import NcGeometry, NcCurrent


current = NcCurrent(args.nc_file, args.typ)
geometry = NcGeometry(args.nc_file, args.typ, "Bilinear")

RED = "\033[1;91m"
match args.coord:
    case "toroidal":
        if current.psi_state != "Good":
            raise SystemExit(f"{RED}The toroidal flux 'ψ' is not a good coordinate")
    case "poloidal":
        if current.psip_state != "Good":
            raise SystemExit(f"{RED}The poloidal flux 'ψp' is not a good coordinate")
    case "radial":
        pass
    case _:
        assert_never(args.coord)

print(geometry)
print(current)

fig = plt.figure(figsize=(14, 5), layout="constrained")
fig.suptitle(r"$Plasma$ $Currents$")
subplot_kw = {"xmargin": 0, "ymargin": 0}
ax = fig.subplots(1, 2, subplot_kw=subplot_kw)

# ==========================================================================

g_array = current.g_array
i_array = current.i_array

g_data_interp = np.zeros(args.interp_points)
i_data_interp = np.zeros(args.interp_points)
dg_data_interp = np.zeros(args.interp_points)
di_data_interp = np.zeros(args.interp_points)

match args.coord:
    case "toroidal":
        coord = r"\psi"
        units = r"[Normalized]"
        x_array = geometry.psi_array
        xi_array = np.linspace(0, geometry.psi_wall, args.interp_points)
        for i, psi in enumerate(xi_array):
            g_data_interp[i] = current.g_of_psi(psi)
            i_data_interp[i] = current.i_of_psi(psi)
            dg_data_interp[i] = current.dg_dpsi(psi)
            di_data_interp[i] = current.di_dpsi(psi)
    case "poloidal":
        coord = r"\psi_p"
        units = r"[Normalized]"
        x_array = geometry.psip_array
        xi_array = np.linspace(0, geometry.psip_wall, args.interp_points)
        for i, psip in enumerate(xi_array):
            g_data_interp[i] = current.g_of_psip(psip)
            i_data_interp[i] = current.i_of_psip(psip)
            dg_data_interp[i] = current.dg_dpsip(psip)
            di_data_interp[i] = current.di_dpsip(psip)
    case "radial":
        coord = r"r"
        units = r"[m]"
        x_array = geometry.r_array
        xi_array = np.linspace(0, geometry.rwall, args.interp_points)
        if geometry.psi_state == "Good":
            for i, r in enumerate(xi_array):
                psi = geometry.psi_of_r(r)
                g_data_interp[i] = current.g_of_psi(psi)
                i_data_interp[i] = current.i_of_psi(psi)
                dg_data_interp[i] = current.dg_dpsi(psi)
                di_data_interp[i] = current.di_dpsi(psi)
        if geometry.psip_state == "Good":
            for i, r in enumerate(xi_array):
                psip = geometry.psip_of_r(r)
                g_data_interp[i] = current.g_of_psip(psip)
                i_data_interp[i] = current.i_of_psip(psip)
                dg_data_interp[i] = current.dg_dpsip(psip)
                di_data_interp[i] = current.di_dpsip(psip)
    case _:
        assert_never(args.coord)


ax[0].set_xlabel(rf"${coord}{units}$")
ax[1].set_xlabel(rf"${coord}{units}$")
ax[0].set_ylabel(rf"$g({coord})$")
ax[1].set_ylabel(rf"$I({coord})$")
ax[0].set_title(rf"$g({coord})$")
ax[1].set_title(rf"$I({coord})$")
ax[0].grid(True)
ax[1].grid(True)

scatter_kw = {"c": "k", "s": 4, "zorder": 2}

if args.plot_data_points:
    ax[0].scatter(x_array, g_array, label=r"$data$ $points$", **scatter_kw)
    ax[1].scatter(x_array, i_array, label=r"$data$ $points$", **scatter_kw)
ax[0].plot(xi_array, g_data_interp, c="r", label=ax[0].get_ylabel())
ax[1].plot(xi_array, i_data_interp, c="r", label=ax[1].get_ylabel())

dax = [ax[0].twinx(), ax[1].twinx()]
dax[0].set_ylabel(rf"$\partial g({coord})/\partial {coord}$")
dax[1].set_ylabel(rf"$\partial I({coord})/\partial {coord}$")

dax[0].plot(xi_array, dg_data_interp, c="b")
dax[1].plot(xi_array, di_data_interp, c="b")
ax[0].plot([], [], c="b", label=rf"$\partial g({coord})/\partial {coord}$")
ax[1].plot([], [], c="b", label=rf"$\partial I({coord})/\partial {coord}$")
ax[0].legend()
ax[1].legend()

plt.show()
raise SystemExit

"""Plots a Currents' g(ψp), dg/dψp, I(ψp) and dI/dψp."""

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "nc_file",
    help="the netCDF file",
)
parser.add_argument(
    "-t",
    "--typ",
    help="the interpolation type (default=steffen)",
    choices=["linear", "cubic", "akima", "akimaperiodic", "steffen"],
    default="steffen",
)
parser.add_argument(
    "-r",
    dest="radial",
    help="plot with respect to r[m]",
    action="store_true",
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


currents = NcCurrent(args.nc_file, args.typ)
geometry = NcGeometry(args.nc_file, args.typ, "bilinear")
print(currents)

fig = plt.figure(figsize=(14, 5), layout="constrained")
fig.suptitle(r"$Plasma$ $Currents$")
subplot_kw = {"xmargin": 0, "ymargin": 0}
ax = fig.subplots(1, 2, subplot_kw=subplot_kw)

# ==========================================================================

g_points = currents.g_data
i_points = currents.i_data

g_data_interp = np.zeros(args.interp_points)
i_data_interp = np.zeros(args.interp_points)
dg_data_interp = np.zeros(args.interp_points)
di_data_interp = np.zeros(args.interp_points)

if args.radial:
    coord = r"r"
    units = r"[m]"
    x_data = geometry.r_data
    xi_data = np.linspace(0, geometry.r_wall, args.interp_points)
    for i, r in enumerate(xi_data):
        g_data_interp[i] = currents.g(geometry.psip(r))
        i_data_interp[i] = currents.i(geometry.psip(r))
        dg_data_interp[i] = currents.dg_dpsip(geometry.psip(r))
        di_data_interp[i] = currents.di_dpsip(geometry.psip(r))
else:
    coord = r"\psi_p"
    units = r""
    x_data = geometry.psip_data
    xi_data = np.linspace(0, geometry.psip_wall, args.interp_points)
    for i, psip in enumerate(xi_data):
        g_data_interp[i] = currents.g(psip)
        i_data_interp[i] = currents.i(psip)
        dg_data_interp[i] = currents.dg_dpsip(psip)
        di_data_interp[i] = currents.di_dpsip(psip)

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
    ax[0].scatter(x_data, g_points, label=r"$data$ $points$", **scatter_kw)
    ax[1].scatter(x_data, i_points, label=r"$data$ $points$", **scatter_kw)
ax[0].plot(xi_data, g_data_interp, c="r", label=ax[0].get_ylabel())
ax[1].plot(xi_data, i_data_interp, c="r", label=ax[1].get_ylabel())

dax = [ax[0].twinx(), ax[1].twinx()]
dax[0].set_ylabel(rf"$\partial g({coord})/\partial {coord}$")
dax[1].set_ylabel(rf"$\partial I({coord})/\partial {coord}$")

dax[0].plot(xi_data, dg_data_interp, c="b")
dax[1].plot(xi_data, di_data_interp, c="b")
ax[0].plot([], [], c="b", label=rf"$\partial g({coord})/\partial {coord}$")
ax[1].plot([], [], c="b", label=rf"$\partial I({coord})/\partial {coord}$")
ax[0].legend()
ax[1].legend()

plt.show()
raise SystemExit

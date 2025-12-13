"""Plots a Qfactors's q(ψp), dψ/dψp, and ψ(ψp)."""

import argparse
from typing import get_args
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
    default="cubic",
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
from dexter import NcGeometry, NcQfactor


qfactor = NcQfactor(args.nc_file, args.typ)
geometry = NcGeometry(args.nc_file, args.typ, "Bilinear")
print(geometry)
print(qfactor)

fig = plt.figure(figsize=(11, 5), layout="constrained")
fig.suptitle("$q-factor$ $Profile$")
subplot_kw = {"xmargin": 0, "ymargin": 0}
ax = fig.subplots(1, 2, subplot_kw=subplot_kw)

# ==========================================================================

psi_points = qfactor.psi_data
q_points = qfactor.q_data

q_data_interp = np.zeros(args.interp_points)
psi_data_interp = np.zeros(args.interp_points)
dpsi_data_interp = np.zeros(args.interp_points)

if args.radial:
    coord = r"r"
    units = r"[m]"
    x_data = geometry.r_data
    xi_data = np.linspace(0, geometry.r_wall, args.interp_points)
    for i, r in enumerate(xi_data):
        q_data_interp[i] = qfactor.q(geometry.psip(r))
        psi_data_interp[i] = qfactor.psi(geometry.psip(r))
        dpsi_data_interp[i] = qfactor.dpsi_dpsip(geometry.psip(r))
else:
    coord = r"\psi_p"
    units = r""
    x_data = geometry.psip_data
    xi_data = np.linspace(0, geometry.psip_wall, args.interp_points)
    for i, psip in enumerate(xi_data):
        q_data_interp[i] = qfactor.q(psip)
        psi_data_interp[i] = qfactor.psi(psip)
        dpsi_data_interp[i] = qfactor.dpsi_dpsip(psip)


ax[0].set_xlabel(rf"${coord}{units}$")
ax[1].set_xlabel(rf"${coord}{units}$")
ax[0].set_ylabel(rf"$q({coord})$")
ax[1].set_ylabel(rf"$\psi({coord})$")
ax[0].set_title(rf"$q({coord})$")
ax[1].set_title(rf"$\psi({coord})$")
ax[0].grid(True)
ax[1].grid(True)

scatter_kw = {"c": "k", "s": 4, "zorder": 2}

if args.plot_data_points:
    ax[0].scatter(x_data, q_points, label=r"$data$ $points$", **scatter_kw)
    ax[1].scatter(x_data, psi_points, label=r"$data$ $points$", **scatter_kw)
ax[0].plot(xi_data, q_data_interp, c="r", label=ax[0].get_ylabel())
ax[1].plot(xi_data, psi_data_interp, c="r", label=ax[1].get_ylabel())
ax[0].plot(xi_data, dpsi_data_interp, c="b", label=rf"$d\psi({coord}) / d {coord}$")
ax[0].legend()
ax[1].legend()

plt.show()
raise SystemExit

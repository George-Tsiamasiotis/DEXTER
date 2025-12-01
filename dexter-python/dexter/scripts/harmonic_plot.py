"""Plots a Harmonic's α(ψp), dα(ψp)/dψp and φ(ψp)."""

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
    "m",
    help="The 'm' mode number",
    type=int,
)
parser.add_argument(
    "n",
    help="The 'n' mode number",
    type=int,
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
from dexter import Geometry, Harmonic


harmonic = Harmonic(args.nc_file, args.typ, m=args.m, n=args.n)
geometry = Geometry(args.nc_file, args.typ, "bilinear")
print(harmonic)

fig = plt.figure(figsize=(14, 5), layout="constrained")
fig.suptitle(f"$m={args.m}, n={args.n}$ $Harmonic$")
subplot_kw = {"xmargin": 0, "ymargin": 0}
ax = fig.subplots(1, 2, subplot_kw=subplot_kw)

# ==========================================================================

a_points = harmonic.a_data
phase_points = harmonic.phase_data

a_data_interp = np.zeros(args.interp_points)
da_data_interp = np.zeros(args.interp_points)
phase_data_interp = np.zeros(args.interp_points)

if args.radial:
    coord = r"r"
    units = r"[m]"
    x_data = geometry.r_data
    xi_data = np.linspace(0, geometry.r_wall, args.interp_points)
    for i, r in enumerate(xi_data):
        a_data_interp[i] = harmonic.a(geometry.psip(r))
        da_data_interp[i] = harmonic.da_dpsip(geometry.psip(r))
        phase_data_interp[i] = harmonic.phase(geometry.psip(r))
else:
    coord = r"\psi_p"
    units = r""
    x_data = geometry.psip_data
    xi_data = np.linspace(0, geometry.psip_wall, args.interp_points)
    for i, psip in enumerate(xi_data):
        a_data_interp[i] = harmonic.a(psip)
        da_data_interp[i] = harmonic.da_dpsip(psip)
        phase_data_interp[i] = harmonic.phase(psip)

ax[0].set_xlabel(rf"${coord}{units}$")
ax[1].set_xlabel(rf"${coord}{units}$")
ax[0].set_ylabel(rf"$\alpha({coord})$")
ax[1].set_ylabel(rf"$\phi({coord})$")
ax[0].set_title(rf"$\alpha({coord})$")
ax[1].set_title(rf"$\phi({coord})$")
ax[0].grid(True)
ax[1].grid(True)

scatter_kw = {"c": "k", "s": 4, "zorder": 2}

if args.plot_data_points:
    ax[0].scatter(x_data, a_points, label=r"$data$ $points$", **scatter_kw)
    ax[1].scatter(x_data, phase_points, label=r"$data$ $points$", **scatter_kw)
ax[0].plot(xi_data, a_data_interp, c="r", label=ax[0].get_ylabel())
ax[1].plot(xi_data, phase_data_interp, c="r", label=ax[1].get_ylabel())

dax = ax[0].twinx()
dax.set_ylabel(rf"$\partial \alpha({coord})/\partial {coord}$")

dax.plot(xi_data, da_data_interp, c="b")
ax[0].plot([], [], c="b", label=rf"$\partial \alpha({coord})/\partial {coord}$")
ax[0].legend()
ax[1].legend()

plt.show()
raise SystemExit

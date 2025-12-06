"""Plots an Equilibrium's Jacobian J(ψp, θ)."""

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "nc_file",
    help="the netCDF file",
)
parser.add_argument(
    "-t",
    "--typ",
    help="the interpolation type (default=bicubic)",
    choices=["bilinear", "bicubic"],
    default="bicubic",
)
parser.add_argument(
    "-l",
    "--levels",
    help="the number of contour lines (default=20)",
    type=int,
    default=20,
)
args = parser.parse_args()

# ==========================================================================

import matplotlib.pyplot as plt
from dexter import Geometry


geometry = Geometry(args.nc_file, "linear", args.typ)
print(geometry)

fig = plt.figure(figsize=(8, 7), layout="constrained")
fig.suptitle(r"$Jacobian$ $J(\psi_p,\theta)$")
subplot_kw = {"aspect": "equal"}
ax = fig.subplots(1, subplot_kw=subplot_kw)

# ==========================================================================

r_data = geometry.rlab_data
z_data = geometry.zlab_data
jacobian_data = geometry.jacobian_data

ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")

contour_kw = {"levels": args.levels, "cmap": "plasma"}

contour = ax.contourf(r_data, z_data, jacobian_data, **contour_kw)
plt.colorbar(contour, ax=ax, cax=None)

ax.plot(r_data[-1], z_data[-1], "k-", linewidth=2)

plt.show()
raise SystemExit

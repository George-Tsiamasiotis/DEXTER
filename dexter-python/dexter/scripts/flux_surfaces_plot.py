"""Plots the equilibrium's poloidal flux surfaces."""

import argparse
from math import sqrt

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
    "-n",
    "--number",
    help="The number of flux surfaces to (try to) plot (default=60)",
    type=int,
    default=20,
)
args = parser.parse_args()

# ==========================================================================

import matplotlib.pyplot as plt
from dexter import NcBfield, NcGeometry


bfield = NcBfield(args.nc_file, args.typ)
geometry = NcGeometry(args.nc_file, "linear", args.typ)
print(geometry)
print(bfield)

fig = plt.figure(figsize=(8, 7), layout="constrained")
fig.suptitle(r"$Poloidal$ $flux$ $surfaces$")
subplot_kw = {"aspect": "equal"}
ax = fig.subplots(1, subplot_kw=subplot_kw)

# ==========================================================================

r_data = geometry.rlab_data
z_data = geometry.zlab_data

ax.set_xlabel(r"$R[m]$")
ax.set_ylabel(r"$Z[m]$")
ax.grid(True)

step = max([1, int(geometry.shape[0] / args.number)])
print(f"Displaying {int(geometry.shape[0]/step)} surfaces.")
for i in range(0, geometry.shape[0], step):
    ax.plot(r_data[i], z_data[i], color="blue", zorder=-1)

geom_center = (geometry.rgeo, geometry.zaxis)
axis_point = (geometry.raxis, geometry.zaxis)
ax.plot(r_data[-1], z_data[-1], "r-", linewidth=2, label="wall")


# Cursor
def format_coord(x, y):
    r = sqrt((axis_point[0] - x) ** 2 + (axis_point[1] - y) ** 2)
    return f"(R, Z) = ({x:.5g}, {y:.5g}), r={r:.5g} "


ax.format_coord = format_coord  # type: ignore

ax.plot(*axis_point, "ko", markersize=4, label="$R_{axis}$")
ax.plot(*geom_center, "ro", markersize=4, label="$R_{geometric}$")
ax.legend()

plt.show()
raise SystemExit

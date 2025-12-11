"""Plots a Bfields's B(ψp, θ), dB/dψp, and dB/dθ."""

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
from dexter import NcBfield, NcGeometry


bfield = NcBfield(args.nc_file, args.typ)
geometry = NcGeometry(args.nc_file, "linear", args.typ)
print(geometry)
print(bfield)

fig = plt.figure(figsize=(15, 6), layout="constrained")
fig.suptitle(r"$Magnetic$ $Field$ $Profile$")
subplot_kw = {"aspect": "equal"}
ax = fig.subplots(1, 3, subplot_kw=subplot_kw)

# ==========================================================================

r_data = geometry.rlab_data
z_data = geometry.zlab_data
b_data = bfield.b_data
db_dpsip_data = bfield.db_dpsip_data
db_dtheta_data = bfield.db_dtheta_data

ax[0].set_title(r"$Magnetic$ $field$ $strength$ $B(\psi_p,\theta)$")
ax[1].set_title(r"$\partial B(\psi_p,\theta)/\partial\theta$")
ax[2].set_title(r"$\partial B(\psi_p,\theta)/\partial\psi_p$")
ax[0].set_xlabel(r"$R[m]$")
ax[1].set_xlabel(r"$R[m]$")
ax[2].set_xlabel(r"$R[m]$")
ax[0].set_ylabel(r"$Z[m]$")
ax[1].set_ylabel(r"$Z[m]$")
ax[2].set_ylabel(r"$Z[m]$")

contour_kw = {"levels": args.levels, "cmap": "gist_heat"}
colorbar_kw = {"location": "bottom"}

contour1 = ax[0].contourf(r_data, z_data, b_data, **contour_kw)
contour2 = ax[1].contourf(r_data, z_data, db_dtheta_data, **contour_kw)
contour3 = ax[2].contourf(r_data, z_data, db_dpsip_data, **contour_kw)
plt.colorbar(contour1, ax=ax[0], cax=None, **colorbar_kw)
plt.colorbar(contour2, ax=ax[1], cax=None, **colorbar_kw)
plt.colorbar(contour3, ax=ax[2], cax=None, **colorbar_kw)

for ax in ax:
    ax.plot(r_data[-1], z_data[-1], "k-", linewidth=2)

plt.show()
raise SystemExit

"""Plots the flux surfaces on a numerical equilibrium."""

import argparse
from dexter import numerical_equilibrium

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "nc_file",
    help="the netCDF file",
)
parser.add_argument(
    "-n",
    "--number",
    help="The number of flux surfaces to (try to) plot (default=20)",
    type=int,
    default=20,
)
args = parser.parse_args()

equilibrium = numerical_equilibrium(args.nc_file, "Cubic", "Bicubic")

equilibrium.plot_flux_surfaces(args.number)

raise SystemExit

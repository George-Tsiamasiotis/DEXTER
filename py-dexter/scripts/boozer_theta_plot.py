"""Plots the `boozer_theta = const` lines on a numerical equilibrium."""

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
    help="The number of flux surfaces to (try to) plot (default=60)",
    type=int,
    default=60,
)
args = parser.parse_args()

equilibrium = numerical_equilibrium(args.nc_file, "Cubic", "Bicubic")

equilibrium.plot_boozer_theta(args.number)

raise SystemExit

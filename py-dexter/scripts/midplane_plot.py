"""Plots the flux surfaces on a numerical equilibrium."""

import argparse
from dexter import numerical_equilibrium

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "nc_file",
    help="the netCDF file",
)
args = parser.parse_args()

equilibrium = numerical_equilibrium(args.nc_file, "Cubic", "Bicubic")

equilibrium.plot_midplane()

raise SystemExit

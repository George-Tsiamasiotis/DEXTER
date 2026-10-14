"""Plots a BfieldObject's `B` and its derivatives on the R-Z plane."""

import argparse
from dexter._utils import get_default_args
from dexter import Machine, Interpolation2dType, plot_bfield

plot_bfield_args = get_default_args(plot_bfield)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "nc_file",
    help="the netCDF file",
    type=str,
)
parser.add_argument(
    "-i",
    "--interpolation-type",
    help="The kind 2D Interpolation. Defaults to 'Bicubic'.",
    choices=Interpolation2dType.__args__,
    type=str,
    default="Bicubic",
)
parser.add_argument(
    "-l",
    "--levels",
    help=f"The number of contour levels. Defaults to {plot_bfield_args["levels"]}.",
    type=int,
    default=plot_bfield_args["levels"],
)
parser.add_argument(
    "-c",
    "--contour-points",
    help=f"Whether or not to plot the points comprising the stationary curve. Defaults to {plot_bfield_args["contour_points"]}.",
    action="store_true",
)
args = parser.parse_args()

machine = Machine.FromNetcdf(args.nc_file, "Linear", args.interpolation_type)

plot_bfield(machine, contour_points=args.contour_points, show=True)
print(machine.bfield)

raise SystemExit

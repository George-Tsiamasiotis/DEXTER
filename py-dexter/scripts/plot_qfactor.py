"""Plots a QfactorObject q(ψ), q(ψp), ψp(ψ)$ and ψ(ψp)."""

import argparse
from dexter._utils import _get_default_args
from dexter import Machine, Interpolation1dType, plot_qfactor

plot_qfactor_args = _get_default_args(plot_qfactor)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "nc_file",
    help="the netCDF file",
    type=str,
)
parser.add_argument(
    "-i",
    "--interpolation-type",
    help="The kind 1D Interpolation. Defaults to 'Akima'.",
    choices=Interpolation1dType.__args__,
    type=str,
    default="Cubic",
)
parser.add_argument(
    "-p",
    "--points",
    help="The number of flux points to evaluate. "
    f"Defaults to {plot_qfactor_args["points"]}.",
    type=int,
    default=plot_qfactor_args["points"],
)
parser.add_argument(
    "-d",
    help="Whether or not to plot the data points.",
    action="store_true",
)
args = parser.parse_args()

machine = Machine.FromNetcdf(args.nc_file, args.interpolation_type, "Bilinear")

plot_qfactor(machine, points=args.points, data=args.d, show=True)

raise SystemExit

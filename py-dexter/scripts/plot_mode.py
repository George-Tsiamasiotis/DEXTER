"""Plots a ModeObject's `α`, `φ` and its derivatives with respect to `ψ` or `ψp`."""

import argparse
from dexter._utils import _get_default_args
from dexter import (
    Machine,
    Perturbation,
    NcFluteMode,
    PhaseMethod,
    MagneticFluxKind,
    Interpolation1dType,
    plot_mode,
)

plot_mode_args = _get_default_args(plot_mode)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "nc_file",
    help="the netCDF file",
    type=str,
)
parser.add_argument(
    "m",
    help="the poloidal mode number `m`.",
    type=int,
)
parser.add_argument(
    "n",
    help="the toroidal mode number `n`.",
    type=int,
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
    "-f",
    "--flux",
    help="The kind of magnetic flux on the x-axis. "
    f"Defaults to {plot_mode_args["flux"]}.",
    choices=MagneticFluxKind.__args__,
    type=str,
    default=plot_mode_args["flux"],
)
parser.add_argument(
    "-p",
    "--points",
    help="The number of flux points to evaluate. "
    f"Defaults to {plot_mode_args["points"]}.",
    type=int,
    default=plot_mode_args["points"],
)
parser.add_argument(
    "-d",
    help="Whether or not to plot the data points.",
    action="store_true",
)
args = parser.parse_args()

mode = NcFluteMode(
    path=args.nc_file,
    interp_type=args.interpolation_type,
    m=args.m,
    n=args.n,
    phase_method="Interpolation",
)

plot_mode(mode, flux=args.flux, points=args.points, data=args.d, show=True)

raise SystemExit

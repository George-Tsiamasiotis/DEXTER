"""Creates contour plots for an NcEquilibrium's magnetic field strength and derivatives."""

import argparse
from dexter import numerical_equilibrium
from dexter.types import UnitSystem

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "nc_file",
    help="the netCDF file",
)
parser.add_argument(
    "-l",
    "--levels",
    help="the number of contour lines (default=20)",
    type=int,
    default=20,
)
parser.add_argument(
    "-u",
    "--units",
    help="The units of the magnetic field strength (default='SI')",
    type=str,
    choices=UnitSystem.__args__,
    default="SI",
)
args = parser.parse_args()

equilibrium = numerical_equilibrium(args.nc_file, "Cubic", "Bicubic")

equilibrium.plot_b(args.levels, show=False)
equilibrium.plot_db(args.levels)

raise SystemExit

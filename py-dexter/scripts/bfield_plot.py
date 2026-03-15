"""Creates contour plots for an NcEquilibrium's magnetic field strength and derivatives."""

import argparse
from dexter import NcEquilibrium

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
args = parser.parse_args()

equilibrium = NcEquilibrium(args.nc_file, "Cubic", "Bicubic")

equilibrium.plot_b(args.levels)
equilibrium.plot_db(args.levels)

raise SystemExit

# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "matplotlib",
#     "dexter",
# ]
# ///
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
    help="the interpolation type (default=cubic)",
    choices=["bilinear", "bicubic"],
    default="bicubic",
)
args = parser.parse_args()

# ==========================================================================


import matplotlib.pyplot as plt
from dexter import Bfield
from dexter.plot import b_plot, db_plot


bfield = Bfield(args.nc_file, args.typ)

fig = plt.figure(figsize=(15, 5), layout="constrained")
fig.suptitle(r"$Magnetic$ $Field$ $Profile$")

ax = fig.subplots(1, 3)
b_plot(ax[0], bfield)
db_plot(ax[1], ax[2], bfield)

plt.show()

raise SystemExit

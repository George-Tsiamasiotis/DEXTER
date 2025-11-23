# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "matplotlib",
#     "dexter",
# ]
# ///
"""Plots a Currents' g(ψp), dg/dψp, I(ψp) and dI/dψp."""

import argparse

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "nc_file",
    help="the netCDF file",
)
parser.add_argument(
    "-t",
    "--typ",
    help="the interpolation type (default=steffen)",
    choices=["linear", "cubic", "akima", "akimaperiodic", "steffen"],
    default="steffen",
)
args = parser.parse_args()

# ==========================================================================


import matplotlib.pyplot as plt
from dexter import Currents
from dexter.plot import g_plot, i_plot


currents = Currents(args.nc_file, args.typ)

fig = plt.figure(figsize=(15, 5), layout="constrained")
fig.suptitle(r"$Plasma$ $Currents$")

ax = fig.subplots(1, 2)
g_plot(ax[0], currents)
i_plot(ax[1], currents)

plt.show()

raise SystemExit

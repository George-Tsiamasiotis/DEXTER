# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "matplotlib",
#     "dexter",
# ]
# ///
"""Plots a Qfactors's q(ψp), dψ/dψp, and ψ(ψp)."""

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
from dexter import Qfactor
from dexter.plot import q_plot, psi_plot


qfactor = Qfactor(args.nc_file, args.typ)

fig = plt.figure(figsize=(11, 5), layout="constrained")
fig.suptitle("$q-factor$ $Profile$")

ax = fig.subplots(1, 2)
q_plot(ax[0], qfactor)
psi_plot(ax[1], qfactor)

plt.show()

raise SystemExit

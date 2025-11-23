# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "matplotlib",
#     "dexter",
# ]
# ///
"""Plots a Harmonic's α(ψp), dα(ψp)/dψp and φ(ψp)."""

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
parser.add_argument(
    "m",
    help="The 'm' mode number",
    type=int,
)
parser.add_argument(
    "n",
    help="The 'n' mode number",
    type=int,
)
args = parser.parse_args()

# ==========================================================================


import matplotlib.pyplot as plt
from dexter import Harmonic
from dexter.plot import alpha_plot, phase_plot


harmonic = Harmonic(args.nc_file, args.typ, args.m, args.n)

fig = plt.figure(figsize=(15, 5), layout="constrained")
fig.suptitle(f"$m={args.m}$, $n={args.n}$ $Harmonic$ $({args.typ})$")

ax = fig.subplots(1, 2)
alpha_plot(ax[0], harmonic)
phase_plot(ax[1], harmonic)

plt.show()

raise SystemExit

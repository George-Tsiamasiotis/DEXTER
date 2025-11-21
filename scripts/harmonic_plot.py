#!/usr/bin/env -S uv run --script --active
# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "matplotlib",
#     "dexter",
# ]
# ///
import sys
import matplotlib
import matplotlib.pyplot as plt
from dexter import Harmonic
from dexter.plot import alpha_plot, phase_plot


matplotlib.use("gtk3agg")

(m, n) = (0, 1)  # default
typ = "steffen"
if len(sys.argv) > 2:  # like C
    m = int(sys.argv[1])
    n = int(sys.argv[2])
if len(sys.argv) > 3:
    typ = str(sys.argv[3])


harmonic = Harmonic("./data.nc", typ=typ, m=m, n=n)

fig = plt.figure(figsize=(15, 5), layout="constrained")
fig.suptitle(f"$m={m}$, $n={n}$ $Harmonic$")

ax = fig.subplots(1, 2)
alpha_plot(ax[0], harmonic)
phase_plot(ax[1], harmonic)

plt.show()

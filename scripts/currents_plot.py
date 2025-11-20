#!/usr/bin/env -S uv run --script --active
# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "matplotlib",
#     "pyncare",
# ]
# ///
import sys
import matplotlib
import matplotlib.pyplot as plt
from pyncare import Currents, g_plot, i_plot


matplotlib.use("gtk3agg")

typ = "cubic"  # default
if len(sys.argv) == 2:  # like C
    typ = str(sys.argv[1])

currents = Currents("./data.nc", typ)

fig = plt.figure(figsize=(15, 5), layout="constrained")
fig.suptitle("Plasma Currents")

ax = fig.subplots(1, 2)
g_plot(ax[0], currents)
i_plot(ax[1], currents)

plt.show()

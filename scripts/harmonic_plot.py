#!/usr/bin/env -S uv run --script --active
# /// script
# requires-python = ">=3.14"
# dependencies = [
#     "matplotlib",
#     "pyncare",
# ]
# ///
import matplotlib
import matplotlib.pyplot as plt
from pyncare import Harmonic, alpha_plot, phase_plot


matplotlib.use("gtk3agg")

harmonic = Harmonic("./data.nc", "cubic", m=0, n=1)

fig = plt.figure(figsize=(15, 5), layout="constrained")
fig.suptitle("Plasma Currents")

ax = fig.subplots(1, 2)
alpha_plot(ax[0], harmonic)
phase_plot(ax[1], harmonic)

plt.show()
